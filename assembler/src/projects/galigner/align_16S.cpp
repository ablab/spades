//***************************************************************************
//* Copyright (c) 2015-2018 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <iostream>
#include <fstream>

#include "llvm/Support/YAMLParser.h"
#include "llvm/Support/YAMLTraits.h"
#include <cxxopts/cxxopts.hpp>

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
#include "pipeline/config_struct.hpp"

#include "io/reads/io_helper.hpp"
#include "common/assembly_graph/core/coverage.hpp"

#include "utils/standard_base.hpp"
#include "utils/stl_utils.hpp"
#include "utils/logger/log_writers.hpp"

#include "modules/alignment/pacbio/gap_dijkstra.hpp"
#include "modules/alignment/long_read_mapper.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "io/reads/multifile_reader.hpp"

#include "modules/alignment/sequence_mapper.hpp"
#include "edlib/edlib.h"

#include "primer_aligner.hpp"

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace llvm {
namespace yaml {

template<> struct MappingTraits<graph_aligner::GapClosingConfig> {
    static void mapping(IO& io, graph_aligner::GapClosingConfig& cfg) {
        io.mapRequired("max_vertex_in_gap", cfg.max_vertex_in_gap);
        io.mapRequired("queue_limit", cfg.queue_limit);
        io.mapRequired("iteration_limit", cfg.iteration_limit);
        io.mapRequired("find_shortest_path", cfg.find_shortest_path);
        io.mapRequired("restore_mapping", cfg.restore_mapping);
        io.mapRequired("penalty_interval", cfg.penalty_interval);
    }
};


template<> struct MappingTraits<rna16S_mapping::RnaAlignerConfig> {
    static void mapping(IO& io, rna16S_mapping::RnaAlignerConfig& cfg) {
        io.mapRequired("path_to_sequences", cfg.path_to_sequences);
        io.mapRequired("path_to_primers", cfg.path_to_primers);

        io.mapRequired("path_limit_stretching", cfg.path_limit_stretching);
        io.mapRequired("gap_min_ed", cfg.gap_min_ed);
        io.mapRequired("res_min_ed", cfg.res_min_ed);
        io.mapRequired("position_diff", cfg.position_diff);
        io.mapRequired("alignment_min_length", cfg.alignment_min_length);
        io.mapRequired("primer_max_ed", cfg.primer_max_ed);
        io.mapRequired("graph_max_ed", cfg.graph_max_ed);

        io.mapRequired("gap_closing", cfg.gap_cfg);
    }
};

}
}

namespace rna16S_mapping {
struct Mapping {

    Mapping(const EdgeId &e, const MappingRange &range): last_edge_(e), last_mapping_(range), score_(std::numeric_limits<int>::max()) {}

    Mapping(const EdgeId &e, const MappingRange &range, int score): last_edge_(e), last_mapping_(range), score_(score) {}

    Mapping(const EdgeId &e, const MappingRange &range, const omnigraph::MappingPath<debruijn_graph::EdgeId> &path, int score)
        : last_edge_(e), last_mapping_(range), path_(path), score_(score) {}

    Mapping(): score_(-1) {}

    omnigraph::MappingPath<debruijn_graph::EdgeId> path_;
    EdgeId last_edge_;
    MappingRange last_mapping_;
    int score_;
};

class SequenceAligner {
private:
    const debruijn_graph::conj_graph_pack &gp_;
    const rna16S_mapping::RnaAlignerConfig cfg_;
    const string &output_file_;

public:
    SequenceAligner(const debruijn_graph::conj_graph_pack &gp,
                    const rna16S_mapping::RnaAlignerConfig &cfg,
                    const string &output_file):
        gp_(gp), cfg_(cfg), output_file_(output_file) {}


    bool IsCanonical(EdgeId e) const {
        return e <= gp_.g.conjugate(e);
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : gp_.g.conjugate(e);
    }

    void PrepareInitialState(omnigraph::MappingPath<debruijn_graph::EdgeId> &path, const string &s, bool forward, string &ss, EdgeId &start_e, int &start_pos, int &seq_start_pos) const {
        if (forward) {
            start_e = path.edge_at(path.size() - 1);
            omnigraph::MappingRange mapping = path.mapping_at(path.size() - 1);
            start_pos = mapping.mapped_range.start_pos;
            ss = s.substr(mapping.initial_range.start_pos, (int) s.size() - mapping.initial_range.start_pos );
            seq_start_pos = mapping.initial_range.start_pos;
            DEBUG("sz=" << path.size() << " start_pos=" << start_pos << " s_sz=" << ss.size());
        } else {
            start_e = gp_.g.conjugate(path.edge_at(0));
            omnigraph::MappingRange mapping = path.mapping_at(0);
            if (gp_.g.length(start_e) + gp_.g.k() - mapping.mapped_range.start_pos < gp_.g.length(start_e)) {
                start_pos = gp_.g.length(start_e) + gp_.g.k() - mapping.mapped_range.start_pos;
                seq_start_pos = mapping.initial_range.start_pos;
            } else {
                start_pos = gp_.g.length(start_e) - 1;
                seq_start_pos = mapping.initial_range.start_pos + (gp_.g.k() - mapping.mapped_range.start_pos + 1);
            }
            string c_ss = s.substr(0, seq_start_pos);

            map<char, char> nucs = {{'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}, {'U', 'A'}
                , {'R', 'Y'}, {'Y', 'R'}
                , {'K', 'M'}, {'M', 'K'}
                , {'S', 'W'}, {'W', 'S'}
                , {'B', 'V'}, {'V', 'B'}
                , {'D', 'H'}, {'H', 'D'}
                , {'N', 'N'}
            };
            ss = "";
            int num = 0;
            for (int i = c_ss.size() - 1; i >= 0; -- i) {
                ss = ss + nucs[c_ss[i]];
                num ++;
            }

            DEBUG("c_ss=" << c_ss.size() << " ss=" << ss.size() << " num=" << num);
        }
    }

    void UpdatePath(omnigraph::MappingPath<debruijn_graph::EdgeId> &path,
                    std::vector<EdgeId> &ans,
                    int start_pos, int end_pos,
                    int seq_start_pos, int seq_end_pos,
                    bool forward) const {
        if (forward) {
            omnigraph::MappingPath<debruijn_graph::EdgeId> cur_sorted;
            for (int i = 0; i < path.size() - 1; ++i) {
                cur_sorted.push_back(path[i].first, path[i].second);
            }
            if (ans.size() ==  1) {
                cur_sorted.push_back(ans[0], omnigraph::MappingRange(Range(seq_start_pos, seq_start_pos + seq_end_pos),
                                     Range(start_pos, end_pos) ));
            } else {
                cur_sorted.push_back(ans[0], omnigraph::MappingRange(Range(seq_start_pos, seq_start_pos + gp_.g.length(ans[0]) ),
                                     Range(0, gp_.g.length(ans[0])) ));
                for (int i = 1; i < ans.size() - 1; ++i) {
                    cur_sorted.push_back(ans[i], omnigraph::MappingRange(Range(0, gp_.g.length(ans[i])), Range(0, gp_.g.length(ans[i]) ) ));
                }
                DEBUG("forward " << seq_end_pos << " endps=" << end_pos << " anssz=" << ans.size() )
                cur_sorted.push_back(ans[ans.size() - 1], omnigraph::MappingRange(Range(0, seq_start_pos + seq_end_pos), Range(0, end_pos ) )); //0, seq_start_pos + seq_end_pos -  gp_.g.k()
            }
            path = cur_sorted;
        } else {
            omnigraph::MappingPath<debruijn_graph::EdgeId> cur_sorted;
            int start = gp_.g.length(ans[ans.size() - 1]) + gp_.g.k() - end_pos;
            int seq_start = seq_start_pos - seq_end_pos;
            int cur_ind = ans.size() - 1;
            while (cur_ind >= 0 && start - (int) gp_.g.length(ans[cur_ind]) > 0) {
                start -= gp_.g.length(ans[cur_ind]);
                cur_ind --;
            }
            if (cur_ind == 0) {
                cur_sorted.push_back(gp_.g.conjugate(ans[cur_ind]), omnigraph::MappingRange(Range(seq_start, seq_start + path[0].second.initial_range.end_pos), Range(start, start + seq_end_pos ) ));
            } else {
                cur_sorted.push_back(gp_.g.conjugate(ans[cur_ind]), omnigraph::MappingRange(Range(seq_start, seq_start + gp_.g.length(ans[cur_ind])), Range(start, gp_.g.length(ans[cur_ind]) ) ));
                for (int i = cur_ind - 1; i > 0; --i) {
                    cur_sorted.push_back(gp_.g.conjugate(ans[i]), omnigraph::MappingRange(Range(0, gp_.g.length(ans[0])), Range(0, gp_.g.length(ans[i])) ));
                }
                cur_sorted.push_back(gp_.g.conjugate(ans[0]), omnigraph::MappingRange(Range(0 , path[0].second.initial_range.end_pos), Range(0, path[0].second.mapped_range.end_pos) ));
            }
            for (int i = 1; i < path.size(); ++i) {
                cur_sorted.push_back(path[i].first, path[i].second);
            }

            path = cur_sorted;
        }
    }

    int GrowEnds(omnigraph::MappingPath<debruijn_graph::EdgeId> &path, const string &s, bool forward) const {
        VERIFY(path.size() > 0);
        string ss;
        int start_pos = -1;
        int seq_start_pos = -1;
        EdgeId start_e = EdgeId();
        PrepareInitialState(path, s, forward, ss, start_e, start_pos, seq_start_pos);

        int s_len = int(ss.size());
        int score = cfg_.gap_min_ed;
        if (s_len > 2000) {
            DEBUG("EdgeDijkstra: sequence is too long " << s_len)
            return 0;
        }
        if (s_len < (int) gp_.g.k()) {
            DEBUG("EdgeDijkstra: sequence is too small " << s_len)
            return 0;
        }
        DEBUG(" EdgeDijkstra: String length " << s_len << " max-len " << score << " start_pos=" << start_pos << " start_e=" << start_e.int_id() << " e_len=" << gp_.g.length(start_e) + gp_.g.k() );
        graph_aligner::DijkstraEndsReconstructor algo(gp_.g, cfg_.gap_cfg, ss, start_e, start_pos, score, true);
        algo.CloseGap();
        score = algo.edit_distance();
        if (score == -1) {
            DEBUG("EdgeDijkstra didn't find anything edge=" << start_e.int_id() << " s_start=" << start_pos << " seq_len=" << ss.size())
            return 0;
        } else {
            DEBUG("EdgeDijkstra found path edge=" << start_e.int_id() << " s_start=" << start_pos << " seq_len=" << ss.size() << " score=" << score )
        }

        std::vector<EdgeId> ans = algo.path();
        int end_pos = algo.path_end_position();
        int seq_end_pos = algo.seq_end_position();
        DEBUG("graph_end=" << end_pos << " seq_end_pos=" << seq_end_pos << " anssz=" << ans.size());
        UpdatePath(path, ans, start_pos, end_pos, seq_start_pos, seq_end_pos, forward);
        return score;
    }

    bool Consistent(const EdgeId &a_e, const MappingRange &a_range,
                    const EdgeId &b_e, const MappingRange &b_range,
                    std::map<VertexId, size_t> &vertex_pathlen) {
        int start_pos = a_range.initial_range.start_pos;
        int end_pos = b_range.initial_range.start_pos;
        int seq_len = -start_pos + end_pos;
        int path_max_length = (int) ((double) seq_len * cfg_.path_limit_stretching);

        //INFO("path_max_length=" << path_max_length << " seq_len=" << seq_len);
        VertexId start_v = gp_.g.EdgeEnd(a_e);
        VertexId end_v = gp_.g.EdgeStart(b_e);
        if (end_pos < start_pos) {
            //INFO("modifying limits because of some bullshit magic, seq length 0 " << start_pos << " " << end_pos)
            end_pos = start_pos;
            return false;
        }
        if (a_e == b_e) {
            if (b_range.mapped_range.start_pos < a_range.mapped_range.start_pos) {
                //if (vertex_pathlen.count(end_v) == 0){
                VERIFY(dist_.count(start_v) > 0)
                if (dist_[start_v].count(end_v) == 0) {
                    DEBUG(" Equaledges, but end before start")
                    return false;
                } else {
                    DEBUG(" Equaledges, but end before start, but ok")
                    return true;
                }
            }
            DEBUG(" Equaledges")
            return true;
        }

        //if (vertex_pathlen.count(end_v) == 0) {
        VERIFY(dist_.count(start_v) > 0)
        if (dist_[start_v].count(end_v) == 0) {
            //INFO(" Cant get")
            return false;
        }
        return true;
    }

    Mapping FillGap(const EdgeId &a_e, const MappingRange &a_range,
                    const EdgeId &b_e, const MappingRange &b_range,
                    const string &s) {
        std::map<VertexId, size_t> vertex_pathlen;
        if (!Consistent(a_e, a_range, b_e, b_range, vertex_pathlen)) {
            INFO("Not consistent");
            return Mapping();
        }

        int start_pos = a_range.initial_range.start_pos;
        int end_pos = b_range.initial_range.start_pos;
        string ss = s.substr(start_pos, min(end_pos, int(s.size()) ) - start_pos);
        int s_len = int(ss.size());
        int path_max_length = cfg_.gap_min_ed;
        int tid = omp_get_thread_num();
        INFO("TID=" << tid << ". Dijkstra: String length " << s_len << " max-len " << path_max_length << " start_pos=" << a_range.initial_range.start_pos << " end_pos=" << b_range.initial_range.start_pos << " start_pos_g=" << a_range.mapped_range.start_pos << " end_pos_g=" << b_range.mapped_range.start_pos);
        if (s_len > 2000 && vertex_pathlen.size() > 100000) {
            DEBUG("Dijkstra on't run: Too big gap or too many paths");
            return Mapping();
        }
        utils::perf_counter perf;
        graph_aligner::DijkstraGapFiller gap_filler(gp_.g, cfg_.gap_cfg,
                ss, a_e, b_e,
                a_range.mapped_range.start_pos,
                b_range.mapped_range.start_pos,
                path_max_length, vertex_pathlen, true);
        gap_filler.CloseGap();
        int score = gap_filler.edit_distance();
        VertexId start_v = gp_.g.EdgeEnd(a_e);
        VertexId end_v = gp_.g.EdgeStart(b_e);

        if (score == std::numeric_limits<int>::max()) {
            //INFO("TID=" << tid << ". Dijkstra didn't find anything")
            INFO("TIME.NO=" << perf.time() << " " << s_len << " " << score << " " << gp_.g.int_id(start_v) << " " << gp_.g.int_id(end_v) << " tid=" << tid);
            score = std::numeric_limits<int>::max();
            return Mapping();
        } else {
            //INFO("TID=" << tid << ". Dijkstra score=" << score << " ssize=" << s_len);
            INFO("TIME.YES=" << perf.time() << " " << s_len << " " << score << " " << gp_.g.int_id(start_v) << " " << gp_.g.int_id(end_v) << " tid=" << tid);
        }
        omnigraph::MappingPath<debruijn_graph::EdgeId> ans = gap_filler.mapping_path();
        if (ans.size() == 1) {
            return Mapping(b_e, MappingRange(Range(b_range.initial_range.start_pos, b_range.initial_range.end_pos),
                                             Range(b_range.mapped_range.start_pos, b_range.mapped_range.end_pos) ), score);
        } else {
            omnigraph::MappingPath<debruijn_graph::EdgeId> path;
            path.push_back(ans.edge_at(0), MappingRange(Range(start_pos + ans.mapping_at(0).initial_range.start_pos, start_pos + ans.mapping_at(0).initial_range.end_pos),
                           Range(ans.mapping_at(0).mapped_range) ));
            for (size_t i = 1; i < ans.size() - 1; ++ i) {
                path.push_back(ans.edge_at(i), MappingRange(Range(0, 0),
                               Range(ans.mapping_at(i).mapped_range) ));
            }
            return Mapping(b_e, b_range, path, score);
        }
    }


    vector<Mapping> FillGapsBetweenPrimers(const std::vector<ReadMapping> &mappings, const io::SingleRead &read) {
        vector<Mapping> resulting_paths;
        vector<Mapping> working_paths;
        VERIFY(mappings[0].e_.size() > 0)
        for (size_t i = 0; i < mappings[0].e_.size(); ++ i) {
            working_paths.push_back(Mapping(mappings[0].e_.at(i), mappings[0].range_.at(i), 0));
        }
        INFO("Start merging")
        utils::perf_counter perf;
        for (size_t i = 1; i < mappings.size(); ++ i) {
            INFO("Iter=" << i << " " << mappings[i].e_.size());
            vector<Mapping> new_working_paths;
            vector<bool> used_working_paths(working_paths.size());
            for (size_t j = 0; j < mappings[i].e_.size(); ++ j) {
                Mapping best_path;
                int best_used = -1;
                for (size_t k = 0; k < working_paths.size(); ++ k) {
                    const MappingRange &a_range = working_paths[k].last_mapping_;
                    const MappingRange &b_range = mappings[i].range_.at(j);
                    const EdgeId &a_e = working_paths[k].last_edge_;
                    const EdgeId &b_e = mappings[i].e_.at(j);
                    Mapping intermediate_path = FillGap(a_e, a_range, b_e, b_range, read.GetSequenceString());

                    if (intermediate_path.score_ != std::numeric_limits<int>::max() && (best_path.score_ == std::numeric_limits<int>::max() || best_path.score_ > intermediate_path.score_ + working_paths[k].score_)) {
                        best_path.path_ = working_paths[k].path_;
                        for (size_t z = 0; z < intermediate_path.path_.size(); ++ z) {
                            best_path.path_.push_back(intermediate_path.path_.edge_at(z), intermediate_path.path_.mapping_at(z));
                        }
                        best_path.last_edge_ = intermediate_path.last_edge_;
                        best_path.last_mapping_ = intermediate_path.last_mapping_;
                        best_path.score_ = intermediate_path.score_ + working_paths[k].score_;
                        best_used = k;
                        //INFO("  Here sz=" << best_path.path_.size() << " score=" << best_path.score_)
                    }
                }

                if (best_path.score_ != std::numeric_limits<int>::max() && best_path.score_ < cfg_.gap_min_ed) {
                    used_working_paths[best_used] = true;
                    new_working_paths.push_back(best_path);
                }
            }
            working_paths = new_working_paths;
        }
        INFO("TIME.FillGaps=" << perf.time());

        for (size_t k = 0; k < working_paths.size(); ++ k) {
            resulting_paths.push_back(working_paths[k]);
        }
        int len_max = 0;
        int ed = 100500;
        std::string ans = "";
        int seq_start = -1;
        int seq_end = -1;
        std::string cur_path = "";
        std::string cur_path_len = "";
        int best_ind = -1;
        for (size_t i = 0; i < resulting_paths.size(); ++ i) {
            Mapping &path = resulting_paths[i];
            path.path_.push_back(path.last_edge_, path.last_mapping_);
            if (path.score_ < cfg_.gap_min_ed) {
                utils::perf_counter perf;
                DEBUG("Grow Ends front")
                int before = path.score_;
                path.score_ += GrowEnds(path.path_, read.GetSequenceString(), false);
                DEBUG("Grow Ends back")
                path.score_ += GrowEnds(path.path_, read.GetSequenceString(), true);
                INFO("TIME.GrowEnds=" << perf.time() << " " << path.score_);
                int len_cur = path.path_.mapping_at(path.path_.size() - 1).initial_range.end_pos - path.path_.mapping_at(0).initial_range.start_pos;
                // if (path.path_.mapping_at(path.path_.size() - 1).initial_range.end_pos - path.path_.mapping_at(0).initial_range.start_pos > 1200) {
                //     ans = "";
                //     for (size_t k = 0; k < path.path_.size(); ++ k) {
                //         ans += std::to_string(path.path_.edge_at(k).int_id())  + ",";
                //     }
                //     INFO("Good path " << path.score_ << " " << before << " " << ans)
                // }
                // INFO("path_sz=" << path.path_.size() << " last_e=" <<path.last_edge_.int_id() << " last_range s=" << path.last_mapping_.initial_range.start_pos << " " << path.last_mapping_.initial_range.end_pos)
                if (path.path_.mapping_at(path.path_.size() - 1).initial_range.end_pos -
                        path.path_.mapping_at(0).initial_range.start_pos > cfg_.alignment_min_length
                        && ed > path.score_) {
                    len_max = path.path_.mapping_at(path.path_.size() - 1).initial_range.end_pos - path.path_.mapping_at(0).initial_range.start_pos;
                    ed = path.score_;
                    best_ind = i;
                }
            }
        }
        if (best_ind >= 0) {
            Mapping &best_path = resulting_paths[best_ind];

            ans = "";
            seq_start = best_path.path_.mapping_at(0).initial_range.start_pos;
            seq_end = best_path.path_.mapping_at(best_path.path_.size() - 1).initial_range.end_pos;
            cur_path = "";
            cur_path_len = "";
            int end_ind = best_path.path_.size();
            int sum = gp_.g.k();
            while (sum - (int) (best_path.path_.mapping_at(end_ind - 1).mapped_range.end_pos - 0) >= 0) {
                sum -= (int) (best_path.path_.mapping_at(end_ind - 1).mapped_range.end_pos - 0);
                end_ind --;
            }
            sum = gp_.g.k();
            for (size_t k = 0; k < end_ind; ++ k) {
                if (sum - (int) (best_path.path_.mapping_at(k).mapped_range.end_pos - best_path.path_.mapping_at(k).mapped_range.start_pos) < 0) {
                    ans += std::to_string(best_path.path_.edge_at(k).int_id())  + ",";
                    cur_path += std::to_string(best_path.path_.edge_at(k).int_id()) + " (" + std::to_string(best_path.path_.mapping_at(k).mapped_range.start_pos) + "," + std::to_string(best_path.path_.mapping_at(k).mapped_range.end_pos) + ") ["
                                + std::to_string(best_path.path_.mapping_at(k).initial_range.start_pos) + "," + std::to_string(best_path.path_.mapping_at(k).initial_range.end_pos) + "], ";
                    cur_path_len += std::to_string(best_path.path_.mapping_at(k).mapped_range.end_pos - best_path.path_.mapping_at(k).mapped_range.start_pos) + ",";
                }
                sum -= (int) (best_path.path_.mapping_at(k).mapped_range.end_pos - best_path.path_.mapping_at(k).mapped_range.start_pos);

            }
            len_max = best_path.path_.mapping_at(best_path.path_.size() - 1).initial_range.end_pos - best_path.path_.mapping_at(0).initial_range.start_pos;
            INFO("Read=" << read.name() << " Primer num=" << mappings.size() << ". paths_num=" << resulting_paths.size() << " maxlen=" << len_max << ". ans=" << ans << ". ed=" << ed)
            if (ed < cfg_.res_min_ed) {

                std::string sum_str = read.name() + "\t" + std::to_string(seq_start) + "\t" + std::to_string(seq_end) + "\t"
                                      + std::to_string(read.size()) +  "\t" + cur_path + "\t" + cur_path_len + "\t" + std::to_string(ed) + "\n";
                #pragma omp critical
                {
                    ofstream myfile;
                    myfile.open(output_file_ + ".tsv", std::ofstream::out | std::ofstream::app);
                    myfile << sum_str;
                    myfile.close();
                }
            }
        }
        return resulting_paths;
    }

    void AlignRead(const io::SingleRead &read) {
        utils::perf_counter perf;
        std::vector<ReadMapping> mappings;
        int num = 0;
        int ind = 0;
        for (auto primer : primers_) {
            int start_pos = -1;
            int end_pos = -1;
            int dist = EditDistance(primer.read_seq_, read.GetSequenceString(), start_pos, end_pos);
            //INFO("dist=" << dist << " start_pos=" << start_pos << " end_pos=" << end_pos << " primer_s=" << primer.start_pos_ << " primer_e="<< primer.end_pos_)
            if (dist <= cfg_.primer_max_ed && abs(start_pos - primer.start_pos_) < cfg_.position_diff && abs(end_pos - primer.end_pos_) < cfg_.position_diff) {
                num += 1;
                INFO("dist=" << dist << " start_pos=" << start_pos << " end_pos=" << end_pos <<
                     " primer_s=" << primer.start_pos_ << " primer_e=" << primer.end_pos_ << " sz=" << primer.range_.size())
                for (size_t i = 0; i < primer.range_.size(); ++ i) {
                    //INFO(ind << " e=" << primer.e_[i].int_id() << " " << primer.range_[i].mapped_range.start_pos << " " << primer.range_[i].mapped_range.end_pos)
                    primer.range_[i] = MappingRange(Range(start_pos, end_pos),
                                                    Range(primer.range_[i].mapped_range.start_pos, primer.range_[i].mapped_range.end_pos));
                }
                mappings.push_back(primer);
            }
            ind ++;
        }
        int tid = omp_get_thread_num();
        INFO("TID=" << tid << " Read=" << read.name() << " Primer num=" << num)
        if (mappings.size() > 0) {
            vector<Mapping> res = FillGapsBetweenPrimers(mappings, read);
        } else {
            vector<Mapping>();
        }
        INFO("TIME.Read2=" << perf.time() << " name=" << read.name() )
    }
};

class Rna16SMapper {
public:
    Rna16SMapper(const debruijn_graph::conj_graph_pack &gp,
                 const rna16S_mapping::RnaAlignerConfig &cfg,
                 const string &output_file):
        gp_(gp), cfg_(cfg), output_file_(output_file) {}

private:
    const debruijn_graph::conj_graph_pack &gp_;
    const rna16S_mapping::RnaAlignerConfig cfg_;
    const string &output_file_;
    
    PrimerAligner primer_aligner;
    SequenceAligner sequence_aligner;
};

void LoadSequences(const string &sequence_fasta, std::vector<io::SingleRead> &wrappedreads) {
    io::ReadStreamList<io::SingleRead> streams;
    streams.push_back(make_shared<io::FileReadStream>(sequence_fasta)); //make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(sequence_fasta)));
    io::SingleStreamPtr sstream = io::MultifileWrap(streams);

    while (!sstream->eof()) {
        io::SingleRead read;
        *sstream >> read;
        wrappedreads.push_back(std::move(read));
    }
    INFO("Loaded sequences from " << sequence_fasta);
}

void Launch(RnaAlignerConfig cfg, const string &output_file, int threads) {
    string tmpdir = fs::make_temp_dir(fs::current_dir(), "tmp");
    debruijn_graph::conj_graph_pack gp(cfg.K, tmpdir, 0);
    debruijn_graph::graphio::ScanGraphPack(cfg.path_to_graphfile, gp);
    INFO("Loaded graph with " << gp.g.size() << " vertices");

    std::vector<io::SingleRead> wrappedprimers;
    std::vector<io::SingleRead> wrappedreads;
    LoadSequences(cfg.path_to_primers, wrappedprimers);
    LoadSequences(cfg.path_to_sequences, wrappedreads);

    ofstream myfile;
    myfile.open(output_file + ".gpa", std::ofstream::out | std::ofstream::app);
    myfile << "H\n";
    myfile.close();

    aligner.AnalyzePrimers(wrappedprimers, threads);

    SequenceAligner aligner(gp, cfg, output_file);
    #pragma omp parallel num_threads(threads)
    #pragma omp for
    for (size_t i = 0 ; i < wrappedreads.size(); ++i) {
        aligner.AlignRead(wrappedreads[i]);
    }

    fs::remove_dir(tmpdir);
}
} // namespace rna16S_mapping

int main(int argc, char **argv) {
    unsigned nthreads;
    std::string output_file;
    std::vector<std::string> input;

    cxxopts::Options options(argv[0], " <K - graph K-mer length> \
                             <Path to graph saves> \
                             <YAML-config sequences and graph description> - Tool for SILVA alignment on graph ");
    options.add_options()
    ("o,outfile", "Output file prefix", cxxopts::value<std::string>(output_file)->default_value("./silva_result"), "prefix")
    ("t,threads", "# of threads to use", cxxopts::value<unsigned>(nthreads)->default_value(std::to_string(min(omp_get_max_threads(), 16))), "num")
    ("h,help", "Print help");

    options.add_options("Input")
    ("positional", "", cxxopts::value<std::vector<std::string>>(input));

    options.parse_positional("positional");
    options.parse(argc, argv);
    if (options.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    if (options.count("positional") < 3) {
        std::cerr << "ERROR: No input proper input" << std::endl << std::endl;
        std::cout << options.help() << std::endl;
        exit(-1);
    }

    std::string cfg = input[2];
    create_console_logger();
    auto buf = llvm::MemoryBuffer::getFile(cfg);
    VERIFY_MSG(buf, "Failed to load config file " + cfg);
    llvm::yaml::Input yin(*buf.get());
    rna16S_mapping::RnaAlignerConfig config;
    yin >> config;
    config.K = std::stoi(input[0]);
    config.path_to_graphfile = input[1];
    omp_set_num_threads(nthreads);

    rna16S_mapping::Launch(config, output_file, nthreads);
    return 0;
}
