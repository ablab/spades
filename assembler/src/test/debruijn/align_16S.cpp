//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/standard_base.hpp"
#include "utils/stl_utils.hpp"
#include "utils/logger/log_writers.hpp"

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
#include "pipeline/config_struct.hpp"

#include "assembly_graph/core/graph.hpp"
#include "modules/alignment/sequence_mapper.hpp"

#include "io/reads/io_helper.hpp"
#include "common/assembly_graph/core/coverage.hpp"

//#include "modules/alignment/pacbio/gap_dijkstra.hpp"
#include "modules/alignment/long_read_mapper.hpp"
#include "modules/alignment/short_read_mapper.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "io/reads/multifile_reader.hpp"

#include <iostream>
#include <fstream>

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

struct ReadMapping {

    void ExtractPositions() {
        start_pos_ = -1;
        end_pos_ = -1;
        string name = read_;
        int open_bracket_ind = -1;
        int close_bracket_ind = -1;
        int comma_ind = -1;
        for (int i = 0; i < name.size(); ++ i) {
            if (name[i] == '[') {
                open_bracket_ind = i;
            }
            if (name[i] == ']') {
                close_bracket_ind = i;
            }
            if (name[i] == ',') {
                comma_ind = i;
            }
        }
        start_pos_ = std::stoi(name.substr(open_bracket_ind + 1, comma_ind - open_bracket_ind));
        end_pos_ = std::stoi(name.substr(comma_ind + 1, close_bracket_ind - comma_ind));
        INFO("Name=" << name << " start_pos_=" << start_pos_ << " end_pos_=" << end_pos_);
    }

    ReadMapping(const std::string &read, const std::string &read_seq, const int &v, const std::vector<EdgeId> &e, const std::vector<MappingRange> &range) 
        : read_(read), read_seq_(read_seq), v_(v), e_(e), range_(range){
        ExtractPositions();
    }

    bool operator < (const ReadMapping &mapping) const {
        return (this->start_pos_ < mapping.start_pos_ || (this->start_pos_ == mapping.start_pos_ && this->end_pos_ < mapping.start_pos_) );
    }

    std::string read_;
    std::string read_seq_;
    int v_;
    std::vector<EdgeId> e_;
    std::vector<MappingRange> range_;
    int start_pos_;
    int end_pos_;
};


struct Mapping {

    Mapping(const EdgeId &e, const MappingRange &range): last_edge_(e), last_mapping_(range), score_(-1){}

    Mapping(const EdgeId &e, const MappingRange &range, int score): last_edge_(e), last_mapping_(range), score_(score){}

    Mapping(const EdgeId &e, const MappingRange &range, const omnigraph::MappingPath<debruijn_graph::EdgeId> &path, int score)
            : last_edge_(e), last_mapping_(range), path_(path), score_(score){}

    Mapping(): score_(-1){}

    omnigraph::MappingPath<debruijn_graph::EdgeId> path_;
    EdgeId last_edge_;
    MappingRange last_mapping_;
    int score_;
};


class SequenceAligner {
private:    
    const conj_graph_pack &gp_;
    //const pacbio::PacBioMappingIndex<Graph> pac_index_;
    const string &output_file_;
    const config::debruijn_config::pacbio_processor &pb_config_;
    int ed_threshold_;
    int res_ed_threshold_;
    int min_length_;
    int primer_threshold_;
    int graph_threshold_;

    std::vector<ReadMapping> primers_;

public:
    SequenceAligner(const conj_graph_pack &gp, 
                 const alignment::BWAIndex::AlignmentMode mode,
                 const config::debruijn_config::pacbio_processor &pb,
                 const string &output_file,
                 int ed_threshold, int res_ed_threshold, int min_length, int primer_threshold, int graph_threshold):
      gp_(gp), output_file_(output_file), pb_config_(pb)
      , ed_threshold_(ed_threshold), res_ed_threshold_(res_ed_threshold)
      , min_length_(min_length), primer_threshold_(primer_threshold)
      , graph_threshold_(graph_threshold) {}

    
    bool IsCanonical(EdgeId e) const {
        return e <= gp_.g.conjugate(e);
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : gp_.g.conjugate(e);
    }

    int EditDistance(const string &a, const string &b, int &start_pos, int &end_pos) {
        //INFO("Before ed")
        int a_len = (int) a.length();
        int b_len = (int) b.length();
        VERIFY(a_len > 0);
        VERIFY(b_len > 0);
        edlib::EdlibEqualityPair additionalEqualities[36] = {{'U', 'T'}
                                                , {'R', 'A'}, {'R', 'G'}
                                                , {'Y', 'C'}, {'Y', 'T'}, {'Y', 'U'}
                                                , {'K', 'G'}, {'K', 'T'}, {'K', 'U'}
                                                , {'M', 'A'}, {'M', 'C'}
                                                , {'S', 'C'}, {'S', 'G'}
                                                , {'W', 'A'}, {'W', 'T'}, {'W', 'U'}
                                                , {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'B', 'U'}
                                                , {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'D', 'U'}
                                                , {'H', 'A'}, {'H', 'C'}, {'H', 'T'}, {'H', 'U'}
                                                , {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
                                                , {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, {'N', 'U'} };
        edlib::EdlibAlignResult result = edlib::edlibAlign(a.c_str(), a_len, b.c_str(), b_len
                                                       , edlib::edlibNewAlignConfig(a_len, edlib::EDLIB_MODE_HW, edlib::EDLIB_TASK_LOC,
                                                                             additionalEqualities, 36));
        int score = -1;
        if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
            if (result.numLocations > 0) {
                score = result.editDistance;
                start_pos = result.startLocations[0];
                end_pos = result.endLocations[0];
            } else {
                WARN("EditDistance: something wrong with edlib result");
            }
        }
        edlib::edlibFreeAlignResult(result);
        return score;
    }

    void AlignPrimer(const io::SingleRead &read, int &v, std::vector<EdgeId> &e, std::vector<MappingRange> &range) {
        v = graph_threshold_;
        int k = 10;
        for (auto it = gp_.g.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            EdgeId eid = *it;
            std::string edge_str = gp_.g.EdgeNucls(eid).str();
            int start_pos = 0;
            int end_pos = 1;
            int dist = EditDistance(read.GetSequenceString(), edge_str, start_pos, end_pos);
            if (dist == -1) continue;
            if (v >= dist && start_pos < gp_.g.length(eid)) {
                e.push_back(eid);
                range.push_back(MappingRange(Range(0, read.size()), Range(start_pos, end_pos + 1)) );
            }
        }

        std::string ans = " dist=" +  std::to_string(v) + " num=" +  std::to_string(e.size());
        // for (EdgeId eid: e) {
        //     ans += " " + std::to_string(eid.int_id()) + ",";
        // }
        // ans += "\n";
        // for (MappingRange r: range) {
        //     ans += " " + std::to_string(r.mapped_range.start_pos) + "-" + std::to_string(r.mapped_range.end_pos) + ",";
        // }
        INFO("Primer name=" << read.name() <<" seq=" << read.GetSequenceString() << " " << ans);
    }

    void PreparePrimers(std::vector<io::SingleRead> &wrappedprimers, int threads) {
#pragma omp parallel num_threads(threads)
#pragma omp for
        for (size_t i =0 ; i < wrappedprimers.size(); ++i) {
            int v = 0;
            std::vector<EdgeId> e;
            std::vector<MappingRange> range;
            AlignPrimer(wrappedprimers[i], v, e, range);
            const std::string &name = wrappedprimers[i].name();
            const std::string &seq = wrappedprimers[i].GetSequenceString();
            INFO("i=" << i)
#pragma omp critical
{
            if (e.size() > 0) {
                primers_.push_back(ReadMapping(name, seq, v, e, range));
            }
}   
        }   
        std::sort(primers_.begin(), primers_.end());
    }

    void PrepareInitialState(omnigraph::MappingPath<debruijn_graph::EdgeId> &path, const string &s, bool forward, string &ss, EdgeId &start_e, int &start_pos, int &seq_start_pos) const {
        if (forward){
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
                                    , {'N', 'N'} };
            ss = "";
            int num = 0;
            for (int i = c_ss.size() - 1; i >= 0; -- i) {
                ss = ss + nucs[c_ss[i]];
                num ++;
            }
            
            DEBUG("c_ss=" << c_ss.size() << " ss=" << ss.size() << " num=" << num);
        }
    }

    void UpdatePath(omnigraph::MappingPath<debruijn_graph::EdgeId> &path, std::vector<EdgeId> &ans, int start_pos, int end_pos, int seq_start_pos, int seq_end_pos, bool forward) const {
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
            while (cur_ind >= 0 && start - (int) gp_.g.length(ans[cur_ind]) > 0){
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
        int score = ed_threshold_;
        if (s_len > 2000) {
            DEBUG("EdgeDijkstra: sequence is too long " << s_len)
            return 0;
        }
        if (s_len < (int) gp_.g.k()) {
            DEBUG("EdgeDijkstra: sequence is too small " << s_len)
            return 0;
        }
        DEBUG(" EdgeDijkstra: String length " << s_len << " max-len " << score << " start_pos=" << start_pos << " start_e=" << start_e.int_id() << " e_len=" << gp_.g.length(start_e) + gp_.g.k() );
        pacbio::DijkstraEndsReconstructor algo = pacbio::DijkstraEndsReconstructor(gp_.g, ss, start_e, start_pos, score);
        algo.CloseGap();
        score = algo.GetEditDistance();
        if (score == -1){
            DEBUG("EdgeDijkstra didn't find anything edge=" << start_e.int_id() << " s_start=" << start_pos << " seq_len=" << ss.size())
            return 0;
        } else {
            DEBUG("EdgeDijkstra found path edge=" << start_e.int_id() << " s_start=" << start_pos << " seq_len=" << ss.size() << " score=" << score )
        }

        std::vector<EdgeId> ans = algo.GetPath();
        int end_pos = algo.GetPathEndPosition();
        int seq_end_pos = algo.GetSeqEndPosition();
        DEBUG("graph_end=" << end_pos << " seq_end_pos=" << seq_end_pos << " anssz=" << ans.size());
        UpdatePath(path, ans, start_pos, end_pos, seq_start_pos, seq_end_pos, forward);
        return score;
    }

    bool Consistent(const EdgeId &a_e, const MappingRange &a_range, const EdgeId &b_e, const MappingRange &b_range, std::map<VertexId, size_t> &vertex_pathlen) {
        int start_pos = a_range.initial_range.start_pos;
        int end_pos = b_range.initial_range.start_pos;
        int seq_len = -start_pos + end_pos;
        int path_max_length = (int) ((double) seq_len * pb_config_.path_limit_stretching);
        //INFO("path_max_length=" << path_max_length << " seq_len=" << seq_len);
        VertexId start_v = gp_.g.EdgeEnd(a_e);
        VertexId end_v = gp_.g.EdgeStart(b_e);
        if (end_pos < start_pos) {
            INFO("modifying limits because of some bullshit magic, seq length 0 " << start_pos << " " << end_pos)
            end_pos = start_pos;
            return false;
        }
        
        auto path_searcher_b = omnigraph::DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(gp_.g, path_max_length);
        path_searcher_b.Run(end_v);
        auto path_searcher = omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(gp_.g, path_max_length);
        path_searcher.Run(start_v);
        auto reached_vertices_b = path_searcher_b.ProcessedVertices();
        auto reached_vertices = path_searcher.ProcessedVertices();

        for (auto j_iter = reached_vertices_b.begin(); j_iter != reached_vertices_b.end(); ++j_iter) {
            if (reached_vertices.count(*j_iter) > 0){
                    vertex_pathlen[*j_iter] = path_searcher_b.GetDistance(*j_iter);
            }
        }

        if (a_e == b_e){
            if (b_range.mapped_range.start_pos < a_range.mapped_range.start_pos) {
                if (vertex_pathlen.count(end_v) == 0){
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

        if (vertex_pathlen.count(end_v) == 0) {
            //INFO(" Cant get")
            return false;
        }
        //INFO(" Consistent")
        return true;
    }

    Mapping FillGap(const EdgeId &a_e, const MappingRange &a_range, const EdgeId &b_e, const MappingRange &b_range, const string &s) {
        std::map<VertexId, size_t> vertex_pathlen;
        if (!Consistent(a_e, a_range, b_e, b_range, vertex_pathlen)) {
            //INFO("Not consistent");
            return Mapping();
        }

        int start_pos = a_range.initial_range.start_pos;
        int end_pos = b_range.initial_range.start_pos;
        string ss = s.substr(start_pos, min(end_pos, int(s.size()) ) - start_pos);
        int s_len = int(ss.size());
        int path_max_length = ed_threshold_;
        //INFO(" Dijkstra: String length " << s_len << " max-len " << path_max_length << " start_pos=" <<a_range.initial_range.start_pos << " end_pos=" << b_range.initial_range.start_pos << " start_pos_g=" <<a_range.mapped_range.start_pos << " end_pos_g=" << b_range.mapped_range.start_pos);
        if (s_len > 2000 && vertex_pathlen.size() > 100000){
            DEBUG("Dijkstra on't run: Too big gap or too many paths");
            return Mapping();
        }
        pacbio::DijkstraGapFiller gap_filler = pacbio::DijkstraGapFiller(gp_.g, ss, a_e, b_e, a_range.mapped_range.start_pos, b_range.mapped_range.start_pos, path_max_length, vertex_pathlen);
        gap_filler.CloseGap();
        int score = gap_filler.GetEditDistance();
        if (score == -1){
            DEBUG("Dijkstra didn't find anything")
            score = -1;
            return Mapping();
        } else {
            DEBUG("Dijkstra score=" << score << " ssize=" << s_len);
        }
        omnigraph::MappingPath<debruijn_graph::EdgeId> ans = gap_filler.GetMappingPath();
        if (ans.size() == 1){
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
            for (size_t i = 0; i < mappings[0].e_.size(); ++ i){
                working_paths.push_back(Mapping(mappings[0].e_.at(i), mappings[0].range_.at(i), 0));
            }
            for (size_t i = 1; i < mappings.size(); ++ i){
                //INFO("Iter=" << i << " " << mappings[i].e_.size());
                vector<Mapping> new_working_paths;
                vector<bool> used_working_paths(working_paths.size());
                for (size_t j = 0; j < mappings[i].e_.size(); ++ j) {
                    Mapping best_path;
                    int best_used = -1;
                    int bad = -1;
                    int good = -1;
                    for (size_t k = 0; k < working_paths.size(); ++ k) {
                        const MappingRange &a_range = working_paths[k].last_mapping_;
                        const MappingRange &b_range = mappings[i].range_.at(j);
                        const EdgeId &a_e = working_paths[k].last_edge_;
                        const EdgeId &b_e = mappings[i].e_.at(j);
                        Mapping intermediate_path = FillGap(a_e, a_range, b_e, b_range, read.GetSequenceString());
                        // if (a_e.int_id() == 19108181|| b_e.int_id() == 2669802) {
                        //     INFO( "a_e=" << a_e.int_id() << " b_e=" << b_e.int_id() << " " << best_path.score_ << " " << intermediate_path.score_  << " " << working_paths[k].score_ << " " << a_range.mapped_range.start_pos << " " << a_range.mapped_range.end_pos << " " << b_range.mapped_range.start_pos << " " << b_range.mapped_range.end_pos  )
                        // }
                        // if (a_e.int_id() == 2669802 || b_e.int_id() == 19108181) {
                        //     INFO( "a_e=" << a_e.int_id() << " b_e=" << b_e.int_id() << " " << best_path.score_ << " " << intermediate_path.score_  << " " << working_paths[k].score_ << " " << a_range.mapped_range.start_pos << " " << a_range.mapped_range.end_pos << " " << b_range.mapped_range.start_pos << " " << b_range.mapped_range.end_pos  )
                        // }

                        // if (a_e.int_id() == 19226365 || b_e.int_id() == 19226365) {
                        //     INFO( "a_e=" << a_e.int_id() << " b_e=" << b_e.int_id() << " " << best_path.score_ << " " << intermediate_path.score_  << " " << working_paths[k].score_ << " " << a_range.mapped_range.start_pos << " " << a_range.mapped_range.end_pos << " " << b_range.mapped_range.start_pos << " " << b_range.mapped_range.end_pos  )
                        // }

                        // if (a_e.int_id() == 3001321 || b_e.int_id() == 3001321) {
                        //     INFO( "a_e=" << a_e.int_id() << " b_e=" << b_e.int_id() << " " << best_path.score_ << " " << intermediate_path.score_  << " " << working_paths[k].score_ << " " << a_range.mapped_range.start_pos << " " << a_range.mapped_range.end_pos << " " << b_range.mapped_range.start_pos << " " << b_range.mapped_range.end_pos  )
                        // }

                        if (intermediate_path.score_ != -1 && (best_path.score_ == -1 || best_path.score_ > intermediate_path.score_ + working_paths[k].score_)) {
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
    
                    if (best_path.score_ != -1 && best_path.score_ < ed_threshold_) {
                        used_working_paths[best_used] = true;
                        new_working_paths.push_back(best_path);
                    }
                }
                working_paths = new_working_paths;
            }
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
                if (path.score_ < ed_threshold_) {
                    DEBUG("Grow Ends front")
                    int before = path.score_;
                    path.score_ += GrowEnds(path.path_, read.GetSequenceString(), false);
                    DEBUG("Grow Ends back")
                    path.score_ += GrowEnds(path.path_, read.GetSequenceString(), true);
                    int len_cur = path.path_.mapping_at(path.path_.size() - 1).initial_range.end_pos - path.path_.mapping_at(0).initial_range.start_pos;
                    // if (path.path_.mapping_at(path.path_.size() - 1).initial_range.end_pos - path.path_.mapping_at(0).initial_range.start_pos > 1200) {
                    //     ans = "";
                    //     for (size_t k = 0; k < path.path_.size(); ++ k) {
                    //         ans += std::to_string(path.path_.edge_at(k).int_id())  + ",";
                    //     }
                    //     INFO("Good path " << path.score_ << " " << before << " " << ans)
                    // }
                    // INFO("path_sz=" << path.path_.size() << " last_e=" <<path.last_edge_.int_id() << " last_range s=" << path.last_mapping_.initial_range.start_pos << " " << path.last_mapping_.initial_range.end_pos)
                    if (path.path_.mapping_at(path.path_.size() - 1).initial_range.end_pos - path.path_.mapping_at(0).initial_range.start_pos > min_length_ && ed > path.score_) {
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
                    int sum = 56;
                    while (sum - (int) (best_path.path_.mapping_at(end_ind - 1).mapped_range.end_pos - 0) >= 0) {
                        sum -= (int) (best_path.path_.mapping_at(end_ind - 1).mapped_range.end_pos - 0);
                        end_ind --;
                    }   
                    sum = 56;
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
                    INFO("Read="<< read.name() << " Primer num=" << mappings.size() << ". paths_num=" << resulting_paths.size() << " maxlen=" << len_max << ". ans=" << ans << ". ed=" << ed)
                    if (ed < res_ed_threshold_) {
                        
                        std::string sum_str = read.name() + "\t" + std::to_string(seq_start) + "\t" + std::to_string(seq_end) + "\t" 
                                                             + std::to_string(read.size())+  "\t" + cur_path + "\t" + cur_path_len + "\t" + std::to_string(ed) + "\n";
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
        int threshold = 200;
        int ind = 0;
        for (auto primer: primers_) {
            int start_pos = -1;
            int end_pos = -1;
            int dist = EditDistance(primer.read_seq_, read.GetSequenceString(), start_pos, end_pos);
            //INFO("dist=" << dist << " start_pos=" << start_pos << " end_pos=" << end_pos << " primer_s=" << primer.start_pos_ << " primer_e="<< primer.end_pos_)
            if (dist <= primer_threshold_ && abs(start_pos - primer.start_pos_) < threshold && abs(end_pos - primer.end_pos_) < threshold){
                num += 1;
                INFO("dist=" << dist << " start_pos=" << start_pos << " end_pos=" << end_pos << " primer_s=" << primer.start_pos_ << " primer_e="<< primer.end_pos_)
                for (size_t i = 0; i < primer.range_.size(); ++ i){
                    //INFO(ind << " e=" << primer.e_[i].int_id() << " " << primer.range_[i].mapped_range.start_pos << " " << primer.range_[i].mapped_range.end_pos)
                    primer.range_[i] = MappingRange(Range(start_pos, end_pos), Range(primer.range_[i].mapped_range.start_pos, primer.range_[i].mapped_range.end_pos));
                }
                mappings.push_back(primer);
            }
            ind ++;
        }
        int tid = omp_get_thread_num();
        INFO("TID=" << tid << " Read="<< read.name() << " Primer num=" << num)
        if (mappings.size() > 0) {
            vector<Mapping> res = FillGapsBetweenPrimers(mappings, read);
        } else {
            vector<Mapping>();
        }
        INFO("Read2="<< read.name() << " time=" << perf.time())
    }
};

config::debruijn_config::pacbio_processor InitializePacBioProcessor() {
    config::debruijn_config::pacbio_processor pb;  
    pb.bwa_length_cutoff = 0; //500
    pb.compression_cutoff = 0.6; // 0.6 
    pb.path_limit_stretching = 1.3; //1.3
    pb.path_limit_pressing = 0.7;//0.7
    pb.max_path_in_dijkstra = 15000; //15000
    pb.max_vertex_in_dijkstra = 15000; //2000
//gap_closer
    pb.long_seq_limit = 400; //400
    pb.pacbio_min_gap_quantity = 2; //2
    pb.contigs_min_gap_quantity = 1; //1
    pb.max_contigs_gap_length = 10000; // 10000
    return pb;
}

void LoadReads(const string &sequence_fasta, std::vector<io::SingleRead> &wrappedreads) {
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

void Launch(size_t K, const string &saves_path, const string &primer_fasta, const string &sequence_fasta, const string &mapper_type, const string &output_file, int threads) {
    conj_graph_pack gp(K, "tmp3", 0);
    graphio::ScanGraphPack(saves_path, gp);
    INFO("Loaded graph with " << gp.g.size() << " vertices");

    std::vector<io::SingleRead> wrappedprimers;
    std::vector<io::SingleRead> wrappedreads;
    LoadReads(primer_fasta, wrappedprimers);
        
    LoadReads(sequence_fasta, wrappedreads);
    
    alignment::BWAIndex::AlignmentMode mode;
    if (mapper_type == "pacbio"){
        mode = alignment::BWAIndex::AlignmentMode::PacBio;
    } else if (mapper_type == "nanopore"){
        mode = alignment::BWAIndex::AlignmentMode::Ont2D;
    } else if (mapper_type == "16S"){
        mode = alignment::BWAIndex::AlignmentMode::Rna16S;
    } else {
        mode = alignment::BWAIndex::AlignmentMode::Default;
    }

    config::debruijn_config::pacbio_processor pb = InitializePacBioProcessor();
    int ed_threshold = 40;
    int res_ed_threshold = 200;
    int min_length = 1200;
    int primer_threshold = 1;
    int graph_threshold = 3;
    SequenceAligner aligner(gp, mode, pb, output_file, ed_threshold, res_ed_threshold, min_length, primer_threshold, graph_threshold); 
    INFO("SequenceAligner created");

    ofstream myfile;
    myfile.open(output_file + ".gpa", std::ofstream::out | std::ofstream::app);
    myfile << "H\n";
    myfile.close();

    aligner.PreparePrimers(wrappedprimers, threads);

#pragma omp parallel num_threads(threads)
#pragma omp for
    for (size_t i =0 ; i < wrappedreads.size(); ++i) {
        aligner.AlignRead(wrappedreads[i]);
    }
}
}

int main(int argc, char **argv) {
    if (argc < 7) {
        cout << "Usage: longreads_aligner <K>"
             << " <saves path> <primer file(fasta/fastq)> <sequences file (fasta/fastq)> <mapper type: {pacbio, nanopore, default} > <ouput-prefix> {threads-num(16)}" << endl;
        exit(1);
    }
    int threads = 16;
    if (argc == 8) {
        threads = std::stoi(argv[7]);
    }
    omp_set_num_threads(threads);
    create_console_logger();
    size_t K = std::stoll(argv[1]);
    string saves_path = argv[2];
    INFO("Load graph from " << saves_path);
    string primer_file = argv[3];
    INFO("Load primers from " << primer_file);
    string sequence_file = argv[4];
    INFO("Load sequences from " << sequence_file);
    string mapper_type = argv[5];
    INFO("Mapper type " << mapper_type);
    string output_file = argv[6];
    debruijn_graph::Launch(K, saves_path, primer_file, sequence_file, mapper_type, output_file, threads);
    return 0;
}
