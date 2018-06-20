//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <algorithm>
#include <vector>
#include <set>

#include "assembly_graph/index/edge_multi_index.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"

#include "modules/alignment/edge_index_refiller.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "modules/alignment/gap_info.hpp"

#include "modules/alignment/pacbio/pacbio_read_structures.hpp"
#include "modules/alignment/pacbio/gap_filler.hpp"

namespace pacbio {
enum {
    UNDEF_COLOR = -1,
    DELETED_COLOR = -2
};

struct PathRange {
    size_t seq_start;
    size_t seq_end;

    size_t edge_start;
    size_t edge_end;  
};

struct OneReadMapping {
    std::vector<vector<debruijn_graph::EdgeId>> main_storage;
    std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId>> bwa_paths;
    std::vector<PathRange> read_ranges;
    std::vector<GapDescription> gaps;
    OneReadMapping(const std::vector<vector<debruijn_graph::EdgeId>> &main_storage_,
                   const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId>> &bwa_paths_,
                   const std::vector<GapDescription>& gaps_,
                   const std::vector<PathRange> &read_ranges_) :
            main_storage(main_storage_), bwa_paths(bwa_paths_), gaps(gaps_), read_ranges(read_ranges_){}
};

template<class Graph>
class PacBioMappingIndex {
public:
    typedef std::set<KmerCluster<Graph>> ClustersSet;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

private:
    DECL_LOGGER("PacIndex")

    const Graph &g_;

    static const int LONG_ALIGNMENT_OVERLAP = 300;
    static const size_t SHORT_SPURIOUS_LENGTH = 500;
    mutable std::map<std::pair<VertexId, VertexId>, size_t> distance_cashed_;
    size_t read_count_;
    debruijn_graph::config::pacbio_processor pb_config_;

    alignment::BWAReadMapper<Graph> bwa_mapper_;

    graph_aligner::GapClosingConfig gap_cfg_;

public:

    PacBioMappingIndex(const Graph &g,
                        debruijn_graph::config::pacbio_processor pb_config, 
                        alignment::BWAIndex::AlignmentMode mode, 
                        graph_aligner::GapClosingConfig gap_cfg = graph_aligner::GapClosingConfig())
            : g_(g),
              pb_config_(pb_config),
              bwa_mapper_(g, mode, pb_config.bwa_length_cutoff),
              gap_cfg_(gap_cfg){
        DEBUG("PB Mapping Index construction started");
        DEBUG("Index constructed");
        read_count_ = 0;
    }

    bool similar_in_graph(const MappingInstance &a, const MappingInstance &b,
                          int shift = 0) const {
        if (b.read_position + shift < a.read_position) {
            return similar_in_graph(b, a, -shift);
        } else if (b.read_position == a.read_position) {
            return (abs(int(b.edge_position) + shift - int(a.edge_position)) < 2);
        } else {
//3 to allow small deletion in read on the consecutive edges
            return ((b.edge_position + shift - a.edge_position) * pb_config_.compression_cutoff <= max((b.read_position - a.read_position), 3));
        }
    }

    typename omnigraph::MappingPath<EdgeId> FilterShortAlignments(typename omnigraph::MappingPath<EdgeId> mapped_path) const {
        omnigraph::MappingPath<EdgeId> res;
        size_t length_cutoff = pb_config_.internal_length_cutoff;
        /*
        size_t maxx = 0;
        size_t minn = std::numeric_limits<size_t>::max();
        size_t mini = 0;
        size_t maxi = 0;
        for (size_t i = 0; i < mapped_path.size(); i++) {
            size_t pos = mapped_path[i].second.initial_range.end_pos + mapped_path[i].second.initial_range.start_pos;
            if (minn > pos) {
                minn = pos;
                mini = i;
            }
            if (maxx < pos) {
                maxx = pos;
                maxi = i;
            }
        } */
        for (size_t i = 0; i < mapped_path.size(); i++) {
            size_t rlen = mapped_path[i].second.initial_range.size();
//left and right ends of ranges;
//TODO:: think whether it is right condition
            if (rlen > length_cutoff/* || (rlen > g_.k() && (i == mini || i == maxi))*/){
                vector<pair<size_t, int> > range_limits;
                for (size_t j = 0; j < mapped_path.size(); j++) {
                    if (i != j) {
                        if (mapped_path[i].second.initial_range.Intersect(mapped_path[j].second.initial_range)) {
                            size_t pos_start = std::max (mapped_path[i].second.initial_range.start_pos, mapped_path[j].second.initial_range.start_pos)
                                 - mapped_path[i].second.initial_range.start_pos;
                            size_t pos_end = std::min (mapped_path[i].second.initial_range.end_pos, mapped_path[j].second.initial_range.end_pos)
                                 - mapped_path[i].second.initial_range.start_pos;
                            range_limits.push_back(make_pair(pos_start, 1));
                            range_limits.push_back(make_pair(pos_end, -1));
                        }
                    }
                }
                sort(range_limits.begin(), range_limits.end());
                size_t current_cover = 0;
                if (range_limits.size() > 0) {
//may be negative if some ranges are zero-sized
                    int covered_lays = range_limits[0].second;
                    for (size_t j = 1; j < range_limits.size(); j++) {
                        if (covered_lays > 0) {
                            current_cover += range_limits[j].first - range_limits[j - 1].first;
                        }
                        covered_lays += range_limits[j].second;
                    }
                }
                if (current_cover * 2 <  rlen)
                    res.push_back(mapped_path[i].first, mapped_path[i].second);
                else {
                    DEBUG ("Filtering " << g_.int_id(mapped_path[i].first) << " " << mapped_path[i].second)
                }
            }
        }
        return res;
    }

    typename omnigraph::MappingPath<EdgeId> FilterSpuriousAlignments(typename omnigraph::MappingPath<EdgeId> mapped_path, size_t seq_len) const {
        omnigraph::MappingPath<EdgeId> res;

        for (size_t i = 0; i < mapped_path.size(); i++) {
            omnigraph::MappingRange current = mapped_path[i].second;
//Currently everything in kmers
            size_t expected_additional_left = std::min(current.initial_range.start_pos, current.mapped_range.start_pos);
            size_t expected_additional_right = std::min(seq_len - current.initial_range.end_pos - g_.k(),
                                                       g_.length(mapped_path[i].first) - current.mapped_range.end_pos);

            size_t rlen =
                    mapped_path[i].second.initial_range.end_pos - mapped_path[i].second.initial_range.start_pos;

//FIXME more reasonable condition
//g_.k() to compare sequence lengths, not kmers
            if (rlen < SHORT_SPURIOUS_LENGTH && (rlen + g_.k()) * 2 < expected_additional_left + expected_additional_right) {
                DEBUG ("Skipping spurious alignment " << i << " on edge " << mapped_path[i].first);
            } else
                res.push_back(mapped_path[i].first, mapped_path[i].second);
        }

//        if (res.size() != mapped_path.size()) {
        if (true) {
            DEBUG("Seq len " << seq_len);
            for (const auto &e_mr : mapped_path) {
                EdgeId e = e_mr.first;
                omnigraph::MappingRange mr = e_mr.second;
                DEBUG("Alignments were" << g_.int_id(e) << " e_start=" << mr.mapped_range.start_pos << " e_end=" <<
                     mr.mapped_range.end_pos
                     << " r_start=" << mr.initial_range.start_pos << " r_end=" << mr.initial_range.end_pos);
            }
        }
        return res;
    }

    ClustersSet GetBWAClusters(const io::SingleRead &read) const {
        DEBUG("BWA started")
        ClustersSet res;
        Sequence s = read.sequence();
        if (s.size() < g_.k()) {
            return res;
        }
        omnigraph::MappingPath<EdgeId> m_path = FilterSpuriousAlignments(bwa_mapper_.MapSequence(s), s.size());
        omnigraph::MappingPath<EdgeId> mapped_path = FilterShortAlignments(m_path);

        TRACE(read_count_ << " read_count_");
        TRACE("BWA ended")
        DEBUG(mapped_path.size() <<"  clusters");
        for (const auto &e_mr : mapped_path) {
            EdgeId e = e_mr.first;
            omnigraph::MappingRange mr = e_mr.second;
            DEBUG("ReadName=" << read.name() << " BWA loading edge=" << g_.int_id(e) << " e_start=" << mr.mapped_range.start_pos << " e_end=" << mr.mapped_range.end_pos 
                                                                    << " r_start=" << mr.initial_range.start_pos << " r_end=" << mr.initial_range.end_pos << " qual " << mr.quality <<" len "<< g_.length(e) );
            size_t cut = 0;
            size_t edge_start_pos = mr.mapped_range.start_pos;
            size_t edge_end_pos = mr.mapped_range.end_pos;
            size_t read_start_pos = mr.initial_range.start_pos + cut;
            size_t read_end_pos = mr.initial_range.end_pos;
            if (edge_start_pos >= edge_end_pos || read_start_pos >= read_end_pos) {
                DEBUG ("skipping extra-short alignment");
                continue;
            }
            res.insert(KmerCluster<Graph>(e, edge_start_pos, edge_end_pos, read_start_pos, read_end_pos, mr.quality));
        }
        DEBUG("Ended loading bwa")
        return res;
    }

    std::string DebugEmptyBestScoredPath(VertexId start_v, VertexId end_v, EdgeId prev_edge, EdgeId cur_edge,
                                         size_t prev_last_edge_position, size_t cur_first_edge_position, int seq_len) const {
        size_t result = GetDistance(start_v, end_v, /*update cache*/false);
        std::ostringstream ss;
        ss << "Tangled region between edges " << g_.int_id(prev_edge) << " " << g_.int_id(cur_edge) <<  " is not closed, additions from edges: "
            << int(g_.length(prev_edge)) - int(prev_last_edge_position) <<" " << int(cur_first_edge_position)
            << " and seq "<< seq_len << " and shortest path " << result;
        return ss.str();
    }

    void PrepareInitialState(omnigraph::MappingPath<debruijn_graph::EdgeId> &path, const Sequence &s, bool forward, Sequence &ss, EdgeId &start_e, int &start_pos, int &start_pos_seq) const {
        if (forward){
            start_e = path.edge_at(path.size() - 1);
            omnigraph::MappingRange mapping = path.mapping_at(path.size() - 1);
            start_pos = (int) mapping.mapped_range.end_pos;
            start_pos_seq = mapping.initial_range.end_pos;
            ss = s.Subseq(mapping.initial_range.end_pos, (int) s.size() );
            DEBUG("Forward e=" << start_e.int_id() << " sp=" << start_pos << " seq_sz" << ss.size())
        } else {
            start_e = g_.conjugate(path.edge_at(0));
            omnigraph::MappingRange mapping = path.mapping_at(0);
            start_pos = min((int) g_.length(start_e), (int) g_.length(start_e) + (int) g_.k() - (int) mapping.mapped_range.start_pos);
            start_pos_seq = 0;
            ss = !s.Subseq(0, mapping.initial_range.start_pos);
            DEBUG("Backward e=" << start_e.int_id() << " sp=" << start_pos << " seq_sz" << ss.size())
        }
    }

    void UpdatePath(vector<debruijn_graph::EdgeId> &path, std::vector<EdgeId> &ans, int end_pos, int end_pos_seq, PathRange &range, bool forward) const {
        if (forward) {
            if (end_pos < g_.k()) {
                ans.pop_back();
                end_pos = g_.length(ans[ans.size() - 1]);
            }
            for (int i = 1; i < (int) ans.size() - 1; ++i) {
                path.push_back(ans[i]);
            }
            path.push_back(ans[ans.size() - 1]);
            range.seq_end = end_pos_seq;
            range.edge_end = end_pos;
        } else {
            vector<debruijn_graph::EdgeId> cur_sorted;
            int start = (int) g_.length(ans[ans.size() - 1]) + (int) g_.k() - end_pos;
            int cur_ind = (int) ans.size() - 1;
            while (cur_ind >= 0 && start - (int) g_.length(ans[cur_ind]) > 0){
                start -= (int) g_.length(ans[cur_ind]);
                cur_ind --;
            }
            if (cur_ind > 0){
                cur_sorted.push_back(g_.conjugate(ans[cur_ind]) );
            }
            for (int i = cur_ind - 1; i > 0; --i) {
                cur_sorted.push_back(g_.conjugate(ans[i]));
            }
            for (size_t i = 0; i < path.size(); ++i) {
                cur_sorted.push_back(path[i]);
            }
            path = cur_sorted;
            range.seq_start = 0;
            range.edge_start = start;
        }
    }

    void GrowEnds(omnigraph::MappingPath<debruijn_graph::EdgeId> &bwa_hits, vector<debruijn_graph::EdgeId> &path, const Sequence &s, bool forward, PathRange &range, int &return_code) const {
        VERIFY(path.size() > 0);
        Sequence ss; 
        int start_pos = -1;
        int start_pos_seq = -1;
        return_code = 0;
        EdgeId start_e = EdgeId();
        PrepareInitialState(bwa_hits, s, forward, ss, start_e, start_pos, start_pos_seq);

        int s_len = int(ss.size());
        int score = max(10, s_len/5);
        if (s_len > (int) pb_config_.max_contigs_gap_length) {
            DEBUG("EdgeDijkstra: sequence is too long " << s_len)
            return_code += 1;
            return;
        }
        if (s_len < max((int) g_.length(start_e) + (int) g_.k() - start_pos, (int) g_.k())) {
            DEBUG("EdgeDijkstra: sequence is too small " << s_len)
            return_code += 2;
            return;
        }
        graph_aligner::DijkstraEndsReconstructor algo = graph_aligner::DijkstraEndsReconstructor(g_, gap_cfg_, ss.str(), start_e, start_pos, score);
        algo.CloseGap();
        score = algo.GetEditDistance();
        return_code += algo.GetReturnCode();
        if (score == std::numeric_limits<int>::max()){
            DEBUG("EdgeDijkstra didn't find anything edge=" << start_e.int_id() << " s_start=" << start_pos << " seq_len=" << ss.size())
            return;
        }
        std::vector<EdgeId> ans = algo.GetPath();
        int end_pos = algo.GetPathEndPosition();
        int end_pos_seq = forward? algo.GetSeqEndPosition() + start_pos_seq: 0;
        UpdatePath(path, ans, end_pos, end_pos_seq, range, forward);
    }

    void FillGapsInCluster(const vector<typename ClustersSet::iterator> &cur_cluster,
                           const Sequence &s, 
                           std::vector<vector<debruijn_graph::EdgeId> > &edges,
                           std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &bwa_hits) const {
        omnigraph::MappingPath<debruijn_graph::EdgeId> cur_sorted_hits;
        vector<debruijn_graph::EdgeId> cur_sorted_edges;
        EdgeId prev_edge = EdgeId();

        for (auto iter = cur_cluster.begin(); iter != cur_cluster.end();) {
            EdgeId cur_edge = (*iter)->edgeId;
            if (prev_edge != EdgeId()) {
//Need to find sequence of edges between clusters
                VertexId start_v = g_.EdgeEnd(prev_edge);
                VertexId end_v = g_.EdgeStart(cur_edge);
                auto prev_iter = iter - 1;
                MappingInstance cur_first_index = (*iter)->sorted_positions[(*iter)->first_trustable_index];
                MappingInstance prev_last_index = (*prev_iter)->sorted_positions[(*prev_iter)->last_trustable_index];
                double read_gap_len = (double) (cur_first_index.read_position - prev_last_index.read_position);
//FIXME:: is g_.k() relevant
                double stretched_graph_len = (prev_edge.int_id() != cur_edge.int_id()) ? 
                                            (double) (cur_first_index.edge_position + g_.k()) +
                                            ((int) g_.length(prev_edge) - prev_last_index.edge_position) * pb_config_.path_limit_stretching
                                            : ((int) cur_first_index.edge_position - (int) prev_last_index.edge_position) * pb_config_.path_limit_stretching;
                if ((start_v != end_v || (start_v == end_v && read_gap_len > stretched_graph_len)) 
                    && (prev_edge.int_id() != cur_edge.int_id() || (prev_edge.int_id() == cur_edge.int_id() && stretched_graph_len < 0) ||
                        (prev_edge.int_id() == cur_edge.int_id() && stretched_graph_len > 0 && read_gap_len > stretched_graph_len))) {
                    if (start_v == end_v) {
                        DEBUG("looking for path from vertex to itself, read pos"
                              << cur_first_index.read_position << " " << prev_last_index.read_position
                              << " edge pos: "<< cur_first_index.edge_position << " " << prev_last_index.edge_position
                              <<" edge len " << g_.length(prev_edge));
                        DEBUG("Read gap len " << read_gap_len << " stretched graph path len" <<  stretched_graph_len);
                    }
                    DEBUG(" traversing tangled hregion between "<< g_.int_id(prev_edge)<< " " << g_.int_id(cur_edge));
                    DEBUG(" first pair" << cur_first_index.str() << " edge_len" << g_.length(cur_edge));
                    DEBUG(" last pair" << prev_last_index.str() << " edge_len" << g_.length(prev_edge));
                    std::string s_add, e_add;
                    int seq_end = cur_first_index.read_position;
                    int seq_start = prev_last_index.read_position;
                    auto tmp = g_.EdgeNucls(prev_edge).str();
                    s_add = tmp.substr(prev_last_index.edge_position,
                                       g_.length(prev_edge) - prev_last_index.edge_position);
                    tmp = g_.EdgeNucls(cur_edge).str();
                    e_add = tmp.substr(0, cur_first_index.edge_position);
                    auto limits = GetPathLimits(**prev_iter, **iter, (int) s_add.length(), (int) e_add.length());
                    if (limits.first == -1) {
                        DEBUG ("Failed to find Path limits");
                        bwa_hits.push_back(cur_sorted_hits);
                        edges.push_back(cur_sorted_edges);
                        cur_sorted_edges.clear();
                        cur_sorted_hits.clear();
                        prev_edge = EdgeId();
                        continue;
                    }

                    int s_len = int(s.size());
                    int end_pos = seq_end;
                    if (seq_end < seq_start) {
                        DEBUG ("modifying limits because of some bullshit magic, seq length 0")
                        end_pos = seq_start;
                    }
                    DEBUG("taking subseq" << seq_start <<" "<< end_pos <<" " << s.size());
                    std::string seq_string = s.Subseq(seq_start, min(end_pos, s_len)).str();
                    graph_aligner::GapFiller gap_filler(g_, pb_config_, gap_cfg_);
                    graph_aligner::GapFillerResult res = gap_filler.Run(seq_string, graph_aligner::GraphPosition(prev_edge, prev_last_index.edge_position)
                                                                               , graph_aligner::GraphPosition(cur_edge, cur_first_index.edge_position)
                                                                               , limits.first, limits.second);
                    vector<EdgeId> intermediate_path = res.intermediate_path;
                    if (intermediate_path.size() == 0) {
                        DEBUG(DebugEmptyBestScoredPath(start_v, end_v, prev_edge, cur_edge,
                                      prev_last_index.edge_position, cur_first_index.edge_position, seq_end - seq_start));
                        bwa_hits.push_back(cur_sorted_hits);
                        edges.push_back(cur_sorted_edges);
                        cur_sorted_edges.clear();
                        cur_sorted_hits.clear();
                        prev_edge = EdgeId();
                        continue;
                    }
                    
                    for (EdgeId edge: intermediate_path) {
                        cur_sorted_edges.push_back(edge);
                    }
                }
            }
            MappingInstance cur_first_index = (*iter)->sorted_positions[(*iter)->first_trustable_index];
            MappingInstance cur_last_index = (*iter)->sorted_positions[(*iter)->last_trustable_index];
            cur_sorted_edges.push_back(cur_edge);
            cur_sorted_hits.push_back(cur_edge, omnigraph::MappingRange(Range(cur_first_index.read_position, cur_last_index.read_position), 
                                                                    Range(cur_first_index.edge_position, cur_last_index.edge_position) ));
            prev_edge = cur_edge;
            ++iter;
        }
        if (cur_sorted_edges.size() > 0) {
            edges.push_back(cur_sorted_edges);
            bwa_hits.push_back(cur_sorted_hits);
        }
    }

    bool TopologyGap(EdgeId first, EdgeId second, bool oriented) const {
        omnigraph::TerminalVertexCondition<Graph> condition(g_);
        bool res = condition.Check(g_.EdgeEnd(first)) && condition.Check(g_.EdgeStart(second));
        if (!oriented)
            res |= condition.Check(g_.EdgeStart(first)) && condition.Check(g_.EdgeEnd(second));
        return res;
    }

    std::vector<std::vector<bool>> FillConnectionsTable(const ClustersSet &mapping_descr) const {
        size_t len =  mapping_descr.size();
        TRACE("getting colors, table size "<< len);
        std::vector<std::vector<bool>> cons_table(len);
        for (size_t i = 0; i < len; i++) {
            cons_table[i].resize(len);
            cons_table[i][i] = 0;
        }
        size_t i = 0;
        for (auto i_iter = mapping_descr.begin(); i_iter != mapping_descr.end();
             ++i_iter, ++i) {
            size_t j = i;
            for (auto j_iter = i_iter;
                 j_iter != mapping_descr.end(); ++j_iter, ++j) {
                if (i_iter == j_iter)
                    continue;
                cons_table[i][j] = IsConsistent(*i_iter, *j_iter);
            }
        }
        return cons_table;
    }

    std::vector<int> GetWeightedColors(const ClustersSet &mapping_descr, int &num_colors) const {
        size_t len = mapping_descr.size();
        std::vector<int> colors(len, UNDEF_COLOR);
        std::vector<double> cluster_size(len);
        size_t ii = 0;
        for (const auto &cl : mapping_descr) {
            cluster_size[ii++] = cl.size * cl.quality;
        }
        const auto cons_table = FillConnectionsTable(mapping_descr);

        std::vector<double> max_size(len);
        std::vector<size_t> prev(len);

        int cur_color = 0;
        num_colors = 0;
        while (true) {
            for (size_t i = 0; i < len; ++i) {
                max_size[i] = 0;
                prev[i] = size_t(-1);
            }
            for (size_t i = 0; i < len; ++i) {
                if (colors[i] != UNDEF_COLOR) continue;
                max_size[i] = cluster_size[i];
                for (size_t j = 0; j < i; j++) {
                    if (colors[j] != -1) continue;
                    if (cons_table[j][i] && math::ls(max_size[i], cluster_size[i] + max_size[j])) {
                        max_size[i] = max_size[j] + cluster_size[i];
                        prev[i] = j;
                    }
                }
            }
            double maxx = 0;
            int maxi = -1;
            for (size_t j = 0; j < len; j++) {
                if (math::gr(max_size[j], maxx)) {
                    maxx = max_size[j];
                    maxi = int(j);
                }
            }
            if (maxi == -1) {
                break;
            }
            num_colors ++;
            cur_color = maxi;
            colors[maxi] = cur_color;
            int real_maxi = maxi, min_i = maxi;

            while (prev[maxi] != -1ul) {
                min_i = maxi;
                maxi = int(prev[maxi]);
                colors[maxi] = cur_color;
            }
            while (real_maxi >= min_i) {
                if (colors[real_maxi] == UNDEF_COLOR) {
                    colors[real_maxi] = DELETED_COLOR;
                }
                real_maxi--;
            }
        }
        DEBUG("Num hits clusters=" << num_colors);
        return colors;
    }

    OneReadMapping AddGapDescriptions(const std::vector<typename ClustersSet::iterator> &start_clusters,
                                      const std::vector<typename ClustersSet::iterator> &end_clusters,
                                      const std::vector<vector<debruijn_graph::EdgeId> > &sorted_edges,
                                      const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
                                      const std::vector<PathRange> &read_ranges,
                                      const Sequence &s,
                                      const std::vector<bool> &block_gap_closer) const {
        DEBUG("adding gaps between subreads");
        std::vector<GapDescription> illumina_gaps;
        for (size_t i = 0; i + 1 < sorted_edges.size() ; i++) {
            if (block_gap_closer[i])
                continue;
            size_t j = i + 1;
            EdgeId before_gap = sorted_edges[i][sorted_edges[i].size() - 1];
            EdgeId after_gap = sorted_edges[j][0];
//do not add "gap" for rc-jumping
            if (before_gap != after_gap && before_gap != g_.conjugate(after_gap)) {
                if (TopologyGap(before_gap, after_gap, true)) {
                    if (start_clusters[j]->CanFollow(*end_clusters[i])) {
                        const auto &a = *end_clusters[i];
                        const auto &b = *start_clusters[j];
                        size_t seq_start = a.sorted_positions[a.last_trustable_index].read_position + g_.k();
                        size_t seq_end = b.sorted_positions[b.first_trustable_index].read_position;
                        size_t left_offset = a.sorted_positions[a.last_trustable_index].edge_position;
                        size_t right_offset = b.sorted_positions[b.first_trustable_index].edge_position;
                        auto gap = CreateGapInfoTryFixOverlap(g_, s,
                                                    seq_start, seq_end,
                                                    a.edgeId, left_offset,
                                                    b.edgeId, right_offset);
                        if (gap != GapDescription()) {
                            illumina_gaps.push_back(gap);
                            DEBUG("adding gap between alignments number " << i<< " and " << j);
                        }
                    }

                }
            }
        }
        DEBUG("Resulting hits num=" << sorted_edges.size());
        return OneReadMapping(sorted_edges, sorted_bwa_hits, illumina_gaps, read_ranges);
    }

    void ProcessCluster(const Sequence &s,
                        std::vector<typename ClustersSet::iterator> &cur_cluster,
                        std::vector<typename ClustersSet::iterator> &start_clusters,
                        std::vector<typename ClustersSet::iterator> &end_clusters,
                        std::vector<vector<debruijn_graph::EdgeId> > &sorted_edges,
                        std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
                        std::vector<bool> &block_gap_closer) const {
        std::sort(cur_cluster.begin(), cur_cluster.end(),
                  [](const typename ClustersSet::iterator& a, const typename ClustersSet::iterator& b) {
                      return (a->average_read_position < b->average_read_position);
                  });
        VERIFY(cur_cluster.size() > 0);
        auto cur_cluster_start = cur_cluster.begin();
        for (auto iter = cur_cluster.begin(); iter != cur_cluster.end(); ++iter) {
            auto next_iter = iter + 1;
            if (next_iter == cur_cluster.end() || !IsConsistent(**iter, **next_iter)) {
                if (next_iter != cur_cluster.end()) {
                    DEBUG("clusters splitted:");
                    DEBUG("on "<< (*iter)->str(g_));
                    DEBUG("and " << (*next_iter)->str(g_));
                }
                vector<typename ClustersSet::iterator> splitted_cluster(cur_cluster_start, next_iter);
                std::vector<vector<debruijn_graph::EdgeId> > edges;
                std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > bwa_hits;
                FillGapsInCluster(splitted_cluster, s, edges, bwa_hits);
                for (auto &cur_sorted: edges) {
                    DEBUG("Adding " <<edges.size() << " subreads, cur alignments " << cur_sorted.size());
                    if (cur_sorted.size() > 0) {
                        for(EdgeId eee: cur_sorted) {
                            DEBUG (g_.int_id(eee));
                        }
                        start_clusters.push_back(*cur_cluster_start);
                        end_clusters.push_back(*iter);
                        sorted_edges.push_back(cur_sorted);
                        //Blocking gap closing inside clusters;
                        block_gap_closer.push_back(true);
                    }
                }
                for (auto &cur_sorted: bwa_hits) {
                    if (cur_sorted.size() > 0) {
                        sorted_bwa_hits.push_back(cur_sorted);
                    }
                }
                if (block_gap_closer.size() > 0)
                    block_gap_closer[block_gap_closer.size() - 1] = false;
                cur_cluster_start = next_iter;
            } else {
                DEBUG("connected consecutive clusters:");
                DEBUG("on "<< (*iter)->str(g_));
                DEBUG("and " << (*next_iter)->str(g_));
            }
        }
    }

    OneReadMapping GetReadAlignment(const io::SingleRead &read) const {
        Sequence s = read.sequence();
        ClustersSet mapping_descr = GetBWAClusters(read); //GetOrderClusters(s);
        int num_colors = 0;
        vector<int> colors = GetWeightedColors(mapping_descr, num_colors);
        size_t len =  mapping_descr.size();
        std::vector<vector<debruijn_graph::EdgeId> > sorted_edges;
        vector<omnigraph::MappingPath<debruijn_graph::EdgeId>> sorted_bwa_hits;
        omnigraph::MappingRange read_range;
        vector<bool> block_gap_closer;
        vector<typename ClustersSet::iterator> start_clusters, end_clusters;
        vector<int> used(len);
        auto iter = mapping_descr.begin();
        for (size_t i = 0; i < len; i++, iter ++) {
            used[i] = 0;
            DEBUG(colors[i] <<" " << iter->str(g_));
        }
        //FIXME ferther code is AWFUL
        for (size_t i = 0; i < len; i++) {
            int cur_color = colors[i];
            if (!used[i] && cur_color != DELETED_COLOR) {
                DEBUG("starting new subread");
                std::vector<typename ClustersSet::iterator> cur_cluster;
                used[i] = 1;
                int j = 0;
                for (auto i_iter = mapping_descr.begin(); i_iter != mapping_descr.end(); ++i_iter, ++j) {
                    if (colors[j] == cur_color) {
                        cur_cluster.push_back(i_iter);
                        used[j] = 1;
                    }
                }
                ProcessCluster(s, cur_cluster, start_clusters, end_clusters, sorted_edges, sorted_bwa_hits, block_gap_closer);
            }
        }
        std::vector<PathRange> read_ranges;
        if (sorted_edges.size() == 1 && gap_cfg_.restore_ends) {
            bool forward = true;
            int return_code = 0;
            PathRange cur_range;
            cur_range.seq_start = sorted_bwa_hits[0].mapping_at(0).initial_range.start_pos;
            cur_range.seq_end = sorted_bwa_hits[0].mapping_at(sorted_bwa_hits[0].size() - 1).initial_range.end_pos;
            cur_range.edge_start = sorted_bwa_hits[0].mapping_at(0).mapped_range.start_pos;
            cur_range.edge_end = sorted_bwa_hits[0].mapping_at(sorted_bwa_hits[0].size() - 1).mapped_range.end_pos;
            GrowEnds(sorted_bwa_hits[0], sorted_edges[0], s, !forward, cur_range, return_code);
            DEBUG("Backward return_code_ends=" << return_code)
            GrowEnds(sorted_bwa_hits[0], sorted_edges[0], s, forward, cur_range, return_code);
            DEBUG("Forward return_code_ends=" << return_code)
            read_ranges.push_back(cur_range);
        } else {
            for (auto hits: sorted_bwa_hits){
                PathRange cur_range;
                cur_range.seq_start = hits.mapping_at(0).initial_range.start_pos;
                cur_range.seq_end = hits.mapping_at(hits.size() - 1).initial_range.end_pos;
                cur_range.edge_start = hits.mapping_at(0).mapped_range.start_pos;
                cur_range.edge_end = hits.mapping_at(hits.size() - 1).mapped_range.end_pos;
                read_ranges.push_back(cur_range);
            }
        }
        return AddGapDescriptions(start_clusters, end_clusters, sorted_edges, sorted_bwa_hits, read_ranges, s, block_gap_closer);
    }

    std::pair<int, int> GetPathLimits(const KmerCluster<Graph> &a,
                                      const KmerCluster<Graph> &b,
                                      int s_add_len, int e_add_len) const {
        int start_pos = a.sorted_positions[a.last_trustable_index].read_position;
        int end_pos = b.sorted_positions[b.first_trustable_index].read_position;
        int seq_len = -start_pos + end_pos;
        //int new_seq_len =
//TODO::something more reasonable
        int path_min_len = std::max(int(floor((seq_len - int(g_.k())) * pb_config_.path_limit_pressing)), 0);
        int path_max_len = (int) ((double) (seq_len + g_.k() * 2) * pb_config_.path_limit_stretching);
        if (seq_len < 0) {
            DEBUG("suspicious negative seq_len " << start_pos << " " << end_pos << " " << path_min_len << " " << path_max_len);
            if (path_max_len < 0)
            return std::make_pair(-1, -1);
        }
        path_min_len = std::max(path_min_len - int(s_add_len + e_add_len), 0);
        path_max_len = std::max(path_max_len - int(s_add_len + e_add_len), 0);
        return std::make_pair(path_min_len, path_max_len);
    }

    size_t GetDistance(VertexId start_v, VertexId end_v,
                       bool update_cache = true) const {
        size_t result = size_t(-1);
        auto vertex_pair = std::make_pair(start_v, end_v);
        bool not_found;
        auto distance_it = distance_cashed_.begin();
        //FIXME should permit multiple readers
#pragma omp critical(pac_index)
        {
            distance_it = distance_cashed_.find(vertex_pair);
            not_found = (distance_it == distance_cashed_.end());
        }
        if (not_found) {
            omnigraph::DijkstraHelper<debruijn_graph::Graph>::BoundedDijkstra dijkstra(
                    omnigraph::DijkstraHelper<debruijn_graph::Graph>::CreateBoundedDijkstra(g_,
                                                                                            pb_config_.max_path_in_dijkstra,
                                                                                            pb_config_.max_vertex_in_dijkstra));
            dijkstra.Run(start_v);
            if (dijkstra.DistanceCounted(end_v)) {
                result = dijkstra.GetDistance(end_v);
            }
            if (update_cache) {
#pragma omp critical(pac_index)
                distance_cashed_.insert({vertex_pair, result});
            }
        } else {
            TRACE("taking from cashed");
            result = distance_it->second;
        }

        return result;
    }

    bool IsConsistent(const KmerCluster<Graph> &a,
                      const KmerCluster<Graph> &b) const {
        EdgeId a_edge = a.edgeId;
        EdgeId b_edge = b.edgeId;
        DEBUG("Checking consistency: " << g_.int_id(a_edge) << " and " << g_.int_id(b_edge));
        //FIXME: Is this check useful?
        if (a.sorted_positions[a.last_trustable_index].read_position +
                    (int) pb_config_.max_path_in_dijkstra <
            b.sorted_positions[b.first_trustable_index].read_position) {
            DEBUG("Clusters are too far in read");
            return false;
        }
        size_t result = GetDistance(g_.EdgeEnd(a_edge), g_.EdgeStart(b_edge));
        DEBUG ("Distance: " << result);
        if (result == size_t(-1)) {
            return false;
        }
        result += g_.length(a_edge);
        if (similar_in_graph(a.sorted_positions[1], b.sorted_positions[0], (int)result)) {
            DEBUG(" similar! ")
            return true;
        } else {
//FIXME: reconsider this condition! i.e only one large range? That may allow to decrease the bwa length cutoff
            if (- a.sorted_positions[1].edge_position +  (int)result + b.sorted_positions[0].edge_position  <=
                        b.sorted_positions[0].read_position - a.sorted_positions[1].read_position + 2 * int(g_.k())
                        && (pb_config_.bwa_length_cutoff > 2*g_.k() || pb_config_.internal_length_cutoff > 2*g_.k()) ) {
                DEBUG("overlapping range magic worked, " << - a.sorted_positions[1].edge_position +
                                                           (int)result + b.sorted_positions[0].edge_position
                      << " and " <<  b.sorted_positions[0].read_position - a.sorted_positions[1].read_position + 2 * g_.k());
                DEBUG("Ranges:" << a.str(g_) << " " << b.str(g_)
                                <<" llength and dijkstra shift :" << g_.length(a_edge) << " " << (result - g_.length(a_edge)));
                return true;
            } else {
                DEBUG("Not similar")
                return false;
            }
        }
    }

};
}
