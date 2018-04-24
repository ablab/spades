//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/index/edge_multi_index.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/paths/path_processor.hpp"

#include "modules/alignment/edge_index_refiller.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "modules/alignment/gap_info.hpp"

#include "pipeline/configs/aligner_config.hpp"

#include "sequence/sequence_tools.hpp"

#include "pacbio_read_structures.hpp"

#include <algorithm>
#include <vector>
#include <set>

namespace pacbio {
enum {
    UNDEF_COLOR = -1,
    DELETED_COLOR = -2
};

struct OneReadMapping {
    std::vector<std::vector<debruijn_graph::EdgeId>> main_storage;
    std::vector<GapDescription> gaps;
    OneReadMapping(const std::vector<std::vector<debruijn_graph::EdgeId>>& main_storage_,
                   const std::vector<GapDescription>& gaps_)
            :  main_storage(main_storage_), gaps(gaps_) {}
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

public:

    PacBioMappingIndex(const Graph &g,
                       debruijn_graph::config::pacbio_processor pb_config, alignment::BWAIndex::AlignmentMode mode)
            : g_(g),
              pb_config_(pb_config),
              bwa_mapper_(g, mode, pb_config.bwa_length_cutoff) {
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
            return ((b.edge_position + shift - a.edge_position) * pb_config_.compression_cutoff <= (b.read_position - a.read_position));
        }
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
        if (res.size() != mapped_path.size()) {
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

    ClustersSet GetBWAClusters(const Sequence &s) const {
        DEBUG("BWA started")
        ClustersSet res;
        if (s.size() < g_.k())
            return res;

        omnigraph::MappingPath<EdgeId> mapped_path = FilterSpuriousAlignments(bwa_mapper_.MapSequence(s), s.size());
        TRACE(read_count_ << " read_count_");
        TRACE("BWA ended")
        DEBUG(mapped_path.size() <<"  clusters");
        for (const auto &e_mr : mapped_path) {
            EdgeId e = e_mr.first;
            omnigraph::MappingRange mr = e_mr.second;
            DEBUG("BWA loading edge=" << g_.int_id(e) << " e_start=" << mr.mapped_range.start_pos << " e_end=" << mr.mapped_range.end_pos 
                                                                    << " r_start=" << mr.initial_range.start_pos << " r_end=" << mr.initial_range.end_pos );
            size_t cut = 0;
            size_t edge_start_pos = mr.mapped_range.start_pos;
            size_t edge_end_pos = mr.mapped_range.end_pos;
            size_t read_start_pos = mr.initial_range.start_pos + cut;
            size_t read_end_pos = mr.initial_range.end_pos;
            if (edge_start_pos >= edge_end_pos || read_start_pos >= read_end_pos) {
                DEBUG ("skipping extra-short alignment");
                continue;
            }
            res.insert(KmerCluster<Graph>(e, edge_start_pos, edge_end_pos, read_start_pos, read_end_pos));

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

    vector<vector<EdgeId>> FillGapsInCluster(const vector<typename ClustersSet::iterator> &cur_cluster,
                                     const Sequence &s) const {
        vector<EdgeId> cur_sorted;
        vector<vector<EdgeId>> res;
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
                double stretched_graph_len = (double) (cur_first_index.edge_position + g_.k()) +
                        ((int) g_.length(prev_edge) - prev_last_index.edge_position) * pb_config_.path_limit_stretching;
                if (start_v != end_v || (start_v == end_v && read_gap_len > stretched_graph_len)) {
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
                    string s_add = "";
                    string e_add = "";
                    int seq_end = cur_first_index.read_position;
                    int seq_start = prev_last_index.read_position;
                    string tmp = g_.EdgeNucls(prev_edge).str();
                    s_add = tmp.substr(prev_last_index.edge_position,
                                       g_.length(prev_edge) - prev_last_index.edge_position);
                    tmp = g_.EdgeNucls(cur_edge).str();
                    e_add = tmp.substr(0, cur_first_index.edge_position);
                    pair<int, int> limits = GetPathLimits(**prev_iter, **iter, (int) s_add.length(), (int) e_add.length());
                    if (limits.first == -1) {
                        DEBUG ("Failed to find Path limits");
                        res.push_back(cur_sorted);
                        cur_sorted.clear();
                        prev_edge = EdgeId();
                        continue;
                    }
                    
                    vector<EdgeId> intermediate_path = BestScoredPath(s, start_v, end_v, limits.first, limits.second, seq_start, seq_end, s_add, e_add);
                    if (intermediate_path.size() == 0) {
                        DEBUG(DebugEmptyBestScoredPath(start_v, end_v, prev_edge, cur_edge,
                                      prev_last_index.edge_position, cur_first_index.edge_position, seq_end - seq_start));
                        res.push_back(cur_sorted);
                        cur_sorted.clear();
                        prev_edge = EdgeId();
                        continue;
                    }
                    for (EdgeId edge: intermediate_path) {
                        cur_sorted.push_back(edge);
                    }
                }
            }
            cur_sorted.push_back(cur_edge);
            prev_edge = cur_edge;
            ++iter;
        }
        if (cur_sorted.size() > 0)
            res.push_back(cur_sorted);
        return res;
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

    std::vector<int> GetWeightedColors(const ClustersSet &mapping_descr) const {
        size_t len = mapping_descr.size();
        std::vector<int> colors(len, UNDEF_COLOR);
        std::vector<int> cluster_size(len);
        size_t ii = 0;
        for (const auto &cl : mapping_descr) {
            cluster_size[ii++] = cl.size;
        }
        const auto cons_table = FillConnectionsTable(mapping_descr);

        std::vector<int> max_size(len);
        std::vector<size_t> prev(len);

        int cur_color = 0;
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
                    if (cons_table[j][i] && max_size[i] < cluster_size[i] + max_size[j]) {
                        max_size[i] = max_size[j] + cluster_size[i];
                        prev[i] = j;
                    }
                }
            }
            int maxx = 0;
            int maxi = -1;
            for (size_t j = 0; j < len; j++) {
                if (max_size[j] > maxx) {
                    maxx = max_size[j];
                    maxi = int(j);
                }
            }
            if (maxi == -1) {
                break;
            }
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
        return colors;
    }

    OneReadMapping AddGapDescriptions(const vector<typename ClustersSet::iterator> &start_clusters,
                                      const vector<typename ClustersSet::iterator> &end_clusters,
                                      const std::vector<vector<EdgeId>> &sorted_edges, const Sequence &s,
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
        return OneReadMapping(sorted_edges, illumina_gaps);
    }

    void ProcessCluster(const Sequence &s,
                        std::vector<typename ClustersSet::iterator> &cur_cluster,
                        std::vector<typename ClustersSet::iterator> &start_clusters,
                        std::vector<typename ClustersSet::iterator> &end_clusters,
                        std::vector<std::vector<EdgeId>> &sorted_edges,
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
                auto res = FillGapsInCluster(splitted_cluster, s);
                for (auto &cur_sorted:res) {
                    DEBUG("Adding " <<res.size() << " subreads, cur alignments " << cur_sorted.size());
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
        ClustersSet mapping_descr = GetBWAClusters(s); //GetOrderClusters(s);
        vector<int> colors = GetWeightedColors(mapping_descr);
        size_t len =  mapping_descr.size();
        vector<vector<EdgeId>> sorted_edges;
        vector<bool> block_gap_closer;
        vector<typename ClustersSet::iterator> start_clusters, end_clusters;
        vector<int> used(len);
        auto iter = mapping_descr.begin();
        for (size_t i = 0; i < len; i++, iter ++) {
            used[i] = 0;
            DEBUG(colors[i] <<" " << iter->str(g_));
        }
        for (size_t i = 0; i < len; i++) {
            int cur_color = colors[i];
            if (!used[i] && cur_color != DELETED_COLOR) {
                DEBUG("starting new subread");
                vector<typename ClustersSet::iterator> cur_cluster;
                used[i] = 1;
                int j = 0;
                for (auto i_iter = mapping_descr.begin(); i_iter != mapping_descr.end(); ++i_iter, ++j) {
                    if (colors[j] == cur_color) {
                        cur_cluster.push_back(i_iter);
                        used[j] = 1;
                    }
                }
                ProcessCluster(s, cur_cluster, start_clusters, end_clusters, sorted_edges, block_gap_closer);
            }
        }
        return AddGapDescriptions(start_clusters, end_clusters, sorted_edges, s, block_gap_closer);
    }

    std::pair<int, int> GetPathLimits(const KmerCluster<Graph> &a,
                                      const KmerCluster<Graph> &b,
                                      int s_add_len, int e_add_len) const {
        int start_pos = a.sorted_positions[a.last_trustable_index].read_position;
        int end_pos = b.sorted_positions[b.first_trustable_index].read_position;
        int seq_len = -start_pos + end_pos;
        //int new_seq_len =
//TODO::something more reasonable
        int path_min_len = max(int(floor((seq_len - int(g_.k())) * pb_config_.path_limit_pressing)), 0);
        int path_max_len = (int) ((double) (seq_len + g_.k() * 2) * pb_config_.path_limit_stretching);
        if (seq_len < 0) {
            DEBUG("suspicious negative seq_len " << start_pos << " " << end_pos << " " << path_min_len << " " << path_max_len);
            if (path_max_len < 0)
            return std::make_pair(-1, -1);
        }
        path_min_len = max(path_min_len - int(s_add_len + e_add_len), 0);
        path_max_len = max(path_max_len - int(s_add_len + e_add_len), 0);
        return std::make_pair(path_min_len, path_max_len);
    }

    size_t GetDistance(VertexId start_v, VertexId end_v,
                       bool update_cache = true) const {
        size_t result = size_t(-1);
        pair<VertexId, VertexId> vertex_pair = make_pair(start_v, end_v);
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
        DEBUG("clusters on " << g_.int_id(a_edge) << " and " << g_.int_id(b_edge));
        //FIXME: Is this check useful?
        if (a.sorted_positions[a.last_trustable_index].read_position +
                    (int) pb_config_.max_path_in_dijkstra <
            b.sorted_positions[b.first_trustable_index].read_position) {
            DEBUG("Clusters are too far in read");
            return false;
        }
        size_t result = GetDistance(g_.EdgeEnd(a_edge), g_.EdgeStart(b_edge));
        DEBUG (result);
        if (result == size_t(-1)) {
            return false;
        }
        result += g_.length(a_edge);
        if (similar_in_graph(a.sorted_positions[1], b.sorted_positions[0], (int)result)) {
            DEBUG(" similar! ")
        } else {
//FIXME: reconsider this condition! i.e only one large range? That may allow to decrease the bwa length cutoff
            if (a.size > LONG_ALIGNMENT_OVERLAP && b.size > LONG_ALIGNMENT_OVERLAP &&
                 - a.sorted_positions[1].edge_position +  (int)result + b.sorted_positions[0].edge_position  <=
                        b.sorted_positions[0].read_position - a.sorted_positions[1].read_position + 2 * int(g_.k())) {
                DEBUG("overlapping range magic worked, " << - a.sorted_positions[1].edge_position +
                                                           (int)result + b.sorted_positions[0].edge_position
                      << " and " <<  b.sorted_positions[0].read_position - a.sorted_positions[1].read_position + 2 * g_.k());
                DEBUG("Ranges:" << a.str(g_) << " " << b.str(g_)
                                <<" llength and dijkstra shift :" << g_.length(a_edge) << " " << (result - g_.length(a_edge)));
            } else {
                return false;
            }
        }
        return true;
    }

    std::string PathToString(const vector<EdgeId>& path) const {
        std::string res = "";
        for (auto iter = path.begin(); iter != path.end(); iter++) {
            size_t len = g_.length(*iter);
            std::string tmp = g_.EdgeNucls(*iter).First(len).str();
            res = res + tmp;
        }
        return res;
    }

    std::vector<EdgeId> BestScoredPath(const Sequence &s, VertexId start_v, VertexId end_v,
                                       int path_min_length, int path_max_length,
                                       int start_pos, int end_pos,
                                       std::string &s_add, std::string &e_add) const {
        TRACE(" Traversing tangled region. Start and end vertices resp: " << g_.int_id(start_v) <<" " << g_.int_id(end_v));
        omnigraph::PathStorageCallback<Graph> callback(g_);
        int pres = ProcessPaths(g_,
                                path_min_length, path_max_length,
                                start_v, end_v,
                                callback);
        DEBUG("PathProcessor result: " << pres << " limits " << path_min_length << " " << path_max_length);
        std::vector<std::vector<EdgeId> > paths = callback.paths();
        TRACE("taking subseq" << start_pos <<" "<< end_pos <<" " << s.size());
        int s_len = int(s.size());
        if (end_pos < start_pos) {
            DEBUG ("modifying limits because of some bullshit magic, seq length 0")
            end_pos = start_pos;
        }
        std::string seq_string = s.Subseq(start_pos, min(end_pos + 1, s_len)).str();
        size_t best_path_ind = paths.size();
        int best_score = std::numeric_limits<int>::max();
        if (paths.size() == 0) {
            DEBUG("need to find best scored path between "<<paths.size()<<" , seq_len " << seq_string.length());
            DEBUG ("no paths");
            return vector<EdgeId>(0);
        }
        if (seq_string.length() > pb_config_.max_contigs_gap_length) {
            DEBUG("need to find best scored path between "<<paths.size()<<" , seq_len " << seq_string.length());
            DEBUG("Gap is too large");
            return vector<EdgeId>(0);
        }
        bool additional_debug = (paths.size() > 1 && paths.size() < 10);
        for (size_t i = 0; i < paths.size(); i++) {
            DEBUG("path len "<< paths[i].size());
            if (paths[i].size() == 0) {
                DEBUG ("Pathprocessor returns path with size = 0")
            }
            std::string cur_string = s_add + PathToString(paths[i]) + e_add;
            TRACE("cur_string: " << cur_string <<"\n seq_string " << seq_string);
            int cur_score = StringDistance(cur_string, seq_string);
            //DEBUG only
            if (additional_debug) {
                TRACE("candidate path number "<< i << " , len " << cur_string.length());
                TRACE("graph candidate: " << cur_string);
                TRACE("in pacbio read: " << seq_string);
                for (auto j_iter = paths[i].begin(); j_iter != paths[i].end();
                        ++j_iter) {
                    DEBUG(g_.int_id(*j_iter));
                }
                DEBUG("score: "<< cur_score);
            }
            if (cur_score < best_score) {
                best_score = cur_score;
                best_path_ind = i;
            }
        }
        TRACE(best_score);
        if (best_score == std::numeric_limits<int>::max()) {
            if (paths.size() < 10) {
                for (size_t i = 0; i < paths.size(); i++) {
                    DEBUG ("failed with strings " << seq_string << " " << s_add + PathToString(paths[i]) + e_add);
                }
            }
            DEBUG (paths.size() << " paths available");
            return std::vector<EdgeId>(0);
        }
        if (additional_debug) {
            TRACE("best score found! Path " <<best_path_ind <<" score "<< best_score);
        }
        return paths[best_path_ind];
    }
};

}
