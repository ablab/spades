//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/index/edge_multi_index.hpp"
#include "modules/alignment/edge_index_refiller.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/paths/path_processor.hpp"
// FIXME: Layering violation, get rid of this
#include "pipeline/config_struct.hpp"
#include "pacbio_read_structures.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"

#include <algorithm>
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"

namespace pacbio {
enum {
    UNDEF_COLOR = -1,
    DELETED_COLOR = - 2
};

struct OneReadMapping {
    vector<vector<debruijn_graph::EdgeId>> main_storage;
    vector<GapDescription> gaps;
    vector<size_t> real_length;
    OneReadMapping(const vector<vector<debruijn_graph::EdgeId>>& main_storage_,
                   const vector<GapDescription>& gaps_,
                   const vector<size_t>& real_length_) :
            main_storage(main_storage_), gaps(gaps_), real_length(real_length_) {
    }

};

template<class Graph>
class PacBioMappingIndex {
public:
    typedef map<typename Graph::EdgeId, vector<MappingInstance> > MappingDescription;
    typedef pair<typename Graph::EdgeId, vector<MappingInstance> > ClusterDescription;
    typedef set<KmerCluster<Graph> > ClustersSet;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef debruijn_graph::DeBruijnEdgeMultiIndex<typename Graph::EdgeId> Index;
    typedef typename Index::KeyWithHash KeyWithHash;

private:
    DECL_LOGGER("PacIndex")

    const Graph &g_;
    //size_t pacbio_k;
    size_t debruijn_k;
    const static int short_edge_cutoff = 0;
    const static size_t min_cluster_size = 8;
    const static int max_similarity_distance = 500;

//Debug stasts
    int good_follow = 0;
    int half_bad_follow = 0;
    int bad_follow = 0;

    //set<Sequence> banned_kmers;
    //debruijn_graph::DeBruijnEdgeMultiIndex<typename Graph::EdgeId> tmp_index;
    mutable map<pair<VertexId, VertexId>, size_t> distance_cashed;
    size_t read_count;
    bool ignore_map_to_middle;
    debruijn_graph::config::debruijn_config::pacbio_processor pb_config_;

    alignment::BWAReadMapper<Graph> bwa_mapper_;

public:
    MappingDescription GetSeedsFromRead(const Sequence &s) const;


    PacBioMappingIndex(const Graph &g, size_t k, size_t debruijn_k_, bool ignore_map_to_middle, string out_dir, 
                        debruijn_graph::config::debruijn_config::pacbio_processor pb_config, alignment::BWAIndex::AlignmentMode mode)
            : g_(g),
        //      pacbio_k(k),
              debruijn_k(debruijn_k_),
        //      tmp_index((unsigned) pacbio_k, out_dir), 
              ignore_map_to_middle(ignore_map_to_middle), 
              pb_config_(pb_config),
              bwa_mapper_(g, mode){
              //minimap_mapper(g, out_dir){
        DEBUG("PB Mapping Index construction started");
        //debruijn_graph::EdgeIndexRefiller().Refill(tmp_index, g_);
        DEBUG("Index constructed");
        //FillBannedKmers();
        read_count = 0;
    }
    ~PacBioMappingIndex(){
        DEBUG("good/ugly/bad counts:" << good_follow << " "<<half_bad_follow << " " << bad_follow);
    }
    
    //void FillBannedKmers() {
    //    for (int i = 0; i < 4; i++) {
    //        auto base = nucl((unsigned char) i);
    //        for (int j = 0; j < 4; j++) {
    //            auto other = nucl((unsigned char) j);
    //            for (size_t other_pos = 0; other_pos < pacbio_k; other_pos++) {
    //                string s = "";
    //                for (size_t k = 0; k < pacbio_k; k++) {
    //                    if (k != other_pos)
    //                        s += base;
    //                    else
    //                        s += other;
    //                }
    //                banned_kmers.insert(Sequence(s));
    //            }
    //        }
    //    }
    //}

    bool similar(const MappingInstance &a, const MappingInstance &b,
                        int shift = 0) const {
        if (b.read_position + shift < a.read_position) {
            return similar(b, a, -shift);
        } else if (b.read_position == a.read_position) {
            return (abs(int(b.edge_position) + shift - int(a.edge_position)) < 2);
        } else {
            return ((b.edge_position + shift - a.edge_position >= (b.read_position - a.read_position) * pb_config_.compression_cutoff) &&
                ((b.edge_position + shift - a.edge_position) * pb_config_.compression_cutoff <= (b.read_position - a.read_position)));
        }
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

    typename omnigraph::MappingPath<EdgeId> FilterSpuriousAlignments(typename omnigraph::MappingPath<EdgeId> mp_path, size_t seq_len) const {
        omnigraph::MappingPath<EdgeId> res;
        vector<pair<EdgeId,typename omnigraph::MappingRange> > mapped_path;
//FIXME do we need this now?
        for (size_t i = 0; i < mp_path.size(); i++) 
            mapped_path.push_back(make_pair(mp_path[i].first, mp_path[i].second));
        std::sort(mapped_path.begin(), mapped_path.end(),
                     [](const std::pair<EdgeId, typename omnigraph::MappingRange> &a, const std::pair<EdgeId, typename omnigraph::MappingRange> &b) {
            return (a.second.initial_range.end_pos+ a.second.initial_range.start_pos < b.second.initial_range.end_pos+ b.second.initial_range.start_pos); 
        });

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
            if (rlen < 500 && (rlen + g_.k()) * 2 < expected_additional_left + expected_additional_right) {
                INFO ("Skipping spurious alignment " << i << " on edge " << mapped_path[i].first);
            } else
                res.push_back(mapped_path[i].first, mapped_path[i].second);


        }
        if (res.size() != mapped_path.size()) {
            INFO("Seq len " << seq_len);
            for (const auto &e_mr : mapped_path) {
                EdgeId e = e_mr.first;
                omnigraph::MappingRange mr = e_mr.second;
                INFO("Alignents were" << g_.int_id(e) << " e_start=" << mr.mapped_range.start_pos << " e_end=" <<
                     mr.mapped_range.end_pos
                     << " r_start=" << mr.initial_range.start_pos << " r_end=" << mr.initial_range.end_pos);
            }
        }
        return res;
    }

    ClustersSet GetBWAClusters(const Sequence &s) const {
        DEBUG("BWA started")
        ClustersSet res;
        if (s.size() < g_.k()){
            return res;
        }
        omnigraph::MappingPath<EdgeId> mapped_path = FilterSpuriousAlignments(bwa_mapper_.MapSequence(s), s.size());
        TRACE(read_count << " read_count");
        TRACE("BWA ended")
        DEBUG(mapped_path.size() <<"  clusters");
        for (const auto &e_mr : mapped_path) {
            EdgeId e = e_mr.first;
            omnigraph::MappingRange mr = e_mr.second;
            DEBUG("BWA loading edge=" << g_.int_id(e) << " e_start=" << mr.mapped_range.start_pos << " e_end=" << mr.mapped_range.end_pos 
                                                                    << " r_start=" << mr.initial_range.start_pos << " r_end=" << mr.initial_range.end_pos );
            size_t cut = 0;
            size_t edge_start_pos = mr.mapped_range.start_pos;
/*          if (edge_start_pos < g_.k()) {
                cut = g_.k() - edge_start_pos;
                edge_start_pos = g_.k();
            }
*/
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

    vector<vector<EdgeId>> FillGapsInCluster(const vector<typename ClustersSet::iterator> &cur_cluster,
                                     const Sequence &s) const {
        vector<EdgeId> cur_sorted;
        vector<vector<EdgeId>> res;
        EdgeId prev_edge = EdgeId(0);

        for (auto iter = cur_cluster.begin(); iter != cur_cluster.end();) {
            EdgeId cur_edge = (*iter)->edgeId;
            if (prev_edge != EdgeId(0)) {
//Need to find sequence of edges between clusters
                VertexId start_v = g_.EdgeEnd(prev_edge);
                VertexId end_v = g_.EdgeStart(cur_edge);
                auto prev_iter = iter - 1;
                MappingInstance cur_first_index = (*iter)->sorted_positions[(*iter)->first_trustable_index];
                MappingInstance prev_last_index = (*prev_iter)->sorted_positions[(*prev_iter)->last_trustable_index];

                if (start_v != end_v ||
                        (start_v == end_v &&
                     (double) (cur_first_index.read_position - prev_last_index.read_position) >
//FIXME:: is g_.k() relevant
                     (double) (cur_first_index.edge_position + (int) g_.length(prev_edge) - prev_last_index.edge_position) * 1.3 + g_.k() )) {
                    if (start_v == end_v) {
                        DEBUG("looking for path from vertex to itself, read pos"
                              << cur_first_index.read_position << " " << prev_last_index.read_position
                              << " edge pos: "<< cur_first_index.edge_position << " " << prev_last_index.edge_position
                              <<" edge len " << g_.length(prev_edge));
                        DEBUG((double) (cur_first_index.read_position - prev_last_index.read_position) << " " << (double) 
                        (cur_first_index.edge_position + (int) g_.length(prev_edge) - prev_last_index.edge_position) * 1.3 + g_.k());
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
                        DEBUG ("Failed to find Path limints");
                        res.push_back(cur_sorted);
                        cur_sorted.clear();
                        prev_edge = EdgeId(0);
                        continue;
                    }

                    vector<EdgeId> intermediate_path = BestScoredPath(s, start_v, end_v, limits.first, limits.second, seq_start, seq_end, s_add, e_add);
                    if (intermediate_path.size() == 0) {
                        DEBUG("Tangled region between edges "<< g_.int_id(prev_edge) << " " << g_.int_id(cur_edge) << " is not closed, additions from edges: " << int(g_.length(prev_edge)) - int(prev_last_index.edge_position) <<" " << int(cur_first_index.edge_position)  << " and seq "<< - seq_start + seq_end);
                        res.push_back(cur_sorted);
                        cur_sorted.clear();
                        prev_edge = EdgeId(0);
                        continue;
                    }
                    for (auto j_iter = intermediate_path.begin(); j_iter != intermediate_path.end(); j_iter++) {
                        cur_sorted.push_back(*j_iter);
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

    vector<vector<int>> FillConnectionsTable (const ClustersSet &mapping_descr) const{
        size_t len =  mapping_descr.size();
        TRACE("getting colors, table size "<< len);
        vector<vector<int> > cons_table(len);
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

    vector<int> GetWeightedColors(const ClustersSet &mapping_descr) const {
        size_t len = mapping_descr.size();
        vector<int> colors(len);
        vector<int> cluster_size(len);
        vector<int> max_size(len);
        vector<size_t> prev(len);
        size_t i = 0;
        for (i = 0; i < len; i++) {
            prev[i] = -1;
        }
        for (i = 0; i < len; i++) {
            colors[i] = UNDEF_COLOR;
        }
        i = 0;
        for (auto i_iter = mapping_descr.begin(); i_iter != mapping_descr.end();
                ++i_iter, ++i) {
            cluster_size[i] = i_iter->size;
        }

        auto cons_table = FillConnectionsTable(mapping_descr);
        int cur_color = 0;
        while (true) {
            for (i = 0; i < len; i++) {
                max_size[i] = 0;
                prev[i] = -1ul;
            }
            i = 0;
            for (auto i_iter = mapping_descr.begin(); i_iter != mapping_descr.end();
                        ++i_iter, ++i) {
                if (colors[i] != UNDEF_COLOR) continue;
                max_size[i] = cluster_size[i];
                for (size_t j = 0; j < i; j ++) {
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
                real_maxi --;
            }
        }
        return colors;
    }

    //GapDescription CreateGapDescription(const KmerCluster<debruijn_graph::Graph>& a,
    //                                    const KmerCluster<debruijn_graph::Graph>& b,
    //                                    const Sequence& read) const {

    //}

    OneReadMapping AddGapDescriptions(const vector<typename ClustersSet::iterator> &start_clusters,
                                      const vector<typename ClustersSet::iterator> &end_clusters,
                                      const vector<vector<EdgeId>> &sortedEdges, const Sequence &s,
                                      const vector<bool> &block_gap_closer) const {
        DEBUG("adding gaps between subreads");
        vector<GapDescription> illumina_gaps;
        for (size_t i = 0; i + 1 < sortedEdges.size() ; i++) {
            if (block_gap_closer[i])
                continue;
            size_t j = i + 1;
            EdgeId before_gap = sortedEdges[i][sortedEdges[i].size() - 1];
            EdgeId after_gap = sortedEdges[j][0];
//do not add "gap" for rc-jumping
            if (before_gap != after_gap && before_gap != g_.conjugate(after_gap)) {
                if (TopologyGap(before_gap, after_gap, true)) {
                    if (start_clusters[j]->CanFollow(*end_clusters[i])) {
                        const auto &a = *end_clusters[i];
                        const auto &b = *start_clusters[j];
                        size_t seq_start = a.sorted_positions[a.last_trustable_index].read_position + debruijn_k;
                        size_t seq_end = b.sorted_positions[b.first_trustable_index].read_position;
                        //FIXME check index: left_offsed should be potentially equal to edge_length
                        size_t left_offset = a.sorted_positions[a.last_trustable_index].edge_position;
                        size_t right_offset = b.sorted_positions[b.first_trustable_index].edge_position;
                        auto gap = GapDescription::CreateFixOverlap(g_, s, 
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
        return OneReadMapping(sortedEdges, illumina_gaps, vector<size_t>(0));
    }

    void ProcessCluster(vector<typename ClustersSet::iterator> &cur_cluster,
                        vector<typename ClustersSet::iterator> &start_clusters,
                        vector<typename ClustersSet::iterator> &end_clusters,
                        vector<vector<EdgeId>> &sortedEdges, const Sequence &s,
                        vector<bool> &block_gap_closer) const {
        sort(cur_cluster.begin(), cur_cluster.end(),
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
                        sortedEdges.push_back(cur_sorted);
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

    //FIXME change to io::SingleRead
    OneReadMapping GetReadAlignment(Sequence &s) const {
        ClustersSet mapping_descr = GetBWAClusters(s); //GetOrderClusters(s);
        vector<int> colors = GetWeightedColors(mapping_descr);
        size_t len =  mapping_descr.size();
        vector<size_t> real_length;
        vector<vector<EdgeId>> sortedEdges;
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
                ProcessCluster(cur_cluster, start_clusters, end_clusters, sortedEdges, s, block_gap_closer);
            }
        }
        return AddGapDescriptions(start_clusters, end_clusters, sortedEdges, s, block_gap_closer);
    }

    std::pair<int, int> GetPathLimits(const KmerCluster<Graph> &a,
                                      const KmerCluster<Graph> &b,
                                      int s_add_len, int e_add_len) const {
        int start_pos = a.sorted_positions[a.last_trustable_index].read_position;
        int end_pos = b.sorted_positions[b.first_trustable_index].read_position;
        int seq_len = -start_pos + end_pos;
        //int new_seq_len =
//TODO::something more reasonable
        int path_min_len = max(int(floor((seq_len - int(debruijn_k)) * pb_config_.path_limit_pressing)), 0);
        int path_max_len = (int) ((double) (seq_len + (int) debruijn_k * 2) * pb_config_.path_limit_stretching);
        if (seq_len < 0) {
            INFO("suspicious negative seq_len " << start_pos << " " << end_pos << " " << path_min_len << " " << path_max_len);
            if (path_max_len < 0)
            return std::make_pair(-1, -1);
        }
        path_min_len = max(path_min_len - int(s_add_len + e_add_len), 0);
        path_max_len = max(path_max_len - int(s_add_len + e_add_len), 0);
        return std::make_pair(path_min_len, path_max_len);
    }

//0 - No, 1 - Yes
    int IsConsistent(const KmerCluster<Graph> &a,
                     const KmerCluster<Graph> &b) const {
        EdgeId a_edge = a.edgeId;
        EdgeId b_edge = b.edgeId;
        size_t a_id =  g_.int_id(a_edge);
        size_t b_id =  g_.int_id(b_edge);
        DEBUG("clusters on " << a_id << " and " << b_id );
//FIXME
//Is this check useful?
        if (a.sorted_positions[a.last_trustable_index].read_position + (int) pb_config_.max_path_in_dijkstra <
            b.sorted_positions[b.first_trustable_index].read_position) {
            DEBUG ("Clusters are too far in read");
            return 0;
        }
        VertexId start_v = g_.EdgeEnd(a_edge);
        size_t addition = g_.length(a_edge);
        VertexId end_v = g_.EdgeStart(b_edge);
        pair<VertexId, VertexId> vertex_pair = make_pair(start_v, end_v);

        size_t result = size_t(-1);
        bool not_found = true;
        auto distance_it = distance_cashed.begin();
#pragma omp critical(pac_index)
        {
            distance_it = distance_cashed.find(vertex_pair);
            not_found = (distance_it == distance_cashed.end());
        }
        if (not_found) {
            omnigraph::DijkstraHelper<debruijn_graph::Graph>::BoundedDijkstra dijkstra(
                    omnigraph::DijkstraHelper<debruijn_graph::Graph>::CreateBoundedDijkstra(g_, pb_config_.max_path_in_dijkstra, pb_config_.max_vertex_in_dijkstra));
            dijkstra.Run(start_v);
            if (dijkstra.DistanceCounted(end_v)) {
                result = dijkstra.GetDistance(end_v);
            }
#pragma omp critical(pac_index)
            {
                distance_it = distance_cashed.insert({vertex_pair, result}).first;
            }
        } else {
            TRACE("taking from cashed");
        }

        result = distance_it->second;
        DEBUG (result);
        if (result == size_t(-1)) {
            return 0;
        }
        //TODO: Serious optimization possible
        if (similar_in_graph(a.sorted_positions[1], b.sorted_positions[0], result + addition)) {
            DEBUG(" similar! ")
            return 1;
        } else {
            if  (a.size > 300 && b.size > 300 && - a.sorted_positions[1].edge_position +
                                                         result + addition + b.sorted_positions[0].edge_position  <=
                                                 b.sorted_positions[0].read_position - a.sorted_positions[1].read_position + 2 * g_.k()) {
                WARN("overlapping range magic worked");
                WARN("Ranges:" << a.str(g_) << " " << b.str(g_) <<" llength and dijkstra shift :" << addition << " " << result);
                return 1;

            }
            return 0;
        }
    }

    //FIXME use common function
    string PathToString(const vector<EdgeId>& path) const {
        string res = "";
        for (auto iter = path.begin(); iter != path.end(); iter++) {
            size_t len = g_.length(*iter);
            string tmp = g_.EdgeNucls(*iter).First(len).str();
            res = res + tmp;
        }
        return res;
    }

    vector<EdgeId> BestScoredPath(const Sequence &s, VertexId start_v, VertexId end_v,
                                  int path_min_length, int path_max_length,
                                  int start_pos, int end_pos, string &s_add,
                                  string &e_add) const {
        TRACE(" Traversing tangled region. Start and end vertices resp: " << g_.int_id(start_v) <<" " << g_.int_id(end_v));
        omnigraph::PathStorageCallback<Graph> callback(g_);
        ProcessPaths(g_,
                    path_min_length, path_max_length,
                    start_v, end_v,
                    callback);
        vector<vector<EdgeId> > paths = callback.paths();
        TRACE("taking subseq" << start_pos <<" "<< end_pos <<" " << s.size());
        int s_len = int(s.size());
        string seq_string = s.Subseq(start_pos, min(end_pos + 1, s_len)).str();
        size_t best_path_ind = paths.size();
        int best_score = STRING_DIST_INF;
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
        for (size_t i = 0; i < paths.size(); i++) {
            string cur_string = s_add + PathToString(paths[i]) + e_add;
            TRACE("cur_string: " << cur_string <<"\n seq_string " << seq_string);
            if (paths.size() > 1 && paths.size() < 10) {
                TRACE("candidate path number "<< i << " , len " << cur_string.length());
                TRACE("graph candidate: " << cur_string);
                TRACE("in pacbio read: " << seq_string);
                for (auto j_iter = paths[i].begin(); j_iter != paths[i].end();
                        ++j_iter) {
                    DEBUG(g_.int_id(*j_iter));
                }
            }
            int cur_score = StringDistance(cur_string, seq_string);
            if (paths.size() > 1 && paths.size() < 10) {
                DEBUG("score: "<< cur_score);
            }
            if (cur_score < best_score) {
                best_score = cur_score;
                best_path_ind = i;
            }
        }
        TRACE(best_score);
        if (best_score == STRING_DIST_INF)
            return vector<EdgeId>(0);
        if (paths.size() > 1 && paths.size() < 10) {
            TRACE("best score found! Path " <<best_path_ind <<" score "<< best_score);
        }
        return paths[best_path_ind];
    }

    // Short read alignment
    //omnigraph::MappingPath<EdgeId> GetShortReadAlignment(const Sequence &s) const {
    //    ClustersSet mapping_descr = GetOrderClusters(s);
    //    map<EdgeId, KmerCluster<Graph> > largest_clusters;

    //    //Selecting the biggest cluster for each edge
    //    for (auto iter = mapping_descr.begin(); iter != mapping_descr.end(); ++iter) {

    //        auto first_cluster = iter->sorted_positions[iter->first_trustable_index];
    //        auto last_cluster = iter->sorted_positions[iter->last_trustable_index];
    //        int read_range = last_cluster.read_position - first_cluster.read_position;
    //        int edge_range = last_cluster.edge_position - first_cluster.edge_position;
    //        int cluster_szie = iter->last_trustable_index - iter->first_trustable_index;
    //        if (cluster_szie > 2 * read_range || edge_range < 0 || 2 * edge_range < read_range || edge_range > 2 * read_range) {
    //            //skipping cluster
    //            continue;
    //        }

    //        auto edge_cluster = largest_clusters.find(iter->edgeId);
    //        if (edge_cluster != largest_clusters.end()) {
    //            if (edge_cluster->second.last_trustable_index - edge_cluster->second.first_trustable_index
    //                    < iter->last_trustable_index - iter->first_trustable_index) {

    //                edge_cluster->second = *iter;
    //            }
    //        } else {
    //            largest_clusters.insert(make_pair(iter->edgeId, *iter));
    //        }
    //    }

    //    omnigraph::MappingPath<EdgeId> result;
    //    for (auto iter = largest_clusters.begin(); iter != largest_clusters.end(); ++iter) {
    //        auto first_cluster = iter->second.sorted_positions[iter->second.first_trustable_index];
    //        auto last_cluster = iter->second.sorted_positions[iter->second.last_trustable_index];
    //        omnigraph::MappingRange range(Range(first_cluster.read_position, last_cluster.read_position),
    //                                      Range(first_cluster.edge_position, last_cluster.edge_position));
    //        result.join({iter->second.edgeId, range});
    //    }

    //    return result;
    //}

    //std::pair<EdgeId, size_t> GetUniqueKmerPos(const RtSeq& kmer) const {
    //    KeyWithHash kwh = tmp_index.ConstructKWH(kmer);

    //    if (tmp_index.valid(kwh.key())) {
    //        auto keys = tmp_index.get(kwh);
    //        if (keys.size() == 1) {
    //            return make_pair(keys[0].edge_id, keys[0].offset);
    //        }
    //    }
    //    return std::make_pair(EdgeId(0), -1u);
    //}


};

//template<class Graph>
//typename PacBioMappingIndex<Graph>::MappingDescription PacBioMappingIndex<Graph>::GetSeedsFromRead(const Sequence &s) const {
//    MappingDescription res;
//    //WARNING: removed read_count from here to make const methods
//    int local_read_count = 0;
//    ++local_read_count;
//    if (s.size() < pacbio_k)
//        return res;
//
//    //RtSeq kmer = s.start<RtSeq>(pacbio_k);
//    KeyWithHash kwh = tmp_index.ConstructKWH(s.start<RtSeq>(pacbio_k));
//
//    for (size_t j = pacbio_k; j < s.size(); ++j) {
//        kwh = kwh << s[j];
//        if (!tmp_index.valid(kwh.key())) {
////          INFO("not valid kmer");
//            continue;
//        }
//        auto keys = tmp_index.get(kwh);
//        TRACE("Valid key, size: "<< keys.size());
//
//        int quality = (int) keys.size();
//        if (quality > 1000) {
//            DEBUG ("Ignoring repretive kmer")
//            continue;
//        }
//        for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
//
//            TRACE("and quality:" << quality);
//            int offset = (int)iter->offset;
//            int s_stretched = int ((double)s.size() * 1.2 + 50);
//            int edge_len = int(g_.length(iter->edge_id));
//            //No alignment in vertex, and further than s+eps bp from edge ends;
//            bool correct_alignment = offset > int(debruijn_k - pacbio_k) && offset < edge_len;
//            if (ignore_map_to_middle) {
//                correct_alignment &= (offset < int(debruijn_k - pacbio_k) + s_stretched || offset > edge_len - s_stretched);
//            }
//            if (correct_alignment) {
//                res[iter->edge_id].push_back(MappingInstance((int) iter->offset, (int) (j - pacbio_k + 1), quality));
//            }
//        }
//    }
//
//    for (auto iter = res.begin(); iter != res.end(); ++iter) {
//        sort(iter->second.begin(), iter->second.end());
//        DEBUG("read count "<< local_read_count);
//        DEBUG("edge: " << g_.int_id(iter->first) << "size: " << iter->second.size());
//        for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); j_iter++) {
//            DEBUG(j_iter->str());
//        }
//    }
//
//    return res;
//}

}
