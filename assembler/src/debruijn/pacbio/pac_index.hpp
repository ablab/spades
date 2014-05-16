/*
 * pac_index.hpp
 *
 *  Created on: Jan 21, 2013
 *      Author: lab42
 */
#pragma once

#include "indices/edge_multi_index.hpp"
#include "indices/edge_index_builders.hpp"
#include <algorithm>
#include "pacbio_read_structures.hpp"

namespace pacbio {

template<class Graph>
struct MappingDescription {

};

template<class Graph>
class PacBioMappingIndex {
public:
    typedef map<typename Graph::EdgeId, vector<MappingInstance> > MappingDescription;
    typedef pair<typename Graph::EdgeId, vector<MappingInstance> > ClusterDescription;
    typedef set<KmerCluster<Graph> > ClustersSet;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef     debruijn_graph::DeBruijnEdgeMultiIndex<typename Graph::EdgeId> Index;
    typedef typename Index::KeyWithHash KeyWithHash;

private:
    DECL_LOGGER("PacIndex")

    const Graph &g_;
    size_t pacbio_k;
    size_t debruijn_k;
    const static int short_edge_cutoff = 0;
    const static size_t min_cluster_size = 8;
    const static int max_similarity_distance = 500;
    double compression_cutoff;
    double domination_cutoff;
    set<Sequence> banned_kmers;
    debruijn_graph::DeBruijnEdgeMultiIndex<typename Graph::EdgeId> tmp_index;
    map<pair<VertexId, VertexId>, vector<size_t> > distance_cashed;
    size_t read_count;
    bool ignore_map_to_middle;
public:
    MappingDescription Locate(const Sequence &s) const;

    PacBioMappingIndex(const Graph &g, size_t k, size_t debruijn_k_, bool ignore_map_to_middle)
            : g_(g),
              pacbio_k(k),
              debruijn_k(debruijn_k_),
              tmp_index((unsigned) pacbio_k, cfg::get().output_dir), ignore_map_to_middle(ignore_map_to_middle) {
        DEBUG("PB Mapping Index construction started");

        typedef typename debruijn_graph::EdgeIndexHelper<debruijn_graph::DeBruijnEdgeMultiIndex<typename Graph::EdgeId>>::GraphPositionFillingIndexBuilderT Builder;

        Builder().BuildIndexFromGraph(tmp_index, g_);
        FillBannedKmers();
        compression_cutoff = cfg::get().pb.compression_cutoff;  // 0.6
        domination_cutoff = cfg::get().pb.domination_cutoff;  //1.5
        //INFO(tmp_index.size());
        read_count = 0;
    }

    void FillBannedKmers() {
        for (int i = 0; i < 4; i++) {
            auto base = nucl((unsigned char) i);
            for (int j = 0; j < 4; j++) {
                auto other = nucl((unsigned char) j);
                for (size_t other_pos = 0; other_pos < pacbio_k; other_pos++) {
                    string s = "";
                    for (size_t k = 0; k < pacbio_k; k++) {
                        if (k != other_pos)
                            s += base;
                        else
                            s += other;
                    }
                    banned_kmers.insert(Sequence(s));
                }
            }
        }
    }

    bool similar(const MappingInstance &a, const MappingInstance &b,
                        int shift = 0) const {
        if (b.read_position + shift < a.read_position) {
            return similar(b, a, -shift);
        } else if (b.read_position == a.read_position) {
            return (abs(int(b.edge_position) + shift - int(a.edge_position)) < 2);
        } else {
            return ((b.edge_position + shift - a.edge_position
                    >= (b.read_position - a.read_position) * compression_cutoff)
                    && (b.edge_position + shift - a.edge_position
                            <= (b.read_position - a.read_position)
                                    / compression_cutoff));
        }
    }

    void dfs_cluster(vector<int> &used, vector<MappingInstance> &to_add,
                     const int cur_ind,
                     const typename MappingDescription::iterator iter) const {
        size_t len = iter->second.size();
        for (size_t k = 0; k < len; k++) {
            if (!used[k] && similar(iter->second[cur_ind], iter->second[k])) {
                to_add.push_back(iter->second[k]);
                used[k] = 1;
                dfs_cluster(used, to_add, (int) k, iter);
            }
        }
    }

    void dfs_cluster_norec(vector<int> &used, vector<MappingInstance> &to_add,
                     const size_t cur_ind,
                     const typename MappingDescription::iterator iter, vector<vector<size_t> > &similarity_list) const {
        std::deque<size_t> stack;
        stack.push_back(cur_ind);
        used[cur_ind] = 1;
        while (stack.size() > 0) {
            size_t k = stack.back();
            stack.pop_back();
            to_add.push_back(iter->second[k]);

            for (size_t i = 0; i < similarity_list[k].size(); i++) {
                if (!used[similarity_list[k][i]]) {
                    stack.push_back(similarity_list[k][i]);
                    used[similarity_list[k][i]] = 1;
                }
            }
        }
    }

    ClustersSet GetClusters(const Sequence &s) const {
        MappingDescription descr = Locate(s);
        ClustersSet res;
        TRACE(read_count << " read_count");

        DEBUG(descr.size() <<"  clusters");
        for (auto iter = descr.begin(); iter != descr.end(); ++iter) {
            size_t edge_id = g_.int_id(iter->first);
            DEBUG(edge_id);
            sort(iter->second.begin(), iter->second.end(), ReadPositionComparator());
            set<vector<MappingInstance> > edge_cluster_set;
            size_t len = iter->second.size();
            vector<vector<size_t> > similarity_list(len);
            int cnt = 0;
            for (size_t i = 0; i < len; i++){
                for (size_t j = i + 1; j < len; j++){
                    if (iter->second[i].read_position + max_similarity_distance < iter->second[j].read_position) {
                        break;
                    }
                    if (similar(iter->second[i], iter->second[j])) {
                        similarity_list[i].push_back(j);
                        similarity_list[j].push_back(i);
                        cnt ++;
                        if (cnt % 10000 == 0) {
                            DEBUG(cnt);
                        }
                    }
                }
            }
            DEBUG(len <<"  kmers in cluster");
            vector<int> used(len);
            for (size_t i = 0; i < len; i++) {
                if (!used[i]) {
                    used[i] = 1;
                    vector<MappingInstance> to_add;
                    //to_add.push_back(iter->second[i]);
                    //dfs_cluster(used, to_add, (int) i, iter);
                    dfs_cluster_norec(used, to_add, (int) i, iter, similarity_list);
                    sort(to_add.begin(), to_add.end(), ReadPositionComparator());
                    DEBUG(to_add.size()<<" subcluster size");
                    size_t count = 1;
                    size_t longest_len = 0;
                    auto cur_start = to_add.begin();
                    auto best_start = to_add.begin();
                    int cur_edge_len = (int)g_.length(iter->first);
                    int j = 0;
                        DEBUG("new_cluster sz " << to_add.size());
                    for (auto j_iter = to_add.begin();
                            j_iter < to_add.end() - 1; j_iter++, j++) {
//Do not spilt clusters in the middle, only beginning is interesting.
                        if ((j * 5 < (int)to_add.size() || (j + 1) * 5 > (int)to_add.size() * 4) &&
                            !similar(*j_iter, *(j_iter + 1))
                            && (j < 5000 && cur_edge_len - j < 5000)) {
                            if (longest_len < count) {
                                longest_len = count;
                                best_start = cur_start;
                            }
                            count = 1;
                            cur_start = j_iter + 1;
                        } else {
                            count++;
                        }
                    }
                    if (longest_len < count) {
                        longest_len = count;
                        best_start = cur_start;
                    }
                    vector<MappingInstance> filtered(best_start,
                                                     best_start + longest_len);
                    if ((count < to_add.size()
                            && to_add.size() > min_cluster_size)) {
                        DEBUG("in cluster size " << to_add.size() << ", " << to_add.size() - longest_len<<"were removed as trash")
                    }
                    TRACE("adding cluster "" edge "<< edge_id << " len " <<to_add.size() )
                    res.insert(KmerCluster<Graph>(iter->first, filtered));
                }
            }
        }
        FilterClusters(res);
        return res;
    }

    //filter clusters that are too small or fully located on a vertex or dominated by some other cluster.
    void FilterClusters(ClustersSet &clusters) const {
        for (auto i_iter = clusters.begin(); i_iter != clusters.end();) {
            size_t edge_id = g_.int_id(i_iter->edgeId);

            int len = (int) g_.length(i_iter->edgeId);
            auto sorted_by_edge = i_iter->sorted_positions;
            sort(sorted_by_edge.begin(), sorted_by_edge.end());
            double good = 0;
            DEBUG("filtering cluster of size " << sorted_by_edge.size());
            DEBUG(edge_id <<" : edgeId");
            for (auto iter = sorted_by_edge.begin();
                    iter < sorted_by_edge.end(); iter++) {
                if (iter->IsUnique())
                    good++;
//              good += 1.0 / (iter->quality * iter->quality);
            }
            DEBUG("good " << good);

            if (good < min_cluster_size || (len < short_edge_cutoff)) {
                if (len < short_edge_cutoff) {
                    DEBUG("Life is too long, and edge is too short!");
                }
                auto tmp_iter = i_iter;
                tmp_iter++;
                clusters.erase(i_iter);
                i_iter = tmp_iter;
            } else {
                if (sorted_by_edge[0].edge_position >= len
                        || sorted_by_edge[i_iter->size - 1].edge_position
                                <= int(debruijn_k) - int(pacbio_k)) {
                    DEBUG("All anchors in vertex");
                    auto tmp_iter = i_iter;
                    tmp_iter++;
                    clusters.erase(i_iter);
                    i_iter = tmp_iter;
                } else {
                    i_iter++;
                }
            }
        }
        for (auto i_iter = clusters.begin(); i_iter != clusters.end();) {
            size_t edge_id = g_.int_id(i_iter->edgeId);
            auto sorted_by_edge = i_iter->sorted_positions;

            DEBUG("filtering  with cluster edge, stage 2 "<< edge_id << " len " << sorted_by_edge.size() << " clusters still alive: "<< clusters.size());
            for (auto j_iter = clusters.begin(); j_iter != clusters.end();) {
                if (i_iter != j_iter) {
                    if (dominates(*i_iter, *j_iter)) {
                        TRACE("cluster is dominated");
                        auto tmp_iter = j_iter;
                        tmp_iter++;
                        TRACE("cluster on edge " << g_.int_id(j_iter->edgeId));
                        TRACE("erased - dominated");
                        clusters.erase(j_iter);
                        j_iter = tmp_iter;
                    } else {
                        j_iter++;
                    }
                } else {
                    j_iter++;
                }
            }
            DEBUG("cluster size "<< i_iter->sorted_positions.size() << "survived filtering");
            i_iter++;
        }
    }

    // is "non strictly dominates" required?
    inline bool dominates(const KmerCluster<Graph> &a,
                          const KmerCluster<Graph> &b) const {
        size_t a_size = a.size;
        size_t b_size = b.size;
        if ((double) a_size < (double) b_size * domination_cutoff
                || a.sorted_positions[a.first_trustable_index].read_position
                        > b.sorted_positions[b.first_trustable_index].read_position
                || a.sorted_positions[a.last_trustable_index].read_position
                        < b.sorted_positions[b.last_trustable_index].read_position) {
            return false;
        } else {
            return true;
        }
    }

    vector<EdgeId> FillGapsInCluster(vector<pair<size_t, typename ClustersSet::iterator> > &cur_cluster,
                                     const Sequence &s) {
        vector<EdgeId> cur_sorted;
        EdgeId prev_edge = EdgeId(0);
        if (read_count == 151) {
            DEBUG("edges in cluster:")
            for (auto iter = cur_cluster.begin(); iter != cur_cluster.end();
                    ++iter) {
                EdgeId cur_edge = iter->second->edgeId;
                DEBUG(g_.int_id(cur_edge) << " " << iter->second->sorted_positions.size());
            }
        }
        for (auto iter = cur_cluster.begin(); iter != cur_cluster.end();
                ++iter) {
            EdgeId cur_edge = iter->second->edgeId;
            if (prev_edge != EdgeId(0)) {
//Need to find sequence of edges between clusters
                VertexId start_v = g_.EdgeEnd(prev_edge);
                VertexId end_v = g_.EdgeStart(cur_edge);
                auto prev_iter = iter - 1;
                MappingInstance cur_first_index =
                        iter->second->sorted_positions[iter->second
                                ->first_trustable_index];
                MappingInstance prev_last_index = prev_iter->second
                        ->sorted_positions[prev_iter->second
                        ->last_trustable_index];

                if (start_v != end_v ||
                    (start_v == end_v &&
                     (double) (cur_first_index.read_position - prev_last_index.read_position) >
                     (double) (cur_first_index.edge_position + (int) g_.length(prev_edge) - prev_last_index.edge_position) * 1.3)) {
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
                    pair<int, int> limits = GetPathLimits(*(prev_iter->second),
                                                          *(iter->second),
                                                          (int) s_add.length(),
                                                          (int) e_add.length());
                    if (limits.first == -1)
                        return vector<EdgeId>(0);

                    vector<EdgeId> intermediate_path = BestScoredPath(s, start_v, end_v, limits.first, limits.second, seq_start, seq_end, s_add, e_add);
                    if (intermediate_path.size() == 0) {
                        DEBUG("Tangled region between edgees "<< g_.int_id(prev_edge) << " " << g_.int_id(cur_edge) << " is not closed, additions from edges: " << int(g_.length(prev_edge)) - int(prev_last_index.edge_position) <<" " << int(cur_first_index.edge_position) - int(debruijn_k - pacbio_k ) << " and seq "<< - seq_start + seq_end);
                        if (cfg::get().pb.additional_debug_info) {
                            DEBUG(" escpected gap length: " << -int(g_.length(prev_edge)) + int(prev_last_index.edge_position) - int(cur_first_index.edge_position) + int(debruijn_k - pacbio_k ) - seq_start + seq_end);
                            PathStorageCallback<Graph> callback(g_);
                            PathProcessor<Graph> path_processor(g_, 0, 4000,
                                                                start_v, end_v,
                                                                callback);
                            path_processor.Process();
                            vector<vector<EdgeId> > paths = callback.paths();
                            stringstream s_buf;
                            for (auto p_iter = paths.begin();
                                    p_iter != paths.end(); p_iter++) {
                                size_t tlen = 0;
                                for (auto path_iter = p_iter->begin();
                                        path_iter != p_iter->end();
                                        path_iter++) {
                                    tlen += g_.length(*path_iter);
                                }
                                s_buf << tlen << " ";
                            }
                            DEBUG(s_buf.str());
                        }
                        return intermediate_path;
                    }
                    for (auto j_iter = intermediate_path.begin(); j_iter != intermediate_path.end(); j_iter++) {
                        cur_sorted.push_back(*j_iter);
                    }
                }
            }
            cur_sorted.push_back(cur_edge);
            prev_edge = cur_edge;
        }
        return cur_sorted;
    }

    bool TopologyGap(EdgeId first, EdgeId second, bool oriented) const {
        bool res = (g_.IsDeadStart(g_.EdgeStart(first)) && g_.IsDeadEnd(g_.EdgeEnd(second)));
        if (!oriented)
            res |= g_.IsDeadEnd(g_.EdgeEnd(first)) && g_.IsDeadStart(g_.EdgeStart(second));
        return res;
    }

    vector<int> GetColors(ClustersSet &mapping_descr, Sequence &s) {
        int len = (int) mapping_descr.size();
	DEBUG("getting colors, table size "<< len);
        vector<vector<int> > table(len);
        for (int i = 0; i < len; i++) {
            table[i].resize(len);
            table[i][i] = 1;
        }
        int i = 0;
        vector<int> colors(len);
        vector<int> cluster_size(len);
        for (int i = 0; i < len; i++) {
            colors[i] = i;
        }
        for (auto i_iter = mapping_descr.begin(); i_iter != mapping_descr.end();
                ++i_iter, ++i) {
            cluster_size[i] = i_iter->size;
        }
        i = 0;
        if (len > 1) {
            TRACE(len << "clusters");
        }
        for (auto i_iter = mapping_descr.begin(); i_iter != mapping_descr.end();
                ++i_iter, ++i) {
            int j = 0;
            for (auto j_iter = mapping_descr.begin();
                    j_iter != mapping_descr.end(); ++j_iter, ++j) {
                if (i_iter == j_iter)
                    continue;
                table[i][j] = IsConsistent(s, *i_iter, *j_iter);
                if (table[i][j]) {
                    int j_color = colors[j];
                    for (int k = 0; k < int(len); k++) {
                        if (colors[k] == j_color) {
                            colors[k] = colors[i];
                            TRACE("recoloring");
                            cluster_size[i] += cluster_size[k];
                            cluster_size[k] = 0;
                        }
                    }
                }
            }
        }
        return colors;
    }

    OneReadMapping<Graph> GetReadAlignment(Sequence &s) {
        ClustersSet mapping_descr = GetClusters(s);
        DEBUG("clusters got");
        int len = (int) mapping_descr.size();
        vector<size_t> real_length;

        vector<int> colors = GetColors(mapping_descr, s);
        vector<vector<EdgeId> > sortedEdges;
        vector<typename ClustersSet::iterator> start_clusters, end_clusters;
        vector<GapDescription<Graph> > illumina_gaps;
        vector<int> used(len);
	auto iter = mapping_descr.begin();
        for (int i = 0; i < len; i++, iter ++) {
            used[i] = 0; 
	    DEBUG(colors[i] <<" " << iter->str(g_));
	}
        for (int i = 0; i < len; i++) {
            if (!used[i]) {
                DEBUG("starting new subread");
                vector<pair<size_t, typename ClustersSet::iterator> > cur_cluster;
                used[i] = 1;
                int j = 0;
                int cur_color = colors[i];
                for (auto i_iter = mapping_descr.begin();
                        i_iter != mapping_descr.end(); ++i_iter, ++j) {
                    if (colors[j] == cur_color) {
                        cur_cluster.push_back(
                                make_pair(
                                        i_iter->sorted_positions.begin()
                                                ->read_position,
                                        i_iter));
                        used[j] = 1;
                    }
                }
                sort(cur_cluster.begin(), cur_cluster.end(),
                     pair_iterator_less<typename ClustersSet::iterator>());
                VERIFY(cur_cluster.size() > 0);
                auto cur_cluster_start = cur_cluster.begin();
                for (auto iter = cur_cluster.begin(); iter != cur_cluster.end();
                        ++iter) {
                    auto next_iter = iter + 1;
                    if (next_iter == cur_cluster.end()
                            || !IsConsistent(s, *(iter->second),
                                             *(next_iter->second))) {
                        if (next_iter != cur_cluster.end()) {
                            DEBUG("cluster splitted...");
                        }
                        vector<pair<size_t, typename ClustersSet::iterator> > splitted_cluster(
                                cur_cluster_start, next_iter);
                        vector<EdgeId> cur_sorted = FillGapsInCluster(
                                splitted_cluster, s);
                        if (cur_sorted.size() > 0) {
                            start_clusters.push_back(cur_cluster_start->second);
                            end_clusters.push_back(iter->second);
                            sortedEdges.push_back(cur_sorted);
                        }
                        cur_cluster_start = next_iter;
                    }
                }
            }
        }
	DEBUG("adding gaps between subreads");
        int alignments = int(sortedEdges.size());
        for (int i = 0; i < alignments; i++) {
            for (int j = 0; j < alignments; j++) {
                EdgeId before_gap = sortedEdges[j][sortedEdges[j].size() - 1];
                EdgeId after_gap = sortedEdges[i][0];
//do not add "gap" for rc-jumping
                if (before_gap != after_gap
                        && before_gap != g_.conjugate(after_gap)) {
                    if (i != j && TopologyGap(before_gap, after_gap, true)) {
                        if (start_clusters[j]->CanFollow(*end_clusters[i])) {
                            illumina_gaps.push_back(
                                    GapDescription<Graph>(*end_clusters[i],
                                                          *start_clusters[j], s,
                                                          (int) pacbio_k));
                        }

                    }
                }
            }
        }
        return OneReadMapping<Graph>(sortedEdges, illumina_gaps, real_length);
    }

    std::pair<int, int> GetPathLimits(const KmerCluster<Graph> &a,
                                      const KmerCluster<Graph> &b,
                                      int s_add_len, int e_add_len) {
        int start_pos = a.sorted_positions[a.last_trustable_index].read_position;
        int end_pos = b.sorted_positions[b.first_trustable_index].read_position;
        int seq_len = -start_pos + end_pos;
        PathStorageCallback<Graph> callback(g_);
//TODO::something more reasonable
        int path_min_len = max(int(floor((seq_len - int(debruijn_k)) * cfg::get().pb.path_limit_pressing)), 0);
        int path_max_len = (int) ((double) (seq_len + (int) debruijn_k) * cfg::get().pb.path_limit_stretching);
        if (seq_len < 0) {
            DEBUG("suspicious negative seq_len " << start_pos << " " << end_pos << " " << path_min_len << " " << path_max_len);
            return std::make_pair(-1, -1);
        }
        path_min_len = max(path_min_len - int(s_add_len + e_add_len), 0);
        path_max_len = max(path_max_len - int(s_add_len + e_add_len), 0);
        return std::make_pair(path_min_len, path_max_len);
    }

//0 - No, 1 - Yes
    int IsConsistent(Sequence &s, const KmerCluster<Graph> &a,
                     const KmerCluster<Graph> &b) {
        EdgeId a_edge = a.edgeId;
        EdgeId b_edge = b.edgeId;
        if (abs(a.sorted_positions[a.last_trustable_index].read_position - b.sorted_positions[b.first_trustable_index].read_position) > 5000) {
            return 0;
        }
        VertexId start_v = g_.EdgeEnd(a_edge);
        size_t addition = g_.length(a_edge);
        VertexId end_v = g_.EdgeStart(b_edge);
        pair<VertexId, VertexId> vertex_pair = make_pair(start_v, end_v);
        vector<size_t> result;
        if (distance_cashed.find(vertex_pair) == distance_cashed.end()) {
            DistancesLengthsCallback<Graph> callback(g_);
            PathProcessor<Graph> path_processor(g_, 0, s.size() / 3, start_v,
                                                end_v, callback);
            path_processor.Process();
            result = callback.distances();
            distance_cashed[vertex_pair] = result;
        }
        result = distance_cashed[vertex_pair];
        //TODO: Serious optimization possible
        for (size_t i = 0; i < result.size(); i++) {
            for (auto a_iter = a.sorted_positions.begin();
                    a_iter != a.sorted_positions.end(); ++a_iter) {
                if (a_iter - a.sorted_positions.begin() > 500 &&  a.sorted_positions.end() - a_iter >500) continue;
                int cnt = 0;
                for (auto b_iter = b.sorted_positions.begin();
                        b_iter != b.sorted_positions.end() && cnt <500; ++b_iter, cnt ++) {
                    if (similar(*a_iter, *b_iter,
                                (int) (result[i] + addition))) {
                        return 1;
                    }
                }
                cnt = 0;
                if (b.sorted_positions.size() > 500) {
                    for (auto b_iter = b.sorted_positions.end() - 1;
                                            b_iter != b.sorted_positions.begin() && cnt < 500; --b_iter, cnt ++) {
                        if (similar(*a_iter, *b_iter,
                                    (int) (result[i] + addition))) {
                            return 1;
                        }
                    }
                }
            }
        }
        return 0;

    }

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
                                  string &e_add) {
        DEBUG(" Traversing tangled region. Start and end vertices resp: " << g_.int_id(start_v) <<" " << g_.int_id(end_v));
        PathStorageCallback<Graph> callback(g_);
        PathProcessor<Graph> path_processor(g_,
                                            path_min_length, path_max_length,
                                            start_v, end_v,
                                            callback);
        path_processor.Process();
        vector<vector<EdgeId> > paths = callback.paths();
        DEBUG("taking subseq" << start_pos <<" "<< end_pos <<" " << s.size());
        int s_len = int(s.size());
        string seq_string = s.Subseq(start_pos, min(end_pos + 1, s_len)).str();
        size_t best_path_ind = paths.size();
        size_t best_score = 1000000000;
        DEBUG("need to find best scored path between "<<paths.size()<<" , seq_len " << seq_string.length());
        if (paths.size() == 0)
            return vector<EdgeId>(0);
        for (size_t i = 0; i < paths.size(); i++) {
            string cur_string = s_add + PathToString(paths[i]) + e_add;
            if (paths.size() > 1 && paths.size() < 10) {
                TRACE("candidate path number "<< i << " , len " << cur_string.length());
                TRACE("graph candidate: " << cur_string);
                TRACE("in pacbio read: " << seq_string);
                for (auto j_iter = paths[i].begin(); j_iter != paths[i].end();
                        ++j_iter) {
                    DEBUG(g_.int_id(*j_iter));
                }
            }
            size_t cur_score = StringDistance(cur_string, seq_string);
            if (paths.size() > 1 && paths.size() < 10) {
                DEBUG("score: "<< cur_score);
            }
            if (cur_score < best_score) {
                best_score = cur_score;
                best_path_ind = i;
            }
        }
        if (best_score == 1000000000)
            return vector<EdgeId>(0);
        if (paths.size() > 1 && paths.size() < 10) {
            DEBUG("best score found! Path " <<best_path_ind <<" score "<< best_score);
        }
        return paths[best_path_ind];
    }

    // Short read alignment
    MappingPath<EdgeId> GetShortReadAlignment(const Sequence &s) const {
        ClustersSet mapping_descr = GetClusters(s);
        map<EdgeId, KmerCluster<Graph> > largest_clusters;

        //Selecting the biggest cluster for each edge
        for (auto iter = mapping_descr.begin(); iter != mapping_descr.end(); ++iter) {

            auto first_cluster = iter->sorted_positions[iter->first_trustable_index];
            auto last_cluster = iter->sorted_positions[iter->last_trustable_index];
            int read_range = last_cluster.read_position - first_cluster.read_position;
            int edge_range = last_cluster.edge_position - first_cluster.edge_position;
            int cluster_szie = iter->last_trustable_index - iter->first_trustable_index;
            if (cluster_szie > 2 * read_range || edge_range < 0 || 2 * edge_range < read_range || edge_range > 2 * read_range) {
                //skipping cluster
                continue;
            }

            auto edge_cluster = largest_clusters.find(iter->edgeId);
            if (edge_cluster != largest_clusters.end()) {
                if (edge_cluster->second.last_trustable_index - edge_cluster->second.first_trustable_index
                        < iter->last_trustable_index - iter->first_trustable_index) {

                    edge_cluster->second = *iter;
                }
            }
            else {
                largest_clusters.insert(make_pair(iter->edgeId, *iter));
            }
        }

        MappingPath<EdgeId> result;
        for (auto iter = largest_clusters.begin(); iter != largest_clusters.end(); ++iter) {
            auto first_cluster = iter->second.sorted_positions[iter->second.first_trustable_index];
            auto last_cluster = iter->second.sorted_positions[iter->second.last_trustable_index];
            MappingRange range(Range(first_cluster.read_position, last_cluster.read_position),
                    Range(first_cluster.edge_position, last_cluster.edge_position));
            result.join(MappingPath<EdgeId>(vector<EdgeId>(1, iter->second.edgeId), vector<MappingRange>(1, range)));
        }

        return result;
    }

    pair<EdgeId, size_t> GetUniqueKmerPos(const runtime_k::RtSeq& kmer) const {
        KeyWithHash kwh = tmp_index.ConstructKWH(kmer);

        if (tmp_index.valid(kwh.key())) {
            auto keys = tmp_index.get(kwh);
            if (keys.size() == 1) {
                return make_pair(keys[0].edge_id, keys[0].offset);
            }
        }
        return make_pair(EdgeId(0), -1u);
    }


};

template<class Graph>
typename PacBioMappingIndex<Graph>::MappingDescription PacBioMappingIndex<Graph>::Locate(const Sequence &s) const {
    MappingDescription res;
    //WARNING: removed read_count from here to make const methods
    int local_read_count = 0;
    ++local_read_count;
    if (s.size() < pacbio_k)
        return res;

    //runtime_k::RtSeq kmer = s.start<runtime_k::RtSeq>(pacbio_k);
    KeyWithHash kwh = tmp_index.ConstructKWH(s.start<runtime_k::RtSeq>(pacbio_k));

    for (size_t j = pacbio_k; j < s.size(); ++j) {
        kwh = kwh << s[j];
        if (!tmp_index.valid(kwh.key())) {
//          INFO("not valid kmer");
            continue;
        }
        auto keys = tmp_index.get(kwh);
        TRACE("Valid key, size: "<< keys.size());

        for (auto iter = keys.begin(); iter != keys.end(); ++iter) {

            int quality = (int) keys.size();
            TRACE("and quality:" << quality);
            if (banned_kmers.find(Sequence(kwh.key())) != banned_kmers.end())
                continue;
            int offset = (int)iter->offset;
            int s_stretched = int ((double)s.size() * 1.2 + 50);
            int edge_len = int(g_.length(iter->edge_id));
            //No alignment in vertex, and further than s+eps bp from edge ends;
            bool correct_alignment = offset > int(debruijn_k - pacbio_k) && offset < edge_len;
            if (ignore_map_to_middle) {
                correct_alignment &= (offset < int(debruijn_k - pacbio_k) + s_stretched || offset > edge_len - s_stretched);
            }
            if (correct_alignment) {
                res[iter->edge_id].push_back(MappingInstance((int) iter->offset, (int) (j - pacbio_k + 1), quality));
            }
        }
    }

    for (auto iter = res.begin(); iter != res.end(); ++iter) {
        sort(iter->second.begin(), iter->second.end());
        DEBUG("read count "<< local_read_count);
        DEBUG("edge: " << g_.int_id(iter->first) << "size: " << iter->second.size());
        for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); j_iter++) {
            DEBUG(j_iter->str());
        }
    }

    return res;
}

}
