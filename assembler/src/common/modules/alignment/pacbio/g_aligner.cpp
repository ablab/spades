//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "modules/alignment/pacbio/g_aligner.hpp"
#include "modules/alignment/pacbio/gap_filler.hpp"
#include "modules/alignment/pacbio/pacbio_read_structures.hpp"

namespace sensitive_aligner {
using namespace std;

using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;
using debruijn_graph::Graph;
void GAligner::FillGapsInCluster(const vector<QualityRange> &cur_cluster,
                                 const Sequence &s,
                                 vector<vector<debruijn_graph::EdgeId> > &edges,
                                 vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &bwa_hits) const {
    omnigraph::MappingPath<debruijn_graph::EdgeId> cur_sorted_hits;
    vector<debruijn_graph::EdgeId> cur_sorted_edges;
    EdgeId prev_edge = EdgeId();

    for (auto iter = cur_cluster.begin(); iter != cur_cluster.end();) {
        EdgeId cur_edge = iter->edgeId;
        if (prev_edge != EdgeId()) {
//Need to find sequence of edges between clusters
            VertexId start_v = g_.EdgeEnd(prev_edge);
            VertexId end_v = g_.EdgeStart(cur_edge);
            auto prev_iter = iter - 1;
            MappingInstance cur_first_index = iter->sorted_positions[iter->first_trustable_index];
            MappingInstance prev_last_index = prev_iter->sorted_positions[prev_iter->last_trustable_index];
            double read_gap_len = (double) (cur_first_index.read_position - prev_last_index.read_position);
//FIXME:: is g_.k() relevant
            double stretched_graph_len = (double) (cur_first_index.edge_position + g_.k()) +
                                         ((int) g_.length(prev_edge) - prev_last_index.edge_position) *
                                         pb_config_.path_limit_stretching;
            if (start_v != end_v || (start_v == end_v && read_gap_len > stretched_graph_len)) {
                if (start_v == end_v) {
                    DEBUG("looking for path from vertex to itself, read pos"
                          << cur_first_index.read_position << " " << prev_last_index.read_position
                          << " edge pos: " << cur_first_index.edge_position << " "
                          << prev_last_index.edge_position
                          << " edge len " << g_.length(prev_edge));
                    DEBUG("Read gap len " << read_gap_len << " stretched graph path len" << stretched_graph_len);
                }
                DEBUG(" traversing tangled hregion between " << g_.int_id(prev_edge) << " " << g_.int_id(cur_edge));
                DEBUG(" first pair" << cur_first_index.str() << " edge_len" << g_.length(cur_edge));
                DEBUG(" last pair" << prev_last_index.str() << " edge_len" << g_.length(prev_edge));
                pair<int, int> limits = GetPathLimits(*prev_iter, *iter, (int) cur_first_index.edge_position
                                                      , (int) g_.length(prev_edge) - prev_last_index.edge_position);

                if (limits.first == -1) {
                    DEBUG ("Failed to find Path limits");
                    bwa_hits.push_back(cur_sorted_hits);
                    edges.push_back(cur_sorted_edges);
                    cur_sorted_edges.clear();
                    cur_sorted_hits.clear();
                    prev_edge = EdgeId();
                    continue;
                }

                int seq_end = cur_first_index.read_position;
                int seq_start = prev_last_index.read_position;
                if (seq_end < seq_start) {
                    DEBUG ("modifying limits because of some bullshit magic, seq length 0")
                    seq_end = seq_start;
                }
                DEBUG("taking subseq" << seq_start << " " << seq_end << " " << s.size());
                string seq_string = s.Subseq(seq_start, min(seq_end, (int) s.size() )).str();
                GapFillerResult res = gap_filler_.Run(seq_string,
                                                      GraphPosition(prev_edge, prev_last_index.edge_position),
                                                      GraphPosition(cur_edge, cur_first_index.edge_position),
                                                      limits.first, limits.second);
                vector<EdgeId> intermediate_path = res.full_intermediate_path;
                if (res.return_code.status != 0) {
                    bwa_hits.push_back(cur_sorted_hits);
                    edges.push_back(cur_sorted_edges);
                    cur_sorted_edges.clear();
                    cur_sorted_hits.clear();
                    prev_edge = EdgeId();
                    continue;
                } else {
                    cur_sorted_edges.pop_back();
                    for (EdgeId edge : intermediate_path) {
                        cur_sorted_edges.push_back(edge);
                    }
                    cur_sorted_edges.pop_back();
                    if (intermediate_path.size() == 1) {
                        DEBUG("Merged e1=" << prev_edge.int_id() << " e2=" << cur_edge.int_id()
                              << " s=" << prev_last_index.edge_position  << " e=" << cur_first_index.edge_position);
                    }
                }
            }
        }
        MappingInstance cur_first_index = iter->sorted_positions[iter->first_trustable_index];
        MappingInstance cur_last_index = iter->sorted_positions[iter->last_trustable_index];
        cur_sorted_edges.push_back(cur_edge);
        cur_sorted_hits.push_back(cur_edge, omnigraph::MappingRange(
                                      Range(cur_first_index.read_position, cur_last_index.read_position),
                                      Range(cur_first_index.edge_position, cur_last_index.edge_position)));
        prev_edge = cur_edge;
        ++iter;
    }
    if (cur_sorted_edges.size() > 0) {
        edges.push_back(cur_sorted_edges);
        bwa_hits.push_back(cur_sorted_hits);
    }
}


pair<int, int> GAligner::GetPathLimits(const QualityRange &a,
                                       const QualityRange &b,
                                       int s_add_len, int e_add_len) const {
    int start_pos = a.sorted_positions[a.last_trustable_index].read_position;
    int end_pos = b.sorted_positions[b.first_trustable_index].read_position;
    int seq_len = -start_pos + end_pos;
//TODO::something more reasonable
    int path_min_len = max(int(floor((seq_len - int(g_.k())) * pb_config_.path_limit_pressing)), 0);
    int path_max_len = (int) ((double) (seq_len + g_.k() * 2) * pb_config_.path_limit_stretching);
    if (seq_len < 0) {
        DEBUG("suspicious negative seq_len " << start_pos << " " << end_pos << " " << path_min_len << " " << path_max_len);
        if (path_max_len < 0)
            return make_pair(-1, -1);
    }
    path_min_len = max(path_min_len - int(s_add_len + e_add_len), 0);
    path_max_len = max(path_max_len - int(s_add_len + e_add_len), 0);
    return make_pair(path_min_len, path_max_len);
}

bool GAligner::TopologyGap(EdgeId first, EdgeId second, bool oriented) const {
    omnigraph::TerminalVertexCondition<Graph> condition(g_);
    bool res = condition.Check(g_.EdgeEnd(first)) && condition.Check(g_.EdgeStart(second));
    if (!oriented)
        res |= condition.Check(g_.EdgeStart(first)) && condition.Check(g_.EdgeEnd(second));
    return res;
}


OneReadMapping GAligner::GetReadAlignment(const io::SingleRead &read) const {
    auto paths  = pac_index_.GetChainingPaths(read);
    size_t len = paths.size();
    vector<vector<debruijn_graph::EdgeId> > sorted_edges;
    vector<omnigraph::MappingPath<debruijn_graph::EdgeId>> sorted_bwa_hits;
    Sequence s = read.sequence();
    omnigraph::MappingRange read_range;
    vector<bool> block_gap_closer;
    vector<QualityRange> start_clusters, end_clusters;
    vector<int> used(len);
    for (size_t i = 0; i < len; i++) {
        ProcessCluster(s, paths[i], start_clusters, end_clusters, sorted_edges, sorted_bwa_hits,
                       block_gap_closer);

    }
    vector<PathRange> read_ranges;
    int max_path = 0;
    size_t max_path_ind = 0;
    for (size_t i = 0; i < sorted_bwa_hits.size(); ++ i) {
        if ((int)(sorted_bwa_hits[i].mapping_at(sorted_bwa_hits[i].size() - 1).initial_range.end_pos -
                  sorted_bwa_hits[i].mapping_at(0).initial_range.start_pos) > max_path) {
            max_path = (int) (sorted_bwa_hits[i].mapping_at(sorted_bwa_hits[i].size() - 1).initial_range.end_pos -
                              sorted_bwa_hits[i].mapping_at(0).initial_range.start_pos);
            max_path_ind = i;
        }
    }
    int shortest_len = (int) s.size() / 2;
    for (size_t i = 0; i < sorted_bwa_hits.size(); ++ i) {
        PathRange cur_range(MappingPoint(sorted_bwa_hits[i].mapping_at(0).initial_range.start_pos
                                         , sorted_bwa_hits[i].mapping_at(0).mapped_range.start_pos),
                            MappingPoint(sorted_bwa_hits[i].mapping_at(sorted_bwa_hits[i].size() - 1).initial_range.end_pos
                                         , sorted_bwa_hits[i].mapping_at(sorted_bwa_hits[i].size() - 1).mapped_range.end_pos));
        if (restore_ends_  && max_path_ind == i &&
                (sorted_edges.size() == 1 || max_path > min(500, shortest_len) )) {
            unsigned return_code = RestoreEndsF(s, (int) s.size(), sorted_edges[i], cur_range);
            for (size_t j = sorted_bwa_hits.size() - 1; j > i  && return_code != 0; -- j) {
                int end = (int) sorted_bwa_hits[j].mapping_at(0).initial_range.start_pos;
                if (end > (int) cur_range.path_end.seq_pos) {
                    return_code = RestoreEndsF(s, end, sorted_edges[i], cur_range);
                }
            }
            return_code = RestoreEndsB(s, 0, sorted_edges[i], cur_range);
            if (i > 0) {
                for (size_t j = 0; j < i && return_code != 0; ++ j) {
                    int start = (int) sorted_bwa_hits[j].mapping_at(sorted_bwa_hits[j].size() - 1).initial_range.end_pos;
                    if (start < (int) cur_range.path_start.seq_pos) {
                        return_code = RestoreEndsB(s, start, sorted_edges[i], cur_range);
                    }
                }
            }
        }
        read_ranges.push_back(cur_range);
    }
    return AddGapDescriptions(start_clusters, end_clusters, sorted_edges, sorted_bwa_hits, read_ranges, s,
                              block_gap_closer);
}

int GAligner::RestoreEndsB(const Sequence &s,
                           int start,
                           vector<debruijn_graph::EdgeId> &sorted_edges,
                           PathRange &cur_range) const {
    bool forward = true;
    GraphPosition start_pos(sorted_edges[0], cur_range.path_start.edge_pos);
    Sequence ss = s.Subseq(start, cur_range.path_start.seq_pos);
    GapFillerResult res_backward = gap_filler_.Run(ss, start_pos, !forward, sorted_edges, cur_range);
    DEBUG("Backward return_code_ends=" << res_backward.return_code.status)
    return res_backward.return_code.status;
}

int GAligner::RestoreEndsF(const Sequence &s,
                           int end,
                           vector<debruijn_graph::EdgeId> &sorted_edges,
                           PathRange &cur_range) const {
    bool forward = true;
    GraphPosition end_pos(sorted_edges[sorted_edges.size() - 1]
                          , cur_range.path_end.edge_pos);
    Sequence ss = s.Subseq(cur_range.path_end.seq_pos, end);
    GapFillerResult res_forward = gap_filler_.Run(ss, end_pos, forward, sorted_edges, cur_range);
    DEBUG("Forward return_code_ends=" << res_forward.return_code.status)
    return res_forward.return_code.status;
}

void GAligner::ProcessCluster(const Sequence &s,
                              vector<QualityRange> &cur_cluster,
                              vector<QualityRange> &start_clusters,
                              vector<QualityRange> &end_clusters,
                              vector<vector<debruijn_graph::EdgeId> > &sorted_edges,
                              vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
                              vector<bool> &block_gap_closer) const {
    sort(cur_cluster.begin(), cur_cluster.end(),
    [](const QualityRange & a, const QualityRange & b) {
        return (a.average_read_position < b.average_read_position);
    });
    VERIFY(cur_cluster.size() > 0);
    auto cur_cluster_start = cur_cluster.begin();
    auto cur_cluster_end = cur_cluster.end() - 1;
    vector<vector<debruijn_graph::EdgeId> > edges;
    vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > bwa_hits;
    FillGapsInCluster(cur_cluster, s, edges, bwa_hits);
    for (auto &cur_sorted : edges) {
        DEBUG("Adding " << edges.size() << " subreads, cur alignments " << cur_sorted.size());
        if (cur_sorted.size() > 0) {
            for (EdgeId eee : cur_sorted) {
                DEBUG (g_.int_id(eee));
            }
            start_clusters.push_back(*cur_cluster_start);
            end_clusters.push_back(*cur_cluster_end);
            sorted_edges.push_back(cur_sorted);
            //Blocking gap closing inside clusters;
            block_gap_closer.push_back(true);
        }
    }
    for (auto &cur_sorted : bwa_hits) {
        if (cur_sorted.size() > 0) {
            sorted_bwa_hits.push_back(cur_sorted);
        }
    }
    //quite strange if bwa hits were not connected because of consistency check failed
    if (block_gap_closer.size() > 0)
        block_gap_closer[block_gap_closer.size() - 1] = false;

}

OneReadMapping GAligner::AddGapDescriptions(const vector<QualityRange> &start_clusters,
        const vector<QualityRange> &end_clusters,
        const vector<vector<debruijn_graph::EdgeId> > &sorted_edges,
        const vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
        const vector<PathRange> &read_ranges,
        const Sequence &s,
        const vector<bool> &block_gap_closer) const {
    DEBUG("adding gaps between subreads");
    vector<GapDescription> illumina_gaps;
    for (size_t i = 0; i + 1 < sorted_edges.size(); i++) {
        if (block_gap_closer[i])
            continue;
        size_t j = i + 1;
        EdgeId before_gap = sorted_edges[i][sorted_edges[i].size() - 1];
        EdgeId after_gap = sorted_edges[j][0];
//do not add "gap" for rc-jumping
        if (before_gap != after_gap && before_gap != g_.conjugate(after_gap)) {
            if (TopologyGap(before_gap, after_gap, true)) {
                if (start_clusters[j].CanFollow(end_clusters[i])) {
                    const auto &a = end_clusters[i];
                    const auto &b = start_clusters[j];
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
                        DEBUG("adding gap between alignments number " << i << " and " << j);
                    }
                }

            }
        }
    }
    DEBUG("Resulting hits num=" << sorted_edges.size());
    return OneReadMapping(sorted_edges, sorted_bwa_hits, illumina_gaps, read_ranges);
}
}