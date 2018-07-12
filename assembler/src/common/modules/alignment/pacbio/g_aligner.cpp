#include "modules/alignment/pacbio/g_aligner.hpp"
#include "modules/alignment/pacbio/gap_filler.hpp"
#include "modules/alignment/pacbio/pacbio_read_structures.hpp"

namespace sensitive_aligner {

    using debruijn_graph::EdgeId;
    using debruijn_graph::VertexId;
    using debruijn_graph::Graph;
    void GAligner::FillGapsInCluster(const vector<QualityRange> &cur_cluster,
                           const Sequence &s,
                           std::vector<vector<debruijn_graph::EdgeId> > &edges,
                           std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &bwa_hits) const {
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
                double stretched_graph_len = (prev_edge.int_id() != cur_edge.int_id()) ?
                                             (double) (cur_first_index.edge_position + g_.k()) +
                                             ((int) g_.length(prev_edge) - prev_last_index.edge_position) *
                                             pb_config_.path_limit_stretching
                                                                                       :
                                             ((int) cur_first_index.edge_position -
                                              (int) prev_last_index.edge_position) * pb_config_.path_limit_stretching;
                if ((start_v != end_v || (start_v == end_v && read_gap_len > stretched_graph_len))
                    && (prev_edge.int_id() != cur_edge.int_id() ||
                        (prev_edge.int_id() == cur_edge.int_id() && stretched_graph_len < 0) ||
                        (prev_edge.int_id() == cur_edge.int_id() && stretched_graph_len > 0 &&
                         read_gap_len > stretched_graph_len))) {
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
                    string s_add = "";
                    string e_add = "";
                    int seq_end = cur_first_index.read_position;
                    int seq_start = prev_last_index.read_position;
                    string tmp = g_.EdgeNucls(prev_edge).str();
                    s_add = tmp.substr(prev_last_index.edge_position,
                                       g_.length(prev_edge) - prev_last_index.edge_position);
                    tmp = g_.EdgeNucls(cur_edge).str();
                    e_add = tmp.substr(0, cur_first_index.edge_position);
                    pair<int, int> limits = pac_index_.GetPathLimits(*prev_iter, *iter, (int) s_add.length(),
                                                          (int) e_add.length());
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
                    DEBUG("taking subseq" << seq_start << " " << end_pos << " " << s.size());
                    std::string seq_string = s.Subseq(seq_start, min(end_pos, s_len)).str();
                    GapFiller gap_filler(g_, pb_config_, gap_cfg_);
                    GapFillerResult res = gap_filler.Run(seq_string,
                                                         GraphPosition(prev_edge, prev_last_index.edge_position),
                                                         GraphPosition(cur_edge, cur_first_index.edge_position),
                                                         limits.first, limits.second);
                    vector<EdgeId> intermediate_path = res.intermediate_path;
                    if (intermediate_path.size() == 0) {
//                        DEBUG(DebugEmptyBestScoredPath(start_v, end_v, prev_edge, cur_edge,
//                                                       prev_last_index.edge_position, cur_first_index.edge_position,
//                                                       seq_end - seq_start));
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

    bool GAligner::TopologyGap(EdgeId first, EdgeId second, bool oriented) const {
        omnigraph::TerminalVertexCondition<Graph> condition(g_);
        bool res = condition.Check(g_.EdgeEnd(first)) && condition.Check(g_.EdgeStart(second));
        if (!oriented)
            res |= condition.Check(g_.EdgeStart(first)) && condition.Check(g_.EdgeEnd(second));
        return res;
    }


    OneReadMapping GAligner::GetReadAlignment(const io::SingleRead &read) const {
        std::vector<ColoredRange> colors = pac_index_.GetRangedColors(read);
        size_t len = colors.size();
        std::vector<vector<debruijn_graph::EdgeId> > sorted_edges;
        vector<omnigraph::MappingPath<debruijn_graph::EdgeId>> sorted_bwa_hits;
        Sequence s = read.sequence();
        omnigraph::MappingRange read_range;
        vector<bool> block_gap_closer;
        vector<QualityRange> start_clusters, end_clusters;
        vector<int> used(len);
        for (size_t i = 0; i < len; i++) {
            used[i] = 0;
            DEBUG(colors[i].second << " " << colors[i].first.str(g_));
        }
        //FIXME ferther code is AWFUL
        for (size_t i = 0; i < len; i++) {
            int cur_color = colors[i].second;
            if (!used[i] && cur_color != DELETED_COLOR) {
                DEBUG("starting new subread");
                vector<QualityRange> cur_cluster;
                used[i] = 1;
                for (size_t j = 0; j < len; j++) {
                    if (colors[j].second == cur_color) {
                        cur_cluster.push_back(colors[j].first);
                        used[j] = 1;
                    }
                }
                ProcessCluster(s, cur_cluster, start_clusters, end_clusters, sorted_edges, sorted_bwa_hits,
                               block_gap_closer);
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
            EndsFiller ends_filler(g_, pb_config_, gap_cfg_);
            ends_filler.Run(sorted_bwa_hits[0], sorted_edges[0], s, !forward, cur_range, return_code);
            DEBUG("Backward return_code_ends=" << return_code)
            ends_filler.Run(sorted_bwa_hits[0], sorted_edges[0], s, forward, cur_range, return_code);
            DEBUG("Forward return_code_ends=" << return_code)
            read_ranges.push_back(cur_range);
        } else {
            for (auto hits: sorted_bwa_hits) {
                PathRange cur_range;
                cur_range.seq_start = hits.mapping_at(0).initial_range.start_pos;
                cur_range.seq_end = hits.mapping_at(hits.size() - 1).initial_range.end_pos;
                cur_range.edge_start = hits.mapping_at(0).mapped_range.start_pos;
                cur_range.edge_end = hits.mapping_at(hits.size() - 1).mapped_range.end_pos;
                read_ranges.push_back(cur_range);
            }
        }
        return AddGapDescriptions(start_clusters, end_clusters, sorted_edges, sorted_bwa_hits, read_ranges, s,
                                  block_gap_closer);
    }

    void GAligner::ProcessCluster(const Sequence &s,
                             std::vector<QualityRange> &cur_cluster,
                             std::vector<QualityRange> &start_clusters,
                             std::vector<QualityRange> &end_clusters,
                             std::vector<vector<debruijn_graph::EdgeId> > &sorted_edges,
                             std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
                             std::vector<bool> &block_gap_closer) const {
        std::sort(cur_cluster.begin(), cur_cluster.end(),
                  [](const QualityRange &a, const QualityRange &b) {
                      return (a.average_read_position < b.average_read_position);
                  });
        VERIFY(cur_cluster.size() > 0);
        auto cur_cluster_start = cur_cluster.begin();
        for (auto iter = cur_cluster.begin(); iter != cur_cluster.end(); ++iter) {
            auto next_iter = iter + 1;
            if (next_iter == cur_cluster.end() || !pac_index_.IsConsistent(*iter, *next_iter)) {
                if (next_iter != cur_cluster.end()) {
                    DEBUG("clusters splitted:");
                    DEBUG("on " << iter->str(g_));
                    DEBUG("and " << next_iter->str(g_));
                }
                vector<QualityRange> splitted_cluster(cur_cluster_start, next_iter);
                std::vector<vector<debruijn_graph::EdgeId> > edges;
                std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > bwa_hits;
                FillGapsInCluster(splitted_cluster, s, edges, bwa_hits);
                for (auto &cur_sorted: edges) {
                    DEBUG("Adding " << edges.size() << " subreads, cur alignments " << cur_sorted.size());
                    if (cur_sorted.size() > 0) {
                        for (EdgeId eee: cur_sorted) {
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
                DEBUG("on " << iter->str(g_));
                DEBUG("and " << next_iter->str(g_));
            }
        }
    }

    OneReadMapping GAligner::AddGapDescriptions(const std::vector<QualityRange> &start_clusters,
                                                const std::vector<QualityRange> &end_clusters,
                                                const std::vector<vector<debruijn_graph::EdgeId> > &sorted_edges,
                                                const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
                                                const std::vector<PathRange> &read_ranges,
                                                const Sequence &s,
                                                const std::vector<bool> &block_gap_closer) const {
        DEBUG("adding gaps between subreads");
        std::vector<GapDescription> illumina_gaps;
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