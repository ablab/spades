//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * path_utils.hpp
 *
 */

#pragma once

#include "path_processor.hpp"
#include "mapping_path.hpp"
#include "assembly_graph/dijkstra/dijkstra_algorithm.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"

#include "sequence/sequence.hpp"
#include "sequence/sequence_tools.hpp"

namespace debruijn_graph {

    template<class Graph>
    std::vector<typename Graph::EdgeId> GetCommonPathsEnd(
            const Graph &g,
            typename Graph::EdgeId e1,
            typename Graph::EdgeId e2,
            size_t min_dist,
            size_t max_dist,
            typename omnigraph::DijkstraHelper<Graph>::BoundedDijkstra &dijkstra) {

        typedef typename Graph::EdgeId EdgeId;
        typedef typename Graph::VertexId VertexId;

        std::vector<EdgeId> res;
        VERIFY (min_dist >= g.length(e1));
        VERIFY (max_dist >= g.length(e1));
        size_t dist = max_dist - g.length(e1);
        VertexId cur_vertex = g.EdgeStart(e2);
        if (!dijkstra.DistanceCounted(cur_vertex))
            return res;
        size_t cur_dist;
        if ((cur_dist = dijkstra.GetDistance(cur_vertex)) > dist)
            return res;
        size_t suffix_len = 0;
        while (cur_dist > 0) {
            EdgeId prev_edge;
            bool found = false;
            for (auto edge: g.IncomingEdges(cur_vertex)) {
                if ((dijkstra.DistanceCounted(g.EdgeStart(edge))) && (
                        suffix_len + g.length(edge) + dijkstra.GetDistance(g.EdgeStart(edge)) <= dist)) {
                    if (found) {
                        std::reverse(res.begin(), res.end());
                        return res;
                    } else {
                        found = true;
                        prev_edge = edge;
                    }
                }
            }
            if (!found)
                return res;
            else {
                suffix_len += g.length(prev_edge);
                VERIFY(cur_dist >= g.length(prev_edge));
                cur_dist -= g.length(prev_edge);
                cur_vertex = g.EdgeStart(prev_edge);
                res.push_back(prev_edge);
            }
        }
        std::reverse(res.begin(), res.end());
        return res;
    }

    template<class Graph>
    std::vector<std::vector<typename Graph::EdgeId>> GetAllPathsBetweenEdges(
            const Graph &g,
            typename Graph::EdgeId &e1,
            typename Graph::EdgeId &e2, size_t min_dist,
            size_t max_dist) {
        omnigraph::PathStorageCallback<Graph> callback(g);
        ProcessPaths(g,
                     min_dist,
                     max_dist,
                     g.EdgeEnd(e1), g.EdgeStart(e2),
                     callback);
        auto paths = callback.paths();
        return paths;
    }

    template<class graph_pack>
    size_t GetAllPathsQuantity(const graph_pack &origin_gp,
                               const typename graph_pack::graph_t::EdgeId &e1,
                               const typename graph_pack::graph_t::EdgeId &e2, double d, double is_var) {
        omnigraph::PathStorageCallback<typename graph_pack::graph_t> callback(origin_gp.g);
        omnigraph::PathProcessor<typename graph_pack::graph_t>
                path_processor(origin_gp.g,
                               (size_t) d - origin_gp.g.length(e1) - size_t(is_var),
                               (size_t) d - origin_gp.g.length(e1) + size_t(is_var),
                               origin_gp.g.EdgeEnd(e1),
                               origin_gp.g.EdgeStart(e2),
                               callback);
        path_processor.Process();
        auto paths = callback.paths();
        TRACE(e1.ind_id() << " " << e2.int_id() << " " << paths.size());
        return paths.size();
    }

    template<class Graph>
    Sequence MergeSequences(const Graph &g, const std::vector<typename Graph::EdgeId> &continuous_path) {
        std::vector<Sequence> path_sequences;
        std::vector<uint32_t> overlaps;
        if (continuous_path.size() == 0)
            return Sequence();
        path_sequences.push_back(g.EdgeNucls(continuous_path[0]));
        for (size_t i = 1; i < continuous_path.size(); ++i) {
            auto end = g.EdgeEnd(continuous_path[i - 1]);
            auto start = g.EdgeStart(continuous_path[i]);
            VERIFY(start == end);
            overlaps.push_back(g.data(start).overlap());
            path_sequences.push_back(g.EdgeNucls(continuous_path[i]));
        }
        return MergeOverlappingSequences(path_sequences, overlaps, /*no need to check again*/false);
    }

    template<class Graph>
    Sequence PathSequence(const Graph &g, const omnigraph::Path<typename Graph::EdgeId> &path) {
        Sequence path_sequence = MergeSequences(g, path.sequence());
        size_t start = path.start_pos();
        size_t end = path_sequence.size() - g.length(path[path.size() - 1]) + path.end_pos();
        return path_sequence.Subseq(start, end);
    }

    template<class Graph>
    size_t PathLength(const Graph &g, const omnigraph::Path<typename Graph::EdgeId> &p) {
        if (p.size() == 0)
            return 0;
        const auto &edges = p.sequence();
        size_t len = omnigraph::CumulativeLength(g, edges);
        len -= p.start_pos();
        len -= g.length(edges.back()) - p.end_pos();
        return len;
    }

}
