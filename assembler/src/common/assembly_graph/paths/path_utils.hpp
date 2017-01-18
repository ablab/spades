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

#include "sequence/sequence.hpp"
#include "path_processor.hpp"
#include "mapping_path.hpp"
#include "assembly_graph/dijkstra/dijkstra_algorithm.hpp"

namespace debruijn_graph {


    template<class Graph>
    vector<typename Graph::EdgeId> GetCommonPathsEnd(
            const Graph &g,
            typename Graph::EdgeId e1,
            typename Graph::EdgeId e2,
            size_t min_dist,
            size_t max_dist,
            typename omnigraph::DijkstraHelper<Graph>::BoundedDijkstra &dijkstra) {

        typedef typename Graph::EdgeId EdgeId;
        typedef typename Graph::VertexId VertexId;

        vector<EdgeId> res;
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
            EdgeId prev_edge(0);
            bool found = false;
            for (auto edge: g.IncomingEdges(cur_vertex)) {
                if ((dijkstra.DistanceCounted(g.EdgeStart(edge))) && (
                        suffix_len + g.length(edge) + dijkstra.GetDistance(g.EdgeStart(edge)) <= dist)) {
                    if (found == true) {
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
    vector<vector<typename Graph::EdgeId> > GetAllPathsBetweenEdges(
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
    Sequence MergeSequences(const Graph &g,
                            const vector<typename Graph::EdgeId> &continuous_path) {
        vector<Sequence> path_sequences;
        path_sequences.push_back(g.EdgeNucls(continuous_path[0]));
        for (size_t i = 1; i < continuous_path.size(); ++i) {
            VERIFY(g.EdgeEnd(continuous_path[i - 1]) == g.EdgeStart(continuous_path[i]));
            path_sequences.push_back(g.EdgeNucls(continuous_path[i]));
        }
        return MergeOverlappingSequences(path_sequences, g.k());
    }

    template<class Graph>
    Sequence PathSequence(const Graph &g, const omnigraph::Path<typename Graph::EdgeId> &path) {
        Sequence path_sequence = MergeSequences(g, path.sequence());
        size_t start = path.start_pos();
        size_t end = path_sequence.size() - g.length(path[path.size() - 1]) + path.end_pos();
        return path_sequence.Subseq(start, end);
    }


}
