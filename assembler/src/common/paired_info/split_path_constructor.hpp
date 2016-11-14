//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
* split_path_constructor.hpp
*
*  Created on: Jun 14, 2012
*      Author: avsirotkin
*/

#pragma once

#include <common/assembly_graph/paths/path_utils.hpp>
#include "utils/logger/logger.hpp"
#include "paired_info/paired_info.hpp"
#include "assembly_graph/paths/path_processor.hpp"
#include "paired_info/pair_info_bounds.hpp"

namespace debruijn_graph {

template<class Graph>
class PathInfoClass {
public:
    typedef typename Graph::EdgeId EdgeId;
    typedef omnigraph::de::PairInfo<EdgeId> PairInfo;

    EdgeId base_edge;
    vector<PairInfo> path;

    PathInfoClass() : base_edge(NULL) { };

    PathInfoClass(const EdgeId Edge) : base_edge(Edge) { };

    std::pair<EdgeId, double> operator[](const size_t i) const {
        if (i == 0)
            return std::make_pair(base_edge, 0.0);

        VERIFY(i < path.size() + 1);
        return std::make_pair(path[i - 1].second, path[i - 1].d());
    }

    size_t size() const { return path.size() + 1; }

    void push_back(const PairInfo &pi) { path.push_back(pi); }

    typename std::vector<PairInfo>::const_iterator begin() const { return path.begin(); }

    typename std::vector<PairInfo>::const_iterator end() const { return path.end(); }

    std::string PrintPath(const Graph &graph) const {
        std::ostringstream ss;
        ss << " " << graph.int_id(base_edge) << ": ";
        for (size_t j = 0; j < path.size(); j++) {
            ss << "(" << graph.int_id(path[j].second) << ", " << path[j].d() << "), ";
        }
        return ss.str();
    }
};

template<class Graph>
class SplitPathConstructor {
    typedef typename Graph::EdgeId EdgeId;
    typedef PathInfoClass<Graph> PathInfo;
    typedef omnigraph::de::PairInfo<EdgeId> PairInfo;
    static const size_t MAX_DIJKSTRA_DEPTH = 3000;
public:
    SplitPathConstructor(const Graph &graph) : graph_(graph) { }

    vector<PathInfo> ConvertPIToSplitPaths(EdgeId cur_edge, const omnigraph::de::PairedInfoIndexT<Graph> &pi,
                                           double is, double is_var) const {
        vector<PairInfo> pair_infos; //TODO: this is an adaptor for the old implementation
        for (auto i : pi.Get(cur_edge))
            for (auto j : i.second)
                pair_infos.emplace_back(cur_edge, i.first, j);
        std::sort(pair_infos.begin(), pair_infos.end(),[&](const PairInfo p1, const PairInfo p2){
            return (p1.point.d > p2.point.d || ((p1.point.d == p2.point.d) && (p1.second < p2.second )));
        });
        vector<PathInfo> result;
        if (pair_infos.empty())
            return result;

        vector<bool> pair_info_used(pair_infos.size());
        TRACE("Preparing path_processor for this base edge");
        size_t path_upper_bound = PairInfoPathLengthUpperBound(graph_.k(), (size_t) is, is_var);

        //FIXME is path_upper_bound enough?


        typename omnigraph::DijkstraHelper<Graph>::BoundedDijkstra dijkstra(omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(graph_, path_upper_bound, MAX_DIJKSTRA_DEPTH));
        dijkstra.Run(graph_.EdgeEnd(cur_edge));
        for (size_t i = pair_infos.size(); i > 0; --i) {
            const PairInfo &cur_info = pair_infos[i - 1];
            if (math::le(cur_info.d(), 0.))
                continue;
            if (pair_info_used[i - 1])
                continue;
            DEBUG("SPC: pi " << cur_info);

            vector<EdgeId> common_part = GetCommonPathsEnd(graph_, cur_edge, cur_info.second,
                                                           (size_t) (cur_info.d() - cur_info.var()),
                                                           (size_t) (cur_info.d() + cur_info.var()),
                                                           dijkstra);
            DEBUG("Found common part of size " << common_part.size());
            PathInfoClass<Graph> sub_res(cur_edge);
            if (common_part.size() > 0) {
                size_t total_length = 0;
                for (size_t j = 0; j < common_part.size(); ++j)
                    total_length += graph_.length(common_part[j]);

                DEBUG("Common part " << ToString(common_part));
                for (size_t j = 0; j < common_part.size(); ++j) {
                    PairInfo cur_pi(cur_edge, common_part[j],
                                    cur_info.d() - (double) total_length,
                                    cur_info.weight(),
                                    cur_info.var());

                    sub_res.push_back(cur_pi);
                    total_length -= graph_.length(common_part[j]);
                    for (size_t ind = 0; ind + 1 < i; ++ind) {
                        if (cur_pi == pair_infos[ind])
                            pair_info_used[ind] = true;
                    }
                }
            }

            sub_res.push_back(cur_info);
            result.push_back(sub_res);
            DEBUG(sub_res.PrintPath(graph_));
        }
        return result;
    }

private:
    const Graph &graph_;
};


}
