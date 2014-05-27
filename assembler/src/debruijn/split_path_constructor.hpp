//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * split_path_constructor.hpp
 *
 *  Created on: Jun 14, 2012
 *      Author: avsirotkin
 */

#pragma once

#include "standard.hpp"
#include "logger/logger.hpp"
#include "de/paired_info.hpp"
#include "omni/path_processor.hpp"

namespace debruijn_graph {

template<class Graph>
class PathInfoClass {
  public:
    typedef typename Graph::EdgeId EdgeId;
    typedef omnigraph::de::PairInfo<EdgeId> PairInfo;

    EdgeId base_edge;
    vector<PairInfo> path;
    PathInfoClass(): base_edge(NULL) {};
    PathInfoClass(const EdgeId Edge): base_edge(Edge) {};
    std::pair<EdgeId, double> operator[](const size_t i) const {
        if (i == 0)
            return std::make_pair(base_edge, 0.0);

        VERIFY(i < path.size() + 1);
        return std::make_pair(path[i-1].second, path[i-1].d());
    }
    size_t size() const { return path.size() + 1; }
    void push_back(const PairInfo& pi) { path.push_back(pi); }
    typename std::vector<PairInfo>::const_iterator begin() const { return path.begin(); }
    typename std::vector<PairInfo>::const_iterator end() const { return path.end(); }
    std::string PrintPath(const Graph& graph) const {
        std::ostringstream ss;
        ss<<" "<<graph.int_id(base_edge)<<": ";
        for (size_t j=0; j < path.size(); j++){
            ss<<"("<<graph.int_id(path[j].second)<<", "<<path[j].d()<<"), ";
        }
        return ss.str();
    }
};

template<class Graph>
class SplitPathConstructor {
    typedef typename Graph::EdgeId EdgeId;
    typedef PathInfoClass<Graph> PathInfo;
    typedef omnigraph::de::PairInfo<EdgeId> PairInfo;

  public:
    SplitPathConstructor(const Graph &graph): graph_(graph) {}

    vector<PathInfo> ConvertPIToSplitPaths(const std::vector<PairInfo>& pair_infos, double is, double is_var) const {
        vector<PathInfo> result;
        if (pair_infos.size() == 0)
            return result;

        EdgeId cur_edge = pair_infos[0].first;
        vector<bool> pair_info_used(pair_infos.size());
        TRACE("Preparing path_processor for this base edge");
        size_t path_upper_bound = PairInfoPathLengthUpperBound(graph_.k(), (size_t) is, is_var);

        PathStorageCallback<Graph> callback(graph_);
        PathProcessor<Graph> path_processor(graph_,
                                            path_upper_bound,
                                            path_upper_bound,
                                            graph_.EdgeEnd(cur_edge), graph_.EdgeStart(cur_edge), callback);

        TRACE("Path_processor is done");

        for (size_t i = pair_infos.size(); i > 0; --i) {
            const PairInfo& cur_info = pair_infos[i - 1];
            if (math::le(cur_info.d(), 0.))
                continue;
            if (pair_info_used[i - 1])
                continue;
            DEBUG("SPC: pi " << cur_info);
            vector<EdgeId> common_part = GetCommonPathsEnd(graph_, cur_edge, cur_info.second,
                                                           (size_t) (cur_info.d() - cur_info.var()),
                                                           (size_t) (cur_info.d() - cur_info.var()),
                                                           path_processor);
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
