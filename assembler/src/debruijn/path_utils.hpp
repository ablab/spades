//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * path_utils.hpp
 *
 */

#pragma once

namespace debruijn_graph {

  // TODO: rewrite this function
  template<class Graph>
    vector<typename Graph::EdgeId> GetCommonPathsEnd(
        const Graph& g,
        typename Graph::EdgeId e1,
        typename Graph::EdgeId e2,
        size_t min_dist,
        size_t max_dist,
        PathProcessor<Graph>& path_processor) 
  {
      typedef typename Graph::EdgeId EdgeId;
      typedef vector<EdgeId> Path;

      PathStorageCallback<Graph> callback(g);
      //PathProcessor<Graph> path_processor(g,
                                          //min_dist - g.length(e1),
                                          //max_dist - g.length(e1),
          //g.EdgeEnd(e1), g.EdgeStart(e2), callback);

      path_processor.SetMinLens({min_dist - g.length(e1)});
      path_processor.SetMaxLen(max_dist - g.length(e1));
      path_processor.SetEndPoints({g.EdgeStart(e2)});
      path_processor.SetCallback(&callback);
      path_processor.ResetCallCount();
      int error_code = path_processor.Process();
      vector<Path> paths = callback.paths();

      vector<EdgeId> result;
      if (error_code != 0) {
        DEBUG("Edge " << g.int_id(e1) << " path_processor problem")
        return result;
      }
      if (paths.size() == 0)
        return result;
      if (paths.size() == 1)
        return paths[0];
      size_t j = 0;
      while (j < paths[0].size()) {
        for (size_t i = 1;  i < paths.size(); ++i) {
          if (j == paths[i].size()) {
            vector<EdgeId> result(paths[0].begin()+(paths[0].size() - j), paths[0].end());
            return result;
          } else {
            if (paths[0][paths[0].size()-1-j] != paths[i][paths[i].size()-1-j]) {
              vector<EdgeId> result(paths[0].begin()+(paths[0].size() - j), paths[0].end());
              return result;
            }
          }
        }
        ++j;
      }
      return paths[0];
    }



  template<class Graph>
    vector<vector<typename Graph::EdgeId> > GetAllPathsBetweenEdges(
        const Graph& g,
        typename Graph::EdgeId& e1,
        typename Graph::EdgeId& e2, size_t min_dist,
        size_t max_dist) {
      PathStorageCallback<Graph> callback(g);
      PathProcessor<Graph> path_processor(g,
          min_dist,
          max_dist, //0, *cfg::get().ds.IS - K + size_t(*cfg::get().ds.is_var),
          g.EdgeEnd(e1), g.EdgeStart(e2),
          callback);
      path_processor.Process();
      auto paths = callback.paths();
      return paths;
    }

template<class graph_pack>
size_t GetAllPathsQuantity(const graph_pack& origin_gp,
                           const typename graph_pack::graph_t::EdgeId& e1,
                           const typename graph_pack::graph_t::EdgeId& e2, double d, double is_var) {
  PathStorageCallback<typename graph_pack::graph_t> callback(origin_gp.g);
  PathProcessor<typename graph_pack::graph_t>
      path_processor(origin_gp.g,
                     d - origin_gp.g.length(e1)
                     - size_t(is_var),
                     d - origin_gp.g.length(e1)
                     + size_t(is_var),
                     origin_gp.g.EdgeEnd(e1), 
                     origin_gp.g.EdgeStart(e2),
                     callback);
  path_processor.Process();
  auto paths = callback.paths();
  TRACE(origin_gp.int_ids.ReturnIntId(e1) << " "
        << origin_gp.int_ids.ReturnIntId(e2) << " "
        << paths.size());
  return paths.size();
}

  template<class graph_pack>
    void GenerateMatePairStats(const graph_pack& origin_gp,
        const PairedInfoIndexT<typename graph_pack::graph_t>& clustered_index, double is_var)
    {
      typedef typename graph_pack::graph_t graph_t;
      typedef typename graph_t::EdgeId EdgeId;
      typedef set<Point> Histogram;

      map<size_t, size_t> sizes;
      for (auto e_iter = origin_gp.g.ConstEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
        EdgeId e1 = *e_iter;
        const InnerMap<Graph>& pi = clustered_index.GetEdgeInfo(*e_iter, 0);
        for (auto ext_iter = pi.begin(); ext_iter != pi.end(); ++ext_iter) {
          EdgeId e2 = ext_iter->first;
          const Histogram& hist = ext_iter->second;
          for (auto i_iter = hist.begin(); i_iter != hist.end(); ++i_iter) {
            if (math::ge(i_iter->d, (double) origin_gp.g.length(e1))) 
            {
              size_t tmp = GetAllPathsQuantity(origin_gp, e1, e2, i_iter->d, is_var);
              if (sizes.find(tmp) == sizes.end())
                sizes.insert(make_pair(tmp, 0));
              ++sizes[tmp];
            }
          }
        }
      }
      INFO("Pathset mate pair statistics:");
      for (auto s_iter = sizes.begin(); s_iter != sizes.end(); s_iter++) {
        INFO("- size: " << s_iter->first << "; pathsets: " << s_iter->second);
      }
    }

}
