//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <iostream>
#include "coloring.hpp"
#include "coordinates_handler.hpp"
#include "compare_standard.hpp"
#include "comparison_utils.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"

namespace cap {

template <class Graph>
class ReliableTargetedBoundedDijkstra;
template <class Graph>
class GenomePathsFinder;

/*
 *  SimpleInversionFinder searches for inversions occured in a set of genomes
 *  relative to the others (all_genomes minus inversed_genomes).
 *  Currently only for two genome sequences.
 *
 *  Algorithm just searches for "cycles" of length 4 with alternating colors
 *  of the following kind:
 *    v1 --red-> v3
 *    v2 --red-> v4
 *    v1 -blue-> v4
 *    v2 -blue-> v3
 */

template <class gp_t>
class SimpleInversionFinder {
 public:
  SimpleInversionFinder(gp_t &gp, ColorHandler<Graph> &coloring,
      CoordinatesHandler<Graph> &coordinates_handler,
      const std::string base_pic_file_name,
      const bool mask_inversed = false)
      : gp_(gp),
        g_(gp_.g),
        coloring_(coloring),
        coordinates_handler_(coordinates_handler),
        base_pic_file_name_(base_pic_file_name),
        num_cycles_found_(0),
        found_lengths_(),
        mask_inversed_(mask_inversed) {
  }

  void FindInversionEvents() {
    num_cycles_found_ = 0;

    INFO("Searching for inversions");
    for (auto it = g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
      CheckForCycle(*it);
    }
    INFO("Searching for inversions done. Cycles found: " << num_cycles_found_);

    INFO("Found lengths:");
    const std::vector<size_t> found = found_lengths();
    for (auto it = found.begin(); it != found.end(); ++it) {
      INFO("" << *it);
    }
  }

  std::vector<size_t> found_lengths() {
    std::sort(found_lengths_.begin(), found_lengths_.end());
    return found_lengths_;
  }

 private:
  typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

  typedef std::vector<EdgeId> EdgeList;
  typedef std::vector<EdgeId> Path;
  typedef std::vector<VertexId> VertexList;
  typedef uint64_t u64int;


  gp_t &gp_;
  Graph &g_;
  ColorHandler<Graph> &coloring_;
  CoordinatesHandler<Graph> &coordinates_handler_;

  //ostream &output_stream_;
  std::string base_pic_file_name_;
  size_t num_cycles_found_;
  std::vector<size_t> found_lengths_;

  bool mask_inversed_;

  void CheckForCycle(const VertexId v1) {
    // the search assumes colors 0 and 1 are used
    const std::vector<EdgeList> &solo_color_v1_out = GetSoloColoredEdgeLists(
        g_.OutgoingEdges(v1));

    for (auto it4 = solo_color_v1_out[1].begin();
         it4 != solo_color_v1_out[1].end(); ++it4) {
      const VertexId v4 = g_.EdgeEnd(*it4);
      if (v4 < v1) continue;

      const EdgeList back_list = GetSoloColoredEdgeLists(
          g_.IncomingEdges(v4))[0];

      for (auto it2 = back_list.begin(); it2 != back_list.end(); ++it2) {
        const VertexId v2 = g_.EdgeStart(*it2);
        if (v2 == v1 || v2 < v1)
          continue;

        const EdgeList v4_to_v3_list = GetSoloColoredEdgeLists(
            g_.OutgoingEdges(v2))[1];

        VertexList v3_list = FindEqualVerticesOnEnds(
            solo_color_v1_out[0], v4_to_v3_list);

        while (v3_list.size() > 0 && v3_list.back() == v4) {
          v3_list.pop_back();
        }

        if (v3_list.size() > 0) {
          if (v3_list[0] < v1)
            continue;

          DEBUG("Found inversion! (" << v3_list.size() << " candidates)");

          ProcessFoundCycle(VertexList({v1, v2, v3_list[0], v4}));

          // should return here
        }
      }
    }
  }

  void ProcessFoundCycle(const VertexList &cycle) {
    stringstream v_list_str;
    for (auto it = cycle.begin(); it != cycle.end(); ++it) {
      v_list_str << g_.str(*it) << " ";
    }
    DEBUG("cycle found: " << v_list_str.str());

    const std::string edge_pic_name = base_pic_file_name_ + "_" +
        ToString(num_cycles_found_) + ".dot";
    const std::string path_pic_name = base_pic_file_name_ + "_path_" +
        ToString(num_cycles_found_) + ".dot";

    /*
    PrintColoredGraphAroundEdge(g_, coloring_, edge, gp_.edge_pos,
        edge_pic_name);
        */

    ssize_t length = -1;
    ssize_t l1 = FindAndPrintPath(cycle[0], g_.conjugate(cycle[1]), cycle[2], cycle[3]);
    //int l2 = FindAndPrintPath(g_.conjugate(cycle[1]), cycle[0], cycle[2], cycle[3]);
    ssize_t l2 = FindAndPrintPath(g_.conjugate(cycle[2]), cycle[3], g_.conjugate(cycle[0]), g_.conjugate(cycle[1]));
    if (l1 > 0 && (length < 0 || length > l1))
      length = l1;
    if (l2 > 0 && (length < 0 || length > l2))
      length = l2;
      /*
    FindAndPrintPath(cycle[2], g_.conjugate(cycle[3]), path_pic_name) ||
    FindAndPrintPath(g_.conjugate(cycle[2]), cycle[3], path_pic_name) ||
    FindAndPrintPath(cycle[0], g_.conjugate(cycle[1]), path_pic_name) ||
    FindAndPrintPath(g_.conjugate(cycle[0]), cycle[1], path_pic_name);
    */

    if (length < 0) {
      INFO("found cycle but not path!");
      return;
    }
    num_cycles_found_++;
    found_lengths_.push_back(length);
  }

  template<class EdgeContainer>
  std::vector<EdgeList> GetSoloColoredEdgeLists(const EdgeContainer &all_edges) const {
    std::vector<EdgeList> result(coloring_.max_colors());

    for (auto it = all_edges.begin(); it != all_edges.end(); ++it) {
      int color = GetEdgeSoloColor(*it);
      if (IsSoloColor(color)) {
        VERIFY(color < (int)result.size());
        result[color].push_back(*it);
      }
    }

    return result;
  }

  int GetEdgeSoloColor(const EdgeId edge) const {
    const TColorSet &color = coloring_.Color(edge);
    int result = 0;
    for (unsigned i = 0; i < coloring_.max_colors(); ++i) {
      if (!color[i])
        continue;

      if (result) {
        result = -1;
      } else {
        result = i + 1;
      }
    }

    return result - 1;
  }

  inline bool IsSoloColor(int solo_color) const {
    return solo_color >= 0;
  }

  VertexList FindEqualVerticesOnEnds(
      const EdgeList &l1, const EdgeList &l2) const {
    const VertexList &vl1 = GetSortedEdgeEnds(l1);
    const VertexList &vl2 = GetSortedEdgeEnds(l2);

    VertexList result;

    size_t i = 0, j = 0;
    while (i < vl1.size() && j < vl2.size()) {
      if (vl1[i] == vl2[j]) {
        result.push_back(vl1[i]);
        ++i, ++j;
      } else if (vl1[i] < vl2[j]) {
        ++i;
      } else {
        ++j;
      }
    }

    return result;
  }

  VertexList GetSortedEdgeEnds(const EdgeList &l) const {
    VertexList result;
    result.reserve(l.size());

    for (auto it = l.begin(); it != l.end(); ++it) {
      result.push_back(g_.EdgeEnd(*it));
    }

    std::sort(result.begin(), result.end());
    return result;
  }

  inline ssize_t FindAndPrintPath(const VertexId v1, const VertexId v2, const VertexId v3, const VertexId v4) const {
    TRACE("Finding paths from " << g_.str(v1) << " to " << g_.str(v2));
    const static size_t max_length = 800000;

    std::vector<EdgeId> out_edges;
    for (const auto e : g_.OutgoingEdges(v1)) {
      if (g_.EdgeEnd(e) == v3 || g_.EdgeEnd(e) == v4) {
        out_edges.push_back(e);
        //INFO("out edge " << g_.str(e));
      }
    }
    if (out_edges.size() == 0)
      return -1;

    VERIFY(out_edges.size() == 2);

    GenomePathsFinder<Graph> dfs(g_, coordinates_handler_);
    std::vector<Path> paths = dfs.FindGenomePaths(out_edges, v2, max_length);
    if (paths.size() > 0) {
      if (paths.size() == 1) {
        //INFO("found only one path, strange");
      }
      return GetPathLength(paths[0]);
    }
    //INFO("could not find any path :(");

    /*
    ReliableTargetedBoundedDijkstra<Graph> dijkstra(g_, v2, max_length, 500);
    dijkstra.run(v1);
    if (dijkstra.DistanceCounted(v2)) {
      TRACE("Finding path done; length=" << dijkstra.GetDistance(v2));
      return dijkstra.GetDistance(v2);
    }
    */

    return -1;
  }

  size_t GetPathLength(const Path &p) const {
    size_t res = 0;
    for (const auto &edge : p) {
      res += g_.length(edge);
    }
    return res;
  }

  inline void PrintPath(const EdgeList &path, const std::string out_file) const {
    const static size_t edge_length = 1;
    const static size_t max_vertices = 100;
    MappingPath<EdgeId> mpath = TrivialMappingPath(g_, path);
    //Path<EdgeId> cpath(path, mpath.start_pos(), mpath.end_pos());

    visualization::graph_labeler::LengthIdGraphLabeler<Graph> basic_labeler(g_);
    visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(g_, gp_.edge_pos);
    visualization::graph_labeler::CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);

    WriteComponentsAlongPath(g_, labeler, out_file, edge_length, max_vertices,
      mpath, *ConstructBorderColorer(g_, coloring_));
  }

    size_t PathLength(const EdgeList &path) const {
        size_t res = 0;
        for (auto I = path.begin(); I != path.end(); ++I)
            res += g_.length(*I);
        return res;
    }

  DECL_LOGGER("SimpleInversionFinder")
  ;
};

template <class Graph>
class GenomePathsFinder {
 public:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef std::vector<EdgeId> Path;
  typedef typename CoordinatesHandler<Graph>::PosArray PosArray;

  GenomePathsFinder(const Graph &g, const CoordinatesHandler<Graph> &crd)
      : g_(g),
        crd_(crd) {
  }

  std::vector<Path> FindGenomePaths(const std::vector<EdgeId> &from,
      const VertexId to, const size_t length_bound) {
    std::vector<Path> result;
    target_ = to;
    answer_.clear();

    for (const auto &e : from) {
      Path path_seq;
      const PosArray init_pos_array = crd_.GetEndPosArray(e);
      RunGenomeDFS(e, init_pos_array, length_bound, path_seq);
    }

    return answer_;
  }

  inline std::vector<Path> FindGenomePaths(const VertexId from,
      const VertexId to, const size_t length_bound) {
    return FindGenomePaths(g_.OutgoingEdges(from), to, length_bound);
  }

 private:
  void RunGenomeDFS(const EdgeId edge, const PosArray &cur_pos,
      const long long remaining, Path &path_seq) {
    path_seq.push_back(edge);

    VertexId vertex = g_.EdgeEnd(edge);
    //INFO("dfs: vertex " << g_.str(vertex) << " from edge " << g_.str(edge) << " thread " << int(cur_pos[0].first) << " - " << cur_pos[0].second);
    if (vertex == target_) {
      answer_.push_back(path_seq);
      path_seq.pop_back();
      return;
    }
    if (remaining < 0) {
      path_seq.pop_back();
      return;
    }


    for (EdgeId e : g_.OutgoingEdges(vertex)) {
      const PosArray further_array = crd_.FilterPosArray(cur_pos, e);
      if (further_array.size() == 0)
        continue;

      RunGenomeDFS(e, further_array, remaining - g_.length(e), path_seq);
    }
    path_seq.pop_back();
  }

  const Graph &g_;
  const CoordinatesHandler<Graph> &crd_;
  VertexId target_;
  std::vector<Path> answer_;
};

/*
template<class Graph>
class ReliableTargetedBoundedDijkstra : public Dijkstra<Graph> {
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef Dijkstra<Graph> base;

 public:
  ReliableTargetedBoundedDijkstra(const Graph &g, const VertexId target,
      const size_t bound, const size_t max_vertex_number)
      : base(g),
        target_(target),
        bound_(bound),
        max_vertex_number_(max_vertex_number),
        vertices_number_(0),
        vertex_limit_exceeded_(false) {
  }

  virtual bool CheckProcessVertex(const VertexId vertex, const size_t distance) {
      ++vertices_number_;

      if (vertices_number_ > max_vertex_number_)
          vertex_limit_exceeded_ = true;

      if (vertex == target_) {
        this->set_finished(true);
      }

      return (vertices_number_ < max_vertex_number_) && (distance <= bound_);
  }

  bool VertexLimitExceeded() const {
      return vertex_limit_exceeded_;
  }

 private:
  const VertexId target_;
  const size_t bound_;
  const size_t max_vertex_number_;
  size_t vertices_number_;
  bool vertex_limit_exceeded_;
};
*/

}

