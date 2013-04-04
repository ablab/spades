//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <iostream>
#include "coloring.hpp"
#include "compare_standard.hpp"
#include "comparison_utils.hpp"
#include "omni/path_processor.hpp"

namespace cap {

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
      /*ostream &output_stream,*/ const std::string base_pic_file_name,
      const bool mask_inversed = false)
      : gp_(gp),
        g_(gp_.g),
        coloring_(coloring),
        //output_stream_(output_stream),
        base_pic_file_name_(base_pic_file_name),
        num_cycles_found_(0),
        mask_inversed_(false) {
  }

  void FindInversionEvents() {
    num_cycles_found_ = 0;

    INFO("Searching for inversions");
    for (auto it = g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
      CheckForCycle(*it);
    }
    INFO("Searching for inversions done. Cycles found: " << num_cycles_found_);
  }

 private:
  typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

  typedef std::vector<EdgeId> EdgeList;
  typedef std::vector<VertexId> VertexList;
  typedef uint64_t u64int;


  gp_t &gp_;
  Graph &g_;
  ColorHandler<Graph> &coloring_;

  //ostream &output_stream_;
  std::string base_pic_file_name_;
  size_t num_cycles_found_;

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

          ProcessFoundCycle(VertexList({v1, v2, v3_list[0], v4}), *it4);

          // should return here
        }
      }
    }
  }

  void ProcessFoundCycle(const VertexList &cycle, const EdgeId edge) {
    stringstream v_list_str;
    for (auto it = cycle.begin(); it != cycle.end(); ++it) {
      v_list_str << g_.str(*it) << " ";
    }
    DEBUG("cycle found: " << v_list_str.str());

    const std::string edge_pic_name = base_pic_file_name_ + "_" +
        ToString(num_cycles_found_) + ".dot";
    const std::string path_pic_name = base_pic_file_name_ + "_path_" +
        ToString(num_cycles_found_) + ".dot";

    PrintColoredGraphAroundEdge(g_, coloring_, edge, gp_.edge_pos,
        edge_pic_name);
    
    // Awful!...
    FindAndPrintPath(cycle[2], g_.conjugate(cycle[3]), path_pic_name) ||
    FindAndPrintPath(g_.conjugate(cycle[2]), cycle[3], path_pic_name) ||
    FindAndPrintPath(cycle[0], g_.conjugate(cycle[1]), path_pic_name) ||
    FindAndPrintPath(g_.conjugate(cycle[0]), cycle[1], path_pic_name);

    num_cycles_found_++;
  }

  std::vector<EdgeList> GetSoloColoredEdgeLists(const EdgeList &all_edges) const {
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
    for (size_t i = 0; i < coloring_.max_colors(); ++i) {
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

  inline bool FindAndPrintPath(const VertexId v1, const VertexId v2, const std::string &out_file) const {
    TRACE("Finding paths from " << g_.str(v1) << " to " << g_.str(v2));
    const static size_t max_length = 20000;

    PathStorageCallback<Graph> callback(g_);
    PathProcessor<Graph> pp(g_, 0, max_length, v1, v2, callback);
    pp.Process();
    const vector<EdgeList>& paths = callback.paths();

    if (paths.size() == 0) {
      INFO("Something's wrong!");
      return false;
    }
    //PrintPath(paths[0], out_file);
    TRACE("Finding path done; length=" << PathLength(paths[0])
        << "; written in " << out_file);
    return true;
  }

  inline void PrintPath(const EdgeList &path, const std::string out_file) const {
    const static size_t edge_length = 1;
    const static size_t max_vertices = 100;
    MappingPath<EdgeId> mpath = TrivialMappingPath(g_, path);
    //Path<EdgeId> cpath(path, mpath.start_pos(), mpath.end_pos());

    LengthIdGraphLabeler<Graph> basic_labeler(g_);
    EdgePosGraphLabeler<Graph> pos_labeler(g_, gp_.edge_pos);
    CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);

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

}

