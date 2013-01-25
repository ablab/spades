//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <iostream>
#include "coloring.hpp"

namespace cap {

template <class gp_t>
class SimpleIndelFinder {
  typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
  typedef uint64_t u64int;

  gp_t &gp_;
  Graph &g_;
  ColorHandler<Graph> &coloring_;
  ostream &output_stream_;

  size_t colors_number_;
  size_t coloring_version_;

  map<VertexId, u64int> processor_coloring_;
  VertexId restricted_vertex_;

  bool found_merge_point_;

  // Maximum number of outgoing edges we are interested in
  // If more, we do not consider this case at all
  const static size_t kOutgingEdgesNumberThreshold = 4;
  const static size_t kDfsDepthThreshold = 4;
  const static size_t kProcessorColorShift = 32;
  const static size_t kProcessorColorMask = (1ll << kProcessorColorShift) - 1;

  inline void LogHeader() {
    output_stream_ << "Indel event:" << std::endl;
  }

  inline void LogFooter() {
    output_stream_ << "End indel event" << std::endl << std::endl;
  }

  inline void LogEdge(EdgeId edge) {
    output_stream_ << gp_.edge_pos.str(edge);
  }

  inline void LogThread(const vector<EdgeId> &path_seq) {
    if (!found_merge_point_) {
      LogHeader();
    }
    found_merge_point_ = true;

    // Recover some thread params
    TColorSet thread_colorset = TColorSet::AllColorsSet(colors_number_);
    size_t thread_length = 0;
    for (auto it = path_seq.begin(); it != path_seq.end(); ++it) {
      thread_colorset &= coloring_.Color(*it);
      thread_length += g_.length(*it);
    }

    output_stream_ << "Thread: length " << thread_length << " genomes (";

    bool is_first = true;
    for (size_t i = 0; i < colors_number_; ++i) {
      if (thread_colorset[i] == 0) {
        continue;
      }
      if (!is_first) {
        output_stream_ << ',';
      } else {
        is_first = false;
      }
      output_stream_ << i;
    }

    output_stream_ << ") path (" << std::endl;

    is_first = true;
    for (auto it = path_seq.begin(); it != path_seq.end(); ++it) {
      LogEdge(*it);
      output_stream_ << std::endl;
    }

    output_stream_ << ')' << std::endl;
  }

  inline void OrVertexColor(VertexId vertex, size_t color_set) {
    u64int current_val = processor_coloring_[vertex];
    if ((current_val & kProcessorColorMask) < coloring_version_) {
      // yes, it clears color bits
      current_val = coloring_version_;
    }
    current_val |= u64int(color_set) << kProcessorColorShift;

    processor_coloring_[vertex] = current_val;
  }

  inline size_t GetVertexColor(VertexId vertex) {
    return size_t(processor_coloring_[vertex] >> kProcessorColorShift);
  }


  inline TColorSet GetIncomingColoring(const VertexId vertex) {
    TColorSet incoming_coloring;
    vector<EdgeId> incoming_edges = g_.IncomingEdges(vertex);
    for (auto it = incoming_edges.begin(); it != incoming_edges.end(); ++it) {
      incoming_coloring |= coloring_.Color(*it);
    }
    return incoming_coloring;
  }


  inline bool CheckColorSetExistence(const vector<EdgeId> &edges, const TColorSet &coloring) {
    for (auto it = edges.begin(); it != edges.end(); ++it) {
      if (coloring_.Color(*it) == coloring) {
        return true;
      }
    }
    return false;
  }


  void ColoringDfs(EdgeId edge, size_t color_mask, size_t depth) {
    VertexId vertex = g_.EdgeEnd(edge);
    if (vertex == restricted_vertex_) {
      return;
    }
    OrVertexColor(vertex, color_mask);

    if (depth >= kDfsDepthThreshold) {
      return;
    }
    vector<EdgeId> further_edges = g_.OutgoingEdges(vertex);
    for (auto it = further_edges.begin(); it != further_edges.end(); ++it) {
      ColoringDfs(*it, color_mask, depth + 1);
    }
  }

  // returns if endpoint was found or not
  bool GatheringDfs(EdgeId edge, size_t color_mask_needed, size_t depth, vector<EdgeId> &path_seq) {
    VertexId vertex = g_.EdgeEnd(edge);
    if (vertex == restricted_vertex_) {
      return false;
    }

    path_seq.push_back(edge);

    if (GetVertexColor(vertex) == color_mask_needed) {
      LogThread(path_seq);

      // dirty hack now
      coloring_.PaintEdge(edge, TColorSet::SingleColor(colors_number_));

      path_seq.pop_back();
      return true;
    }

    if (depth >= kDfsDepthThreshold) {
      path_seq.pop_back();
      return false;
    }

    bool found_merge_point = false;
    vector<EdgeId> further_edges = g_.OutgoingEdges(vertex);
    for (auto it = further_edges.begin(); it != further_edges.end(); ++it) {

      bool is_path_to_merge_point = GatheringDfs(*it, color_mask_needed, depth + 1, path_seq);
      found_merge_point |= is_path_to_merge_point;

      if (is_path_to_merge_point) {
        // dirty hack now
        coloring_.PaintEdge(edge, TColorSet::SingleColor(colors_number_));
      }
    }
    path_seq.pop_back();
    return found_merge_point;
  }

  void CheckForIndelEvent(EdgeId starting_edge) {
    const VertexId starting_vertex = g_.EdgeEnd(starting_edge);
    const vector<EdgeId> outgoing_edges = g_.OutgoingEdges(starting_vertex);

    // Check that is is interesting
    size_t outgoing_edges_number = outgoing_edges.size();
    if (outgoing_edges_number <= 1 || outgoing_edges_number > kOutgingEdgesNumberThreshold) {
      // Nothing to do here
      return;
    }

    // Generate initial colorset
    TColorSet initial_coloring = GetIncomingColoring(starting_vertex);

    // Check that there exists split of colors
    if (CheckColorSetExistence(outgoing_edges, initial_coloring)) {
      // Nothing to do here
      return;
    }
 
    // Dfs and try to resolve all splits
    ++coloring_version_;
    restricted_vertex_ = starting_vertex;
    size_t branch_num = 0;
    for (auto it = outgoing_edges.begin(); it != outgoing_edges.end(); ++it) {
      ColoringDfs(*it, 1 << branch_num, 0);
      branch_num++;
    }

    found_merge_point_ = false;
    vector<EdgeId> edge_seq_vector;
    for (auto it = outgoing_edges.begin(); it != outgoing_edges.end(); ++it) {
      GatheringDfs(*it, (1 << outgoing_edges_number) - 1, 0, edge_seq_vector);
    }

    if (found_merge_point_) {
      TRACE("Resolved split of color bunch");
      LogFooter();
    }

  }

 public:
  SimpleIndelFinder(gp_t &gp, ColorHandler<Graph> &coloring, ostream &output_stream)
      : gp_(gp),
        g_(gp_.g),
        coloring_(coloring),
        output_stream_(output_stream),
        colors_number_(coloring.max_colors()),
        coloring_version_(0),
        processor_coloring_() {
  }

  void FindIndelEvents() {
    INFO("Searching for In-Del events");
    for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
      CheckForIndelEvent(*it);
    }
    INFO("Searching for In-Del events ended");
  }
};

}
