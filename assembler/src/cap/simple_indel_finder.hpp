//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <iostream>
#include "coloring.hpp"
#include "omni/graph_processing_algorithm.hpp"

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

  bool mask_indels_;
  EdgeRemover<Graph> edge_remover_;

  size_t colors_number_;
  size_t coloring_version_;

  unordered_map<VertexId, u64int> processor_coloring_;
  unordered_set<EdgeId> preserve_list_;

  vector<EdgeId> sure_delete_list_;
  vector<EdgeId> probable_delete_list_;

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
      if (mask_indels_) {
        char state = 0;
        for (auto it = path_seq.rbegin(); it != path_seq.rend(); ++it) {
          if (!found_merge_point_)
            preserve_list_.insert(*it);
          else {
            if (preserve_list_.count(*it)) {
              state = 0;
            } else {
              if (state == 0) {
                sure_delete_list_.push_back(*it);
              } else {
                probable_delete_list_.push_back(*it);
              }

              state = 1;
            }
          }
        }
      }
      LogThread(path_seq);
      found_merge_point_ = true;

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

  void DeleteEdges() const {
    TRACE("DELETE EDGES");
    DeletingMergeHandler sure_merge_handler(sure_delete_list_);
    ConditionedSmartSetIterator<Graph, EdgeId, DeletingMergeHandler> sure_smart_it(
        g_, sure_delete_list_.begin(), sure_delete_list_.end(),
        sure_merge_handler);

    DeletingMergeHandler pro_merge_handler(probable_delete_list_);
    ConditionedSmartSetIterator<Graph, EdgeId, DeletingMergeHandler> pro_smart_it(
        g_, probable_delete_list_.begin(), probable_delete_list_.end(),
        pro_merge_handler);

    for (; !sure_smart_it.IsEnd(); ++sure_smart_it) {
      edge_remover_.DeleteEdge(*sure_smart_it);
    }

    while (true) {
      bool found = false;
      for (; !pro_smart_it.IsEnd(); ++pro_smart_it) {
        TRACE("Considering edge " << g_.str(*pro_smart_it));
        VertexId edge_end = g_.EdgeEnd(*pro_smart_it);
        if (g_.OutgoingEdgeCount(edge_end) == 0) {
          found = true;
          edge_remover_.DeleteEdge(*pro_smart_it);
          break;
        }
      }

      if (!found)
        break;

      pro_smart_it.reset();
    }
  }

  void CheckForIndelEvent(EdgeId starting_edge) {
    const VertexId starting_vertex = g_.EdgeEnd(starting_edge);
    const vector<EdgeId> &outgoing_edges = g_.OutgoingEdges(starting_vertex);

    TRACE("New indel event");

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
    preserve_list_.clear();
    sure_delete_list_.clear();
    probable_delete_list_.clear();

    vector<EdgeId> edge_seq_vector;
    for (auto it = outgoing_edges.begin(); it != outgoing_edges.end(); ++it) {
      GatheringDfs(*it, (1 << outgoing_edges_number) - 1, 0, edge_seq_vector);
    }

    if (found_merge_point_) {
      TRACE("Resolved split of color bunch");
      LogFooter();

      if (mask_indels_) {
        TRACE("Removing edges... (masking indels)");
        DeleteEdges();
        TRACE("Done");
      }
    }

  }

 public:
  SimpleIndelFinder(gp_t &gp, ColorHandler<Graph> &coloring, ostream &output_stream, const bool mask_indels = false)
      : gp_(gp),
        g_(gp_.g),
        coloring_(coloring),
        output_stream_(output_stream),
        mask_indels_(mask_indels),
        edge_remover_(gp_.g),
        colors_number_(coloring.max_colors()),
        coloring_version_(0),
        processor_coloring_(),
        preserve_list_(),
        sure_delete_list_(),
        probable_delete_list_() {
  }

  void FindIndelEvents() {
    INFO("Searching for In-Del events");
    for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
      CheckForIndelEvent(*it);
    }
    INFO("Searching for In-Del events ended");
  }

 private:
  class DeletingMergeHandler {
   public:
    DeletingMergeHandler(std::vector<EdgeId> to_delete)
        : delete_set_(to_delete.begin(), to_delete.end()) {
    }

    bool operator()(const std::vector<EdgeId> &old_edges, EdgeId new_edge) {
      bool ret = false;
      for (auto it = old_edges.begin(); it != old_edges.end(); ++it) {
        ret |= bool(delete_set_.count(*it));
        delete_set_.erase(*it);
      }
      if (ret) {
        delete_set_.insert(new_edge);
      }
      return ret;
    }

   private:
    std::unordered_set<EdgeId> delete_set_;
  };

	DECL_LOGGER("SimpleIndelFinder")
	;
};

}
