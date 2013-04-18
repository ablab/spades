//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <iostream>
#include "coloring.hpp"
#include "omni/graph_processing_algorithm.hpp"
#include "path_projector.hpp"

namespace cap {

template <class gp_t>
class SimpleIndelFinder {
  typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
  typedef std::vector<EdgeId> Path;
  typedef std::vector<std::pair<unsigned char, size_t> > PosArray;
  typedef uint64_t u64int;

  gp_t &gp_;
  Graph &g_;
  ColorHandler<Graph> &coloring_;
  CoordinatesHandler<Graph> &coordinates_handler_;

  ostream &output_stream_;

  bool mask_indels_;
  PathProjector<Graph> path_projector_;

  size_t colors_number_;
  size_t coloring_version_;

  unordered_map<VertexId, u64int> processor_coloring_;

  VertexId restricted_vertex_;

  bool found_merge_point_;
  bool need_reflow_;

  vector<Path> alternative_paths_;
  vector<size_t> thread_lengths_;
  size_t snps_;
  size_t indels_;
  size_t unknown_snp_;
  size_t unknown_indel_;

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
    thread_lengths_.push_back(thread_length);

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
    u64int &current_val = processor_coloring_[vertex];
    if ((current_val & kProcessorColorMask) < coloring_version_) {
      // yes, it clears color bits
      current_val = coloring_version_;
    }
    current_val |= u64int(color_set) << kProcessorColorShift;

    //processor_coloring_[vertex] = current_val;
  }

  inline size_t GetVertexColor(VertexId vertex) {
    const u64int current_val = processor_coloring_[vertex];
    if ((current_val & kProcessorColorMask) < coloring_version_) {
      return 0;
    }
    return size_t(current_val >> kProcessorColorShift);
  }


  inline TColorSet GetIncomingColoring(const VertexId vertex) const {
    TColorSet incoming_coloring;
    vector<EdgeId> incoming_edges = g_.IncomingEdges(vertex);
    for (auto it = incoming_edges.begin(); it != incoming_edges.end(); ++it) {
      incoming_coloring |= coloring_.Color(*it);
    }
    return incoming_coloring;
  }


  inline bool CheckColorSetExistence(const vector<EdgeId> &edges,
      const TColorSet &coloring) const {
    for (auto it = edges.begin(); it != edges.end(); ++it) {
      if (coloring_.Color(*it) == coloring) {
        return true;
      }
    }
    return false;
  }


  void ColoringDfs(const EdgeId edge, const PosArray &pos_array,
                   const size_t color_mask, const size_t depth) {
    const VertexId vertex = g_.EdgeEnd(edge);
    if (vertex == restricted_vertex_) {
      return;
    }
    OrVertexColor(vertex, color_mask);

    if (depth >= kDfsDepthThreshold) {
      return;
    }
    vector<EdgeId> further_edges = g_.OutgoingEdges(vertex);
    for (auto it = further_edges.begin(); it != further_edges.end(); ++it) {
      const PosArray further_array = coordinates_handler_.FilterPosArray(
          pos_array, *it);
      if (further_array.size() == 0)
        continue;
      ColoringDfs(*it, further_array, color_mask, depth + 1);
    }
  }

  // returns if endpoint was found or not
  bool GatheringDfs(const EdgeId edge, const PosArray &pos_array,
                    const size_t color_mask_needed, const size_t depth,
                    vector<EdgeId> &path_seq) {
    VertexId vertex = g_.EdgeEnd(edge);
    if (vertex == restricted_vertex_) {
      return false;
    }

    path_seq.push_back(edge);

    if (GetVertexColor(vertex) == color_mask_needed) {
      found_merge_point_ = true;

      INFO("found final vertex " << g_.str(vertex));
      alternative_paths_.push_back(path_seq);
      if (coordinates_handler_.GetContiguousThreads(path_seq).size() !=
             pos_array.size()) {
        INFO("" << coordinates_handler_.GetContiguousThreads(path_seq).size() << " " <<
            pos_array.size());
        VERIFY(false);
      }
      // dirty hack now
      //coloring_.PaintEdge(edge, TColorSet::SingleColor(colors_number_));

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
      const PosArray further_array = coordinates_handler_.FilterPosArray(
          pos_array, *it);
      if (further_array.size() == 0)
        continue;

      bool is_path_to_merge_point = GatheringDfs(*it, further_array,
          color_mask_needed, depth + 1, path_seq);
      found_merge_point |= is_path_to_merge_point;
    }
    path_seq.pop_back();
    return found_merge_point;
  }

  void PaintPath(const Path &path) {
    for (const auto e : path)
      coloring_.PaintEdge(e, TColorSet::SingleColor(colors_number_));
  }

  size_t GetPathLength(const Path &p) const {
    size_t res = 0;
    for (const auto &edge : p) {
      res += g_.length(edge);
    }
    return res;
  }

  std::vector<Path> GetShortestConsistentPaths(const std::vector<Path> &paths) {
    std::vector<Path> good_paths;
    for (size_t i = 0; i < paths.size(); ++i) {
      bool found = false;
      for (size_t j = 0; j < paths.size(); ++j) {
        if (i == j) continue;
        if (g_.EdgeEnd(paths[i].back()) == g_.EdgeEnd(paths[j].back())) {
          found = true;
          break;
        }
      }
      if (found) {
        good_paths.push_back(paths[i]);
      }
    }
    size_t min_len = size_t(-1);
    VertexId min_len_v = VertexId(0);
    for (size_t i = 0; i < good_paths.size(); ++i) {
      size_t cur_len = GetPathLength(good_paths[i]);
      if (cur_len < min_len) {
        min_len = cur_len;
        min_len_v = g_.EdgeEnd(good_paths[i].back());
      }
    }
    for (size_t i = 0; i < good_paths.size(); ++i) {
      if (g_.EdgeEnd(good_paths[i].back()) != min_len_v) {
        std::swap(good_paths[i], good_paths.back());
        good_paths.pop_back();
        i--;
      }
    }
    if (good_paths.size() != paths.size()) {
      if (good_paths.size() == 0) {
        for (const auto &path : paths)
          PaintPath(path);
      } else {
        need_reflow_ = true;
      }
    }
    return good_paths;
  }

  void CollapsePaths() {
    alternative_paths_ = GetShortestConsistentPaths(alternative_paths_);
    INFO("after shorting " << alternative_paths_.size() << " paths");
    
    if (alternative_paths_.size() == 0)
      return;
    for (const auto &path : alternative_paths_) {
      if (GetPathLength(path) > 10 * g_.k()) {
        INFO("Too long path: " << GetPathLength(path));
        alternative_paths_.clear();
        need_reflow_ = false;
        return;
      }
    }

    bool success = path_projector_.CollapsePaths(alternative_paths_);
    if (!success) {
      INFO("Could not collapse paths: " << alternative_paths_ << " ("
          << alternative_paths_.size() << " paths)");
      for (const auto &path : alternative_paths_) {
        PaintPath(path);
      }
    }
    alternative_paths_.clear();
  }

  // Deprecated
  /*
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
  */

  size_t EstimateSNPNumber(const size_t branch_len) {
    if (branch_len < 2) return 0;
    return 1 + (branch_len - (g_.k() + 1) + g_.k() - 1) / g_.k();
  }

  void AnalyseThreadLengths() {
    bool is_simple_snp = true;
    bool is_simple_indel = false; 
    bool equal_lengths = true;

    size_t prev = 0;
    size_t min_branch = size_t(-1);
    for (size_t len : thread_lengths_) {
      if (len != g_.k() + 1) {
        is_simple_snp = false;
      }
      if (prev && len != prev) {
        equal_lengths = false;
      }
      if (len == g_.k()) {
        is_simple_indel = true;
      }

      prev = len;
      min_branch = min(min_branch, len);
    }

    if (is_simple_snp)
      snps_++;
    else if (is_simple_indel)
      indels_++;
    else 
    if (equal_lengths) {
      snps_++;
      unknown_snp_ += EstimateSNPNumber(prev) - 1;
    } else {
      indels_++;
      unknown_snp_ += EstimateSNPNumber(min_branch) - 1;
    }

    thread_lengths_.clear();
  }

  void CheckForIndelEvent(const VertexId starting_vertex) {
    const vector<EdgeId> &outgoing_edges = g_.OutgoingEdges(starting_vertex);

    TRACE("New indel event");

    // Check that is is interesting
    size_t outgoing_edges_number = outgoing_edges.size();
    if (outgoing_edges_number <= 1 || outgoing_edges_number > kOutgingEdgesNumberThreshold) {
      // Nothing to do here
      return;
    }

    /*
    // Generate initial colorset
    TColorSet initial_coloring = GetIncomingColoring(starting_vertex);

    // Check that there exists split of colors
    if (CheckColorSetExistence(outgoing_edges, initial_coloring)) {
      // Nothing to do here
      return;
    }
    */
 
    // Dfs and try to resolve all splits
    ++coloring_version_;
    restricted_vertex_ = starting_vertex;
    size_t branch_num = 0;
    for (auto it = outgoing_edges.begin(); it != outgoing_edges.end(); ++it) {
      const PosArray init_pos_array = coordinates_handler_.GetEndPosArray(*it);
      ColoringDfs(*it, init_pos_array, 1 << branch_num, 0);
      branch_num++;
    }

    found_merge_point_ = false;

    vector<EdgeId> edge_seq_vector;
    for (auto it = outgoing_edges.begin(); it != outgoing_edges.end(); ++it) {
      const PosArray init_pos_array = coordinates_handler_.GetEndPosArray(*it);
      GatheringDfs(*it, init_pos_array,
          (1 << outgoing_edges_number) - 1, 0, edge_seq_vector);
    }

    if (found_merge_point_) {
      //AnalyseThreadLengths();

      TRACE("Resolved split of color bunch");
      LogFooter();

      if (mask_indels_) {
        TRACE("Removing edges... (masking indels)");
        CollapsePaths();
        //DeleteEdges();
        TRACE("Done");
      }
    }

  }

 public:
  SimpleIndelFinder(gp_t &gp, ColorHandler<Graph> &coloring,
      CoordinatesHandler<Graph> &coordinates_handler,
      ostream &output_stream, const bool mask_indels = false)
      : gp_(gp),
        g_(gp_.g),
        coloring_(coloring),
        coordinates_handler_(coordinates_handler),
        output_stream_(output_stream),
        mask_indels_(mask_indels),
        path_projector_(g_, coordinates_handler_),
        colors_number_(coloring.max_colors()),
        coloring_version_(0),
        processor_coloring_() {
  }

  void FindIndelEvents() {
    indels_ = 0;
    snps_ = 0;
    unknown_snp_ = 0;
    unknown_indel_ = 0;
    INFO("Searching for In-Del events");
    for (auto it = g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
      do {
        need_reflow_ = false;
        TRACE("in");
        CheckForIndelEvent(*it);
        TRACE("out");
      } while (need_reflow_);
    }
    INFO("Found around " << snps_/2 << "+" << unknown_snp_/2 << "=" <<
        (snps_+unknown_snp_)/2 << " SNPs and " << indels_/2 << "+" <<
        unknown_indel_/2 << " indels");
    INFO("Searching for In-Del events ended");
  }

 private:
	DECL_LOGGER("SimpleIndelFinder")
	;
};

}
