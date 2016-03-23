//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <iostream>
#include "coloring.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"
#include "path_projector.hpp"
#include "graph_traversal_constraints.hpp"

namespace cap {

template <class gp_t>
class SimpleIndelFinder {
  typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
  typedef std::vector<EdgeId> Path;
  typedef uint64_t u64int;

  gp_t &gp_;
  Graph &g_;
  ColorHandler<Graph> &coloring_;
  CoordinatesHandler<Graph> &coordinates_handler_;

  GraphTraversalConstraints<Graph> &graph_traversal_constraints_;

  ostream &output_stream_;

  bool mask_indels_;
  PathProjector<Graph> path_projector_;

  size_t colors_number_;
  size_t coloring_version_;

  unordered_map<VertexId, u64int> processor_coloring_;
  unordered_map<VertexId, size_t> vertices_distances_;

  VertexId restricted_vertex_;

  bool found_merge_point_;
  VertexId best_vertex_;
  size_t best_distance_to_vertex_;
  bool need_reflow_;

  vector<Path> alternative_paths_;
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

  size_t RelaxVertexDist(const VertexId v, const size_t dist) {
    const auto it = vertices_distances_.find(v);
    if (it == vertices_distances_.end()) {
      return vertices_distances_[v] = dist;
    } else {
      if (it->second > dist) {
        it->second = dist;
      }
      return it->second;
    }
  }

  void ColoringDfs(const EdgeId edge,
                   const size_t color_mask, const size_t path_length,
                   const size_t depth) {
    const VertexId vertex = g_.EdgeEnd(edge);
    if (vertex == restricted_vertex_) {
      return;
    }
    OrVertexColor(vertex, color_mask);
    size_t cur_dist = RelaxVertexDist(vertex, path_length);

    if (__builtin_popcountll(GetVertexColor(vertex)) >= 2) {
      TRACE("checking vertex from edge " << g_.str(edge) << " of dist " << cur_dist);
      if (!found_merge_point_) {
        best_vertex_ = vertex;
        best_distance_to_vertex_ = cur_dist;

        found_merge_point_ = true;
        TRACE("found merge point from edge " << g_.str(edge) << " length " << cur_dist);
      } else if (cur_dist < best_distance_to_vertex_) {
        best_vertex_ = vertex;
        best_distance_to_vertex_ = cur_dist;
        TRACE("found merge point from edge " << g_.str(edge) << " length " << cur_dist);
      }
      return;
    }

    if (depth >= kDfsDepthThreshold) {
      return;
    }
    for (EdgeId e : g_.OutgoingEdges(vertex)) {
      graph_traversal_constraints_.PushEdge(e);
      if (graph_traversal_constraints_.PathIsCorrect())
        ColoringDfs(e, color_mask,
            path_length + g_.length(e), depth + 1);
      graph_traversal_constraints_.PopEdge();
    }
  }

  // returns if endpoint was found or not
  bool GatheringDfs(const EdgeId edge,
                    /*const size_t color_mask_needed,*/ const size_t depth,
                    vector<EdgeId> &path_seq) {
    VertexId vertex = g_.EdgeEnd(edge);
    if (vertex == restricted_vertex_) {
      return false;
    }

    path_seq.push_back(edge);

    //if (GetVertexColor(vertex) == color_mask_needed) {
    //if (__builtin_popcount(GetVertexColor(vertex)) >= 2) {
    if (vertex == best_vertex_) {
      /*
      INFO("found final vertex " << g_.str(vertex));
      if (coordinates_handler_.GetContiguousThreads(path_seq).size() !=
             pos_array.size()) {
        INFO("" << coordinates_handler_.GetContiguousThreads(path_seq).size() << " " <<
            pos_array.size());
        VERIFY(false);
      }
      */

      alternative_paths_.push_back(path_seq);
      path_seq.pop_back();
      return true;
    }

    if (depth >= kDfsDepthThreshold) {
      path_seq.pop_back();
      return false;
    }

    for (EdgeId e : g_.OutgoingEdges(vertex)) {
      graph_traversal_constraints_.PushEdge(e);
      if (graph_traversal_constraints_.PathIsCorrect())
        GatheringDfs(e,
            /*color_mask_needed,*/ depth + 1, path_seq);
      graph_traversal_constraints_.PopEdge();
    }
    path_seq.pop_back();
    return false;
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

  void CollapsePaths() {
    if (alternative_paths_.size() == 0)
      return;
    bool has_short_enough = false;
    for (const auto &path : alternative_paths_) {
      if (GetPathLength(path) < 10 * g_.k()) {
          //todo move threshold out of here
        //TRACE("Too long path: " << GetPathLength(path));
        has_short_enough = true;
        break;
      }
    }
    if (!has_short_enough) {
      alternative_paths_.clear();
      return;
    }

    bool success = path_projector_.CollapsePaths(alternative_paths_);
    if (!success) {
      TRACE("Could not collapse paths: " << alternative_paths_ << " ("
          << alternative_paths_.size() << " paths)");
      for (const auto &path : alternative_paths_) {
        PaintPath(path);
      }
    }
    alternative_paths_.clear();
  }

  size_t EstimateSNPNumber(const size_t branch_len) {
    if (branch_len < 2) return 0;
    return 1 + (branch_len - (g_.k() + 1) + g_.k() - 1) / g_.k();
  }

  void AnalyseThreadLengths() {
    std::vector<size_t> thread_lengths;
    for (const auto &path : alternative_paths_) {
      thread_lengths.push_back(GetPathLength(path));
    }
    bool is_simple_snp = true;
    bool is_simple_indel = false;
    bool equal_lengths = true;

    size_t prev = 0;
    size_t min_branch = size_t(-1);
    for (size_t len : thread_lengths) {
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
  }

  void CheckForIndelEvent(const VertexId starting_vertex) {
    TRACE("New indel event");

    // Check that is is interesting
    size_t outgoing_edges_number = g_.OutgoingEdgeCount(starting_vertex);
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
    found_merge_point_ = false;
    restricted_vertex_ = starting_vertex;
    size_t branch_num = 0;
    for (EdgeId e : g_.OutgoingEdges(starting_vertex)) {
      graph_traversal_constraints_.PushEdge(e);
      ColoringDfs(e, 1 << branch_num, g_.length(e), 0);
      branch_num++;
      graph_traversal_constraints_.PopEdge();
    }

    vertices_distances_.clear();

    if (found_merge_point_) {
      vector<EdgeId> edge_seq_vector;
      for (EdgeId e : g_.OutgoingEdges(starting_vertex)) {
        graph_traversal_constraints_.PushEdge(e);
        GatheringDfs(e,
            /*(1 << outgoing_edges_number) - 1,*/ 0, edge_seq_vector);
        graph_traversal_constraints_.PopEdge();
      }

      AnalyseThreadLengths();

      TRACE("Resolved split of color bunch");

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
      GraphTraversalConstraints<Graph> &graph_traversal_constraints,
      ostream &output_stream, const bool mask_indels = false)
      : gp_(gp),
        g_(gp_.g),
        coloring_(coloring),
        coordinates_handler_(coordinates_handler),
        graph_traversal_constraints_(graph_traversal_constraints),
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
        CheckForIndelEvent(*it);
      } while (need_reflow_);
    }
    INFO("Found around " << snps_/2 << "+" << unknown_snp_/2 << "=" <<
        (snps_+unknown_snp_)/2 << " SNPs and " << indels_/2 << "+" <<
        unknown_indel_/2 << " indels");
  }

 private:
    DECL_LOGGER("SimpleIndelFinder")
    ;
};

}
