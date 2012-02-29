// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_MUXSEARCHRESULT_H
#define PATHS_MUXSEARCHRESULT_H

#include "paths/MuxWalkGraph.h"
#include "paths/KmerPath.h"

/// The structure returned by KmerPathMuxSearcher::FindClosures

struct MuxSearchResult {

  MuxSearchResult() : walk_graph_ptr(NULL) { clear(); }
  ~MuxSearchResult() { delete walk_graph_ptr; }

  MuxWalkGraph* walk_graph_ptr;
  int num_closures_found;
  longlong num_states_explored;
  bool hit_search_limit;

  vecKmerPath all_closures;  // This is only filled in by explicit request
  bool hit_closure_limit;

  void clear() {
    // MuxWalkGraph::clear() leaves its hash_maps with zero
    // elements but possibly a large bucket count, trashing
    // memory locality.  Instead, create a new MuxWalkGraph.
    delete walk_graph_ptr;
    walk_graph_ptr = new MuxWalkGraph;

    num_closures_found = 0;
    num_states_explored = 0;
    hit_search_limit = false;
    hit_closure_limit = false;
    all_closures.clear();
  }

  MuxWalkGraph& WalkGraph() { return *walk_graph_ptr; }
  const MuxWalkGraph& WalkGraph() const { return *walk_graph_ptr; }


  // Monte Carlo estimation of the number of cloures (in Mux space)
  int EstimateNumClosures() {
    return WalkGraph().EstimateNumClosures();
  }
  
  void CalculateAllClosures( const vecKmerPath* pathsFw,
			     const vecKmerPath* pathsRc,
			     unsigned int max_ans=0 ) {
    hit_closure_limit = 
      ! WalkGraph().AllPaths( pathsFw, pathsRc, all_closures, max_ans );
  }
};

#endif
