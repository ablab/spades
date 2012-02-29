// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_MUXSEARCHSTATE
#define PATHS_MUXSEARCHSTATE

#include "paths/Mux.h"

/// The various stacks that are used in a Mux search.
/// They are bundled together in an outside struct to 
/// make it easy to pass them to outside helpers.

struct MuxSearchState {
  vec<Mux> path_to_here;       // empty at start of search
  vec< vec<Mux> > to_explore;  // initial value set by SetupOpeners()
  vec<int> dist_stack;         // empty at start of search.
  
  void clear() {
    path_to_here.clear(); to_explore.clear(); dist_stack.clear();
  }
};



#endif
