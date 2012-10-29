//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "omni/omni_utils.hpp"
#include "sequence/sequence_tools.hpp"
#include "runtime_k.hpp"
#include "edge_index.hpp"
#include <cstdlib>

namespace debruijn_graph {
template<class Graph>
class MismatchMasker {

  typedef typename Graph::EdgeId EdgeId;

  const Graph& g_;
  public:
  struct MismatchInfo {
	  size_t position;
// Ratio:
	  double ratio;
	  MismatchInfo(size_t position, double ratio): position(position), ratio(ratio){}
  };
  map<EdgeId, vector<MismatchInfo> > mismatch_map;

  MismatchMasker(const Graph& g) : g_(g){
  }
  void insert(EdgeId edge, size_t position, double ratio){
	  if (mismatch_map.find(edge) == mismatch_map.end()) {

		  vector<MismatchInfo> tmp;
		  mismatch_map.insert(make_pair(edge, tmp));
		  vector<MismatchInfo> rc_tmp;
		  mismatch_map.insert(make_pair(g_.conjugate(edge), rc_tmp));
	  }
	  mismatch_map[edge].push_back(MismatchInfo(position, ratio));
	  mismatch_map[g_.conjugate(edge)].push_back(MismatchInfo(g_.length(edge) + g_.k() - position - 1, ratio));

  }
  string MaskedEdgeNucls(EdgeId edge, double cutoff) {
	  Sequence s_edge = g_.EdgeNucls(edge);
	  string s = s_edge.str();
	  if (mismatch_map[edge].size() > 0) {
		  DEBUG("in edge length " << g_.length(edge)<< " replaced " << mismatch_map[edge].size() << "mismatches");
	  }
	  for(size_t i = 0; i < mismatch_map[edge].size(); i++)
		  if (mismatch_map[edge][i].ratio > cutoff && is_nucl(s[mismatch_map[edge][i].position]))
			  s[mismatch_map[edge][i].position] =  char (s[mismatch_map[edge][i].position] + 'a' - 'A') ;
	  return s;
  }
};

}


