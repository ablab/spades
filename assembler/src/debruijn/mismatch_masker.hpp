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
	  vector<size_t> counts;
	  MismatchInfo(size_t position, double ratio, vector<size_t> count): position(position), ratio(ratio), counts(count.begin(), count.end()){}
  };
  map<EdgeId, vector<MismatchInfo> > mismatch_map;

  MismatchMasker(const Graph& g) : g_(g){
  }

  void insert(EdgeId edge, size_t position, double ratio, vector<size_t> counts){
	  if (mismatch_map.find(edge) == mismatch_map.end()) {

		  vector<MismatchInfo> tmp;
		  mismatch_map.insert(make_pair(edge, tmp));
		  vector<MismatchInfo> rc_tmp;
		  mismatch_map.insert(make_pair(g_.conjugate(edge), rc_tmp));
	  }
	  VERIFY(counts.size() == 4);
	  mismatch_map[edge].push_back(MismatchInfo(position, ratio, counts));
	  vector<size_t> reversed_counts(4, 0);
	  for(size_t i = 0; i < 4; i ++)
		  reversed_counts[complement(i)] = counts[i];
	  VERIFY(reversed_counts.size() == 4);
	  mismatch_map[g_.conjugate(edge)].push_back(MismatchInfo(g_.length(edge) + g_.k() - position - 1, ratio, reversed_counts));

  }
  string MaskedEdgeNucls(EdgeId edge, double cutoff) {
	  Sequence s_edge = g_.EdgeNucls(edge);
	  string s = s_edge.str();
	  if (mismatch_map[edge].size() > 0) {
		  DEBUG("in edge length " << g_.length(edge)<< " replacing " << mismatch_map[edge].size() << "mismatches");
	  }
	  for(size_t i = 0; i < mismatch_map[edge].size(); i++) {
		  DEBUG(mismatch_map[edge][i].position  << " " << mismatch_map[edge][i].ratio <<" " <<  mismatch_map[edge][i].counts.size());

		  if (mismatch_map[edge][i].ratio > cutoff && is_nucl(s[mismatch_map[edge][i].position])) {
			  DEBUG(s[mismatch_map[edge][i].position]  << "  before " );
			  VERIFY(mismatch_map[edge][i].position < s.length() && mismatch_map[edge][i].position >= 0);
			  if (mismatch_map[edge][i].ratio > 1) {
				  DEBUG("replacing...");
				  for(size_t ii = 0; ii < mismatch_map[edge][i].counts.size(); ii++) {
					  DEBUG(" counts "<< ii << " " << mismatch_map[edge][i].counts[ii])
			  	  }
				  size_t max_count = 0;
	  	  	  	  size_t max_i = 0;
	  	  	  	  VERIFY(mismatch_map[edge][i].counts.size() == 4);
				  for(size_t ii = 0; ii < mismatch_map[edge][i].counts.size(); ii++)
					  if (max_count < mismatch_map[edge][i].counts[ii]) {
						  max_count = mismatch_map[edge][i].counts[ii];
						  max_i = ii;
					  }
				  DEBUG("max_i found: "<< max_i);
				  s[mismatch_map[edge][i].position] = char(nucl(max_i)) ;
				  DEBUG("replaced");
			  }
			  //s[mismatch_map[edge][i].position] = char(s[mismatch_map[edge][i].position] +'a' - 'A');
			  //if (mismatch_map[edge][i].ratio > 0.5)
			  s[mismatch_map[edge][i].position] = char(s[mismatch_map[edge][i].position] +'a' - 'A');
			  if ((mismatch_map[edge][i].position >= 1 && ! is_nucl(s[mismatch_map[edge][i].position - 1])) || (mismatch_map[edge][i].position <s.length() -1  && ! is_nucl(s[mismatch_map[edge][i].position + 1]))){
//				 INFO("replacement to 'N' blocked") ;
			  }else {
//				  if (mismatch_map[edge][i].position >= 1 && s[mismatch_map[edge][i].position - 1] == 'N')
//					  INFO("2 succesive N'th. BUT HOW?");
				  s[mismatch_map[edge][i].position] = char('N');

			  }
			  DEBUG(s[mismatch_map[edge][i].position]  << "  " <<'a' - 'A');
		  }
	  }
	  return s;
  }
};

}


