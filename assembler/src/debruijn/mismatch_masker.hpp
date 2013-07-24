//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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
class MismatchMasker: public GraphActionHandler<Graph> {

  typedef typename Graph::EdgeId EdgeId;

  const Graph& g_;
  public:
  struct MismatchInfo {
	  size_t position;
// Ratio:
	  double ratio;
	  vector<size_t> counts;
	  double cutoff;
	  MismatchInfo(size_t position, double ratio, vector<size_t> count, double cutoff = 0): position(position), ratio(ratio), counts(count.begin(), count.end()), cutoff(cutoff){}
  };
  map<EdgeId, vector<MismatchInfo> > mismatch_map;

  MismatchMasker(Graph& g) : GraphActionHandler<Graph> (g, "Mismatch Masker"), g_(g) {
  }

  virtual void HandleDelete(EdgeId edge) {
//  		SetCoverage(edge, 0);
	  mismatch_map.erase(edge);
  }

  virtual void HandleMerge(const vector<EdgeId>& oldEdges, EdgeId newEdge) {
	  int len_shift = 0;
	  for (auto it = oldEdges.begin(); it != oldEdges.end(); ++it) {
		  if (mismatch_map.find(*it) != mismatch_map.end())
			  for (auto mis_t = mismatch_map[*it].begin(); mis_t != mismatch_map[*it].end(); ++mis_t) {
				  if(mismatch_map.find(newEdge) == mismatch_map.end()) {
					  vector<MismatchInfo> tmp_v;
					  mismatch_map.insert(make_pair(newEdge, tmp_v));
				  }
				  mismatch_map[newEdge].push_back(MismatchInfo(mis_t->position + len_shift, mis_t->ratio, mis_t->counts, mis_t->cutoff));
			  }
		  len_shift += (int) g_.length(*it);
  	  }
  }


  void insert(EdgeId edge, size_t position, double ratio, vector<size_t> counts, double cutoff = 0){
	  if (mismatch_map.find(edge) == mismatch_map.end()) {

		  vector<MismatchInfo> tmp;
		  mismatch_map.insert(make_pair(edge, tmp));
		  vector<MismatchInfo> rc_tmp;
		  mismatch_map.insert(make_pair(g_.conjugate(edge), rc_tmp));
	  }
	  VERIFY(counts.size() == 4);
	  mismatch_map[edge].push_back(MismatchInfo(position, ratio, counts, cutoff));
	  vector<size_t> reversed_counts(4, 0);
	  for(size_t i = 0; i < 4; i ++)
		  reversed_counts[complement((unsigned char) i)] = counts[i];
	  VERIFY(reversed_counts.size() == 4);
	  mismatch_map[g_.conjugate(edge)].push_back(MismatchInfo(g_.length(edge) + g_.k() - position - 1, ratio, reversed_counts, cutoff));

  }
  string MaskedEdgeNucls(EdgeId edge, double cutoff) {
	  const Sequence& s_edge = g_.EdgeNucls(edge);
	  string s = s_edge.str();
	  if (mismatch_map[edge].size() > 0) {
		  DEBUG("in edge length " << g_.length(edge)<< " replacing " << mismatch_map[edge].size() << "mismatches");
	  }
	  for(size_t i = 0; i < mismatch_map[edge].size(); i++) {
		  DEBUG(mismatch_map[edge][i].position  << " " << mismatch_map[edge][i].ratio <<" " <<  mismatch_map[edge][i].counts.size());

		  if (mismatch_map[edge][i].ratio > cutoff && is_nucl(s[mismatch_map[edge][i].position])) {
			  DEBUG(s[mismatch_map[edge][i].position]  << "  before " );
			  VERIFY(mismatch_map[edge][i].position < s.length());
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
				  s[mismatch_map[edge][i].position] = (char) nucl((unsigned char) max_i);
				  DEBUG("replaced");
			  }
			  s[mismatch_map[edge][i].position] = char(s[mismatch_map[edge][i].position] +'a' - 'A');

			  if ((mismatch_map[edge][i].position >= 1 && ! is_nucl(s[mismatch_map[edge][i].position - 1])) || (mismatch_map[edge][i].position <s.length() -1  && ! is_nucl(s[mismatch_map[edge][i].position + 1]))){
//			  if (!cfg::get().mask_all && ((mismatch_map[edge][i].position >= 1 && ! is_nucl(s[mismatch_map[edge][i].position - 1])) || (mismatch_map[edge][i].position <s.length() -1  && ! is_nucl(s[mismatch_map[edge][i].position + 1])))){
				  ;
			  }else {
//				  INFO(mismatch_map[edge][i].cutoff <<" " << mismatch_map[edge][i].ratio <<" " << cfg::get().mismatch_ratio )
				  if (mismatch_map[edge][i].cutoff * cfg::get().mismatch_ratio < mismatch_map[edge][i].ratio )
					  s[mismatch_map[edge][i].position] = char('N');

			  }
			  DEBUG(s[mismatch_map[edge][i].position]  << "  " <<'a' - 'A');
		  }
	  }
	  return s;
  }
};

}


