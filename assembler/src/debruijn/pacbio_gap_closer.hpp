/*
 * pac_index.hpp
 *
 *  Created on: Jan 21, 2013
 *      Author: lab42
 */
#pragma once

#include "debruijn_kmer_index.hpp"
#include "graph_pack.hpp"
#include <algorithm>
#include "pacbio_read_structures.hpp"
#include "long_read_storage.hpp"

template<class Graph>
class PacbioGapCloser {
	typedef typename Graph::EdgeId EdgeId;
private:
	DECL_LOGGER("PacIndex");
	Graph &g_;
	map<EdgeId, map<EdgeId, pair<size_t, string> > > new_edges;
	PacbioGapCloser(Graph &g):g_(g){
	}

};
