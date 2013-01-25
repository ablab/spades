//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "compare_standard.hpp"
#include "coloring.hpp"

namespace cap {

template<class Graph>
struct cap_graph_pack {
	typedef Graph graph_t;
	typedef string contig_id_t;
	typedef typename Graph::EdgeId EdgeId;
	Graph g;
	IdTrackHandler<Graph> int_ids;
	ColorHandler<Graph> coloring;
//	map<contig_id_t, vector<EdgeId>> red_paths;
//	map<contig_id_t, vector<EdgeId>> blue_paths;
	EdgesPositionHandler<Graph> edge_pos;

	cap_graph_pack(size_t k) :
			g(k), int_ids(g), coloring(g), edge_pos(g) {

	}
};

}
