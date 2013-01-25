//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * path_set_tools.hpp
 *
 *  Created on: Oct 23, 2011
 *      Author: undead
 */
#pragma once
#include "path_set.hpp"
#include "graph_pack.hpp"

namespace debruijn_graph {
template <typename graph_pack>
string str(const PathSet<typename graph_pack::graph_t::EdgeId>& pathSet, const graph_pack& gp) {
	stringstream pathsString;
	size_t linecounter = 1;
	for(auto iter = pathSet.paths.begin() ; iter != pathSet.paths.end() ; ++iter)
	{
		pathsString << "Path " << linecounter <<":"<< pathSet.length<< " "<<  gp.int_ids.ReturnIntId(pathSet.start) <<"--" ;
		linecounter++;
		for(size_t i = 0 ; i < (*iter).size() ; ++i)
		{
			pathsString << gp.int_ids.ReturnIntId(((*iter)[i])) << " -- " ;
		}
		pathsString<<  gp.int_ids.ReturnIntId(pathSet.end);
		pathsString<<endl;
	}
	stringstream res;
	res << "id: "<< pathSet.id <<" weight "<< pathSet.weight <<" Start = " << gp.int_ids.ReturnIntId(pathSet.start) <<" ....... "<<"End = " << gp.int_ids.ReturnIntId(pathSet.end)<< endl<< pathsString.str() ;
	return res.str();
}




}

