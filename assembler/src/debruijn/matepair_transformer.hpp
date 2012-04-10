//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once
#include "omni/paired_info.hpp"
#include "path_set.hpp"
#include "path_set_tools.hpp"
#include "xmath.h"
using namespace math;

namespace omnigraph{


template<class graph_pack>
class MatePairTransformer{

typedef typename graph_pack::graph_t Graph;
typedef typename Graph::EdgeId EdgeId;
typedef vector<EdgeId> Path;
const graph_pack& gp;
const PairedInfoIndex<Graph>& pair_info_;

public:
MatePairTransformer(const graph_pack& gp, const PairedInfoIndex<Graph>& pair_info):
    gp(gp), pair_info_(pair_info){}


void Transform(PathSetIndexData<EdgeId> & pathset_index)
{
	for (auto it = pair_info_.begin(); it != pair_info_.end(); ++it) {

		vector < PairInfo < EdgeId >> infos = *it;

        for(auto iter = infos.begin() ; iter != infos.end() ; ++iter)
        {
            if (gr(iter->d, 0.)) {
                INFO( *iter);
                INFO(gp.g.length(iter->first)<<" : " << gp.g.length(iter->second));
                if (iter->variance == 0) {

                    PathStorageCallback<Graph> call_back(gp.g);
                    EdgeId first_edge = iter->first;
                    EdgeId second_edge = iter->second;
                    PathProcessor<Graph> path_processor(
                            gp.g,
                            iter->d - gp.g.length(first_edge),
                            iter->d - gp.g.length(first_edge),
                            gp.g.EdgeEnd(first_edge),
                            gp.g.EdgeStart(second_edge),
                            call_back);
                    path_processor.Process();
                    PathSet<EdgeId> pathset(first_edge, second_edge, iter->d + gp.g.length(second_edge) , call_back.paths(), iter->weight);

                    PathSetFilter(pathset);
                    int res_id = pathset_index.AddPathSet(pathset);
                    INFO("Has id " << res_id);
                }
                else
                {
                    INFO("NON ZERO VARIATION " << *iter << "edges" << iter->first << " " << iter->second);
                }
            }
        }
	}

}


///TODO  Remove invalid paths in this path_set
//
void PathSetFilter(PathSet<EdgeId> & path_set)
{
}

};
}
