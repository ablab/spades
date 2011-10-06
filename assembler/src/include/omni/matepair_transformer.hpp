#pragma once
#include "paired_info.hpp"
#include "path_set.hpp"
#include "xmath.h"
using namespace math;

namespace omnigraph{


template<class Graph>
class MatePairTransformer{


typedef typename Graph::EdgeId EdgeId;
typedef vector<EdgeId > Path;
const Graph& g_;
const PairedInfoIndex<Graph>& pair_info_;

public:
MatePairTransformer(const Graph& g, const PairedInfoIndex<Graph>& pair_info):
    g_(g), pair_info_(pair_info){}


void Transform(PathSetIndexData<EdgeId> & pathset_index)
{
	for (auto it = pair_info_.begin(); it != pair_info_.end(); ++it) {

		vector < PairInfo < EdgeId >> infos = *it;

        for(auto iter = infos.begin() ; iter != infos.end() ; ++iter)
        {
            if (gr(iter->d, 0.)) {
                INFO( *iter);
                if (iter->variance == 0) {

                    PathReceiverCallback<Graph> call_back(g_); 
                    EdgeId first_edge = iter->first;
                    EdgeId second_edge = iter->second;
                    PathProcessor<Graph> path_processor(
                            g_,
                            iter->d - g_.length(first_edge), 
                            iter->d - g_.length(first_edge),
                            g_.EdgeEnd(first_edge),
                            g_.EdgeStart(second_edge),
                            call_back);
                    path_processor.Process();
                    PathSet<EdgeId> pathset(first_edge, second_edge, iter->d + g_.length(second_edge) , call_back.paths());
                    PathSetFilter(pathset);
                    pathset_index.AddPathSet(pathset);
                }
                else
                {
                    INFO("NON ZERO VARIATION " << *iter);
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
