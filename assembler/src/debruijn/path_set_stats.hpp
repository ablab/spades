//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once
#include "graph_pack.hpp"

template<class graph_pack>
class PathSetStats {
    
    typedef typename graph_pack::graph_t::EdgeId EdgeId;
    typedef typename graph_pack::graph_t::VertexId VertexId;
    typedef typename graph_pack::graph_t Graph;
    typedef vector<EdgeId > Path;
    const Graph& g_;
    const PathSetIndexData<EdgeId> &rawData_;
    const PathSetIndexData<EdgeId> &filteredData_;
private:
void PathsNumbersStats()
{
    vector<size_t> rawVector;
    vector<size_t> filteredVector;

    for(auto iter = rawData_.begin() ; iter != rawData_.end() ; ++iter)
    {
        PathSet<EdgeId> pathset = *iter;
        rawVector.push_back( pathset.paths.size());
    }
    for(auto iter = filteredData_.begin() ; iter != filteredData_.end() ; ++iter)
    {
        PathSet<EdgeId> pathset = *iter;
        filteredVector.push_back(pathset.paths.size());
    }
    stringstream outputString;
     WriteVector(rawVector, "raw", outputString );
     WriteVector(filteredVector, "filter" ,outputString);
     outputString<< " h = hist ([ raw, filter ]); ha = bar(h); legend(ha, 'Before Processing Pathset', 'After Processing Pathset' ) ; ";

     //todo write this out to *.m

}
void WriteVector(vector<size_t> v, const char* nameVector,stringstream& outputString)
{
    outputString<< nameVector;
    outputString<< " = [ ";
    for(size_t i = 0 ; i < v.size() ; ++i)
    {
        if(i != v.size() -1)
            outputString<< v[i] << ", " ;
        else
            outputString<< v[i] << " ]; " ;
    }
}
public:

	PathSetStats(const graph_pack& gp, const PathSetIndexData<EdgeId> &rawData, const PathSetIndexData<EdgeId> &filteredData):
        g_(gp.g), rawData_(rawData), filteredData_(filteredData)
    {
    }

	void Count() {
        PathsNumbersStats();
    }


};
