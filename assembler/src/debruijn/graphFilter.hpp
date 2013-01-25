//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef GRAPHFILTER_H

#define GRAPHFILTER_H

#include "sequence/sequence.hpp"
//This graph will be a crucial component in multi-library
//TODO tdd this in the next iteration
template<class Graph>
class IGraphFilter
{
public:
    virtual ~IGraphFilter(){};
    /*
     * Filter all \weakly non eulerian paths.
     */
    virtual void Filter(Graph& graph, Graph& originalGraph) = 0;
};


template <class Graph>
class SimpleGraphFilter: public IGraphFilter<Graph>
{
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    bool IsIsolatedEdge(Graph& originalGraph,const Sequence& seq);
public:
    SimpleGraphFilter(){}
    void Filter(Graph& graph, Graph& originalGraph) ;
    ~SimpleGraphFilter(){};
};

template<class Graph>
void SimpleGraphFilter<Graph>::Filter(Graph& graph, Graph& originalGraph)
{
    for (auto iter = graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
        EdgeId edge = *iter;
        VertexId start = graph.EdgeStart(edge);
        VertexId end = graph.EdgeEnd(edge);
        if((graph.OutgoingEdgeCount(start) == 1) && (graph.IncomingEdgeCount(start) == 0) && (graph.IncomingEdgeCount(end) == 1) && (graph.OutgoingEdgeCount(end) == 0))
        {
            //If this edge is not an isolated edge in the de Bruijn graph, remove it.
            if(!IsIsolatedEdge(originalGraph, graph.EdgeNucls(edge)) ){

                graph.DeleteEdge(edge);
                graph.DeleteVertex(start);
                graph.DeleteVertex(end);
            }
            
        }
    }
}
template<class Graph>
bool SimpleGraphFilter<Graph>::IsIsolatedEdge(Graph& graph,const Sequence& seq)
{
    //I dont' know if there is a fast way to find the edgeID given a string, if there is, this func should be rewritten
    for(auto iter = graph.SmartEdgeBegin() ; !iter.IsEnd(); ++iter)
    {
        if( graph.EdgeNucls(*iter) == seq)
        {
            VertexId start = graph.EdgeStart(*iter);
            VertexId end = graph.EdgeEnd(*iter);
            if((graph.OutgoingEdgeCount(start) == 1) && (graph.IncomingEdgeCount(start) == 0) && (graph.IncomingEdgeCount(end) == 1) && (graph.OutgoingEdgeCount(end) == 0))
                return true;

            return false;
        }
    }
    INFO("ERROR!");
    return true;

}
#endif /* end of include guard: GRAPHFILTER_H */
