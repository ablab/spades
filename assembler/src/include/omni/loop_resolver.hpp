//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef LOOP_RESOLVER_H
#define LOOP_RESOLVER_H

#include <math.h>

namespace omnigraph{

template<class Graph>
class LoopResolver{
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph &graph_;
    double allowCoverageVatiation_;
public:
    LoopResolver(Graph& graph, double allowCoverageVatiation)
        :graph_(graph),allowCoverageVatiation_(allowCoverageVatiation){}
    void ResolveLoops()
    {
        for(auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter)
        {
//            graph_.DeleteEdge(*iter);
        
            if(graph_.EdgeStart(*iter) == graph_.EdgeEnd(*iter))
            {
                INFO("LOOP RE:"<< graph_.length(*iter));
                VertexId loopNode = graph_.EdgeStart(*iter);
                //if the loop appears on a \emph{simple} path.
                if((graph_.OutgoingEdgeCount(loopNode) == 2) && (graph_.IncomingEdgeCount(loopNode)==2) )
                {
                    vector<EdgeId> inComingEdges = graph_.IncomingEdges(loopNode);
                    EdgeId beforeLoopEdge = (inComingEdges[0] == *iter)? inComingEdges[1]: inComingEdges[0]; 
                    vector<EdgeId> outGoingEdges = graph_.OutgoingEdges(loopNode);
                    EdgeId afterLoopEdge = (outGoingEdges[0] == *iter) ? outGoingEdges[1]: outGoingEdges[0];
                    double loopCoverage = graph_.coverage(*iter);
                    double beforeLoopCoverage = graph_.coverage(beforeLoopEdge);
                    double afterLoopCoverage = graph_.coverage(afterLoopEdge);

                    double variance = (afterLoopCoverage + beforeLoopCoverage)/2.0 - loopCoverage ; 

                    INFO("LOOP RE D:"<<variance <<":"<< loopCoverage);
                    if( fabs( variance) <= allowCoverageVatiation_*loopCoverage)
                    {
                        INFO("LOOP RESO D:"<<variance <<":"<< loopCoverage);
                        VertexId addedVertex = graph_.AddVertex();
                        VertexId afterLoopEdgeEndVertex = graph_.EdgeEnd(afterLoopEdge);
                        EdgeId addedAfterLoopEdge =  graph_.AddEdge( addedVertex, afterLoopEdgeEndVertex, graph_.EdgeNucls(afterLoopEdge));
                        graph_.coverage_index().SetCoverage(addedAfterLoopEdge, graph_.coverage(afterLoopEdge) * graph_.length(afterLoopEdge) );
                        graph_.DeleteEdge(afterLoopEdge);

                        EdgeId transformedLoopEdge = graph_.AddEdge(loopNode, addedVertex, graph_.EdgeNucls(*iter));
                        graph_.coverage_index().SetCoverage(transformedLoopEdge, graph_.coverage(*iter) * graph_.length(*iter));
                        graph_.DeleteEdge(*iter);

                    }

                }
            }
        }
    }


};

}

#endif /* end of include guard: LOOP_RESOLVER_H */

