//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/dijkstra/neighbours_iterator.hpp"

using namespace debruijn_graph;

namespace dipspades {

class EdgeGluer {
    Graph &graph_;

    void MoveRelatedEdge(EdgeId edge, VertexId new_start, VertexId new_end){
        EdgeId new_edge = graph_.AddEdge(new_start, new_end, graph_.EdgeNucls(edge));
        TRACE("New edge " << graph_.str(new_edge) << "was added");
        graph_.DeleteEdge(edge);
    }

    void MoverUnrelatedEdge(EdgeId edge, VertexId new_start, VertexId new_end){
        EdgeId new_edge = graph_.AddEdge(new_start, new_end, graph_.EdgeNucls(edge));
        TRACE("New edge - " << graph_.str(new_edge) << " old edge - " << graph_.str(edge));
        if(IsEdgeRelated(graph_, new_edge))
            graph_.DeleteEdge(edge);
        else
            graph_.GlueEdges(edge, new_edge);
    }

    void StandardEdgeMoving(EdgeId edge, VertexId new_start, VertexId new_end){
        graph_.AddEdge(new_start, new_end, graph_.EdgeNucls(edge));
        graph_.DeleteEdge(edge);
    }

public:
    EdgeGluer(Graph &graph) : graph_(graph) { }

    void MoveEdgesFromVertexToVertex(VertexId old_vertex, VertexId new_vertex,
            vector<EdgeId> forbidden_edges){

        TRACE("New start - " << graph_.str(new_vertex) << ", old vertex - " << graph_.str(old_vertex));
        TRACE("Incoming edges");
        for(auto in_edges_iter = SmartSetIterator<Graph, EdgeId>(graph_,
                graph_.IncomingEdges(old_vertex).begin(),
                graph_.IncomingEdges(old_vertex).end());
                !in_edges_iter.IsEnd(); ++in_edges_iter){
            if(find(forbidden_edges.begin(), forbidden_edges.end(), *in_edges_iter) ==
                    forbidden_edges.end()){
                TRACE("Edge " << graph_.str(*in_edges_iter) << " is not forbidden");
                if(IsEdgeRelated(graph_, *in_edges_iter)){
                    TRACE("Edge is related");
                    if(IsEdgeLoop(graph_, *in_edges_iter)){
                        TRACE("Edge is loop");
                        StandardEdgeMoving(*in_edges_iter, new_vertex, new_vertex);
                    }
                    else{
                        TRACE("Edge is adjacent to conjugate");
                        StandardEdgeMoving(*in_edges_iter, graph_.conjugate(new_vertex), new_vertex);
                    }
                }
                else{
                    TRACE("Edge is not related");
                    StandardEdgeMoving(*in_edges_iter, graph_.EdgeStart(*in_edges_iter), new_vertex);
                }
            }
        }

        TRACE("Outgoing edges");
        for(auto out_edges_iter = SmartSetIterator<Graph, EdgeId>(graph_,
                graph_.OutgoingEdges(old_vertex).begin(),
                graph_.OutgoingEdges(old_vertex).end());
                !out_edges_iter.IsEnd(); ++out_edges_iter){
            if(find(forbidden_edges.begin(), forbidden_edges.end(), *out_edges_iter) ==
                    forbidden_edges.end()){
                TRACE("Edge " << graph_.str(*out_edges_iter) << " is not forbidden");
                if(IsEdgeRelated(graph_, *out_edges_iter)){
                    TRACE("Edge is related");
                    if(IsEdgeLoop(graph_, *out_edges_iter)){
                        TRACE("Edge is loop");
                        StandardEdgeMoving(*out_edges_iter, new_vertex, new_vertex);
                    }
                    else{
                        TRACE("Edge is adjacent to conjugate");
                        StandardEdgeMoving(*out_edges_iter, new_vertex, graph_.conjugate(new_vertex));
                    }
                }
                else{
                    TRACE("Edge is not related");
                    StandardEdgeMoving(*out_edges_iter, new_vertex, graph_.EdgeEnd(*out_edges_iter));
                }
            }
        }
    }

private:
    DECL_LOGGER("EdgeGluer");
};

}
