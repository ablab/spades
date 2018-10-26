//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
//#include "core.hpp"
#include "observable_graph.hpp"

namespace omnigraph {

template<class DataMaster>
class ConstructionHelper {
    //typedef GraphCore<DataMaster> Graph;
    typedef ObservableGraph<DataMaster> Graph;
    typedef typename Graph::DataMasterT DataMasterT;
    typedef typename Graph::VertexData VertexData;
    typedef typename Graph::EdgeData EdgeData;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::VertexIt VertexIt;
    typedef typename Graph::edge_const_iterator edge_const_iterator;

    Graph &graph_;

public:

    ConstructionHelper(Graph &graph)
            : graph_(graph) {}

    Graph &graph() {
        return graph_;
    }

    EdgeId AddEdge(const EdgeData &data, EdgeId id = 0) {
        return graph_.AddEdge(data, id);
    }

    void LinkIncomingEdge(VertexId v, EdgeId e) {
        VERIFY(graph_.EdgeEnd(e) == VertexId());
        graph_.cvertex(v)->AddOutgoingEdge(graph_.conjugate(e));
        graph_.edge(e)->SetEndVertex(v);
    }

    void LinkOutgoingEdge(VertexId v, EdgeId e) {
        VERIFY(graph_.EdgeEnd(graph_.conjugate(e)) == VertexId());
        graph_.vertex(v)->AddOutgoingEdge(e);
        graph_.cedge(e)->SetEndVertex(graph_.conjugate(v));
    }

    void LinkEdges(EdgeId e1, EdgeId e2) {
        VertexId v = graph_.EdgeEnd(e1), w = graph_.EdgeStart(e2);
        DeleteLink(w, e2);
        LinkOutgoingEdge(v, e2);
    }

    void DeleteLink(VertexId v, EdgeId e) {
        bool res = graph_.vertex(v)->RemoveOutgoingEdge(e);
        VERIFY(res);
        graph_.cedge(e)->SetEndVertex(VertexId());
    }

    void DeleteUnlinkedEdge(EdgeId e) {
        EdgeId rc = graph_.conjugate(e);
        graph_.DestroyEdge(e, rc);
    }

    void DeleteUnlinkedVertex(VertexId v) {
        graph_.DestroyVertex(v);
    }

    VertexId CreateVertex(const VertexData &data, VertexId id = 0) {
        return graph_.CreateVertex(data, id);
    }
};

}
