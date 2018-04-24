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

    ConstructionHelper(Graph &graph) : graph_(graph) {
    }

    Graph &graph() {
        return graph_;
    }

    EdgeId AddEdge(const EdgeData &data) {
        return AddEdge(data, graph_.GetGraphIdDistributor());
    }

    EdgeId AddEdge(const EdgeData &data, restricted::IdDistributor &id_distributor) {
        return graph_.AddEdge(data, id_distributor);
    }

    void LinkIncomingEdge(VertexId v, EdgeId e) {
        VERIFY(graph_.EdgeEnd(e) == VertexId());
        graph_.conjugate(v)->AddOutgoingEdge(graph_.conjugate(e));
        e->SetEndVertex(v);
    }

    void LinkOutgoingEdge(VertexId v, EdgeId e) {
        VERIFY(graph_.EdgeEnd(graph_.conjugate(e)) == VertexId());
        v->AddOutgoingEdge(e);
        graph_.conjugate(e)->SetEndVertex(graph_.conjugate(v));
    }

    void LinkEdges(EdgeId e1, EdgeId e2) {
        VertexId v = graph_.EdgeEnd(e1), w = graph_.EdgeStart(e2);
        DeleteLink(w, e2);
        LinkOutgoingEdge(v, e2);
    }

    void DeleteLink(VertexId v, EdgeId e) {
        bool res = v->RemoveOutgoingEdge(e);
        VERIFY(res);
        graph_.conjugate(e)->SetEndVertex(VertexId());
    }

    void DeleteUnlinkedEdge(EdgeId e) {
        EdgeId rc = graph_.conjugate(e);
        if (e != rc) {
            delete rc.get();
        }
        delete e.get();
    }

    void DeleteUnlinkedVertex(VertexId v) {
        VertexId rc = graph_.conjugate(v);
        delete rc.get(); // These guys do check that everything is unlinked.
        delete v.get();
    }

    VertexId CreateVertex(const VertexData &data) {
        return CreateVertex(data, graph_.GetGraphIdDistributor());
    }

    VertexId CreateVertex(const VertexData &data, restricted::IdDistributor &id_distributor) {
        return graph_.CreateVertex(data, id_distributor);
    }

    template<class Iter>
    void AddVerticesToGraph(Iter begin, Iter end) {
        for(; begin != end; ++begin) {
            graph_.AddVertexToGraph(*begin);
        }
    }
};

}
