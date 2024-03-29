//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
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


    void MaybeAddVertex(Graph &clone,
                        VertexData data,
                        VertexId id1, VertexId id2) const {
        if (clone.contains(id1))
            return;
        TRACE("Vertex " << id1 << " ~ " << id2 << " .");
        auto new_id = clone.AddVertex(std::move(data), id1, id2);
        VERIFY(new_id == id1);
        VERIFY(clone.conjugate(new_id) == id2);
    }

public:
    ConstructionHelper(Graph &graph)
            : graph_(graph) {}

    Graph &graph() {
        return graph_;
    }

    const DataMaster& master() const noexcept { return graph_.master_; }
    DataMaster& master() noexcept { return graph_.master_; }

    Graph clone() const {
        return component(graph_.e_begin(), graph_.e_end());
    }

    template<class It>
    Graph component(It start, It end) const {
        Graph clone(graph_.master_);
        clone.reserve(graph_.vreserved(), graph_.ereserved());

        for (; start != end; ++start) {
            EdgeId e = *start;
            if (clone.contains(e))
                continue;

            VertexId v1 = graph_.EdgeStart(e), v2 = graph_.EdgeEnd(e);
            VertexData vdata = graph_.data(v1);
            MaybeAddVertex(clone, vdata, v2, v2);

            EdgeId new_id = clone.AddEdge(v1, v2,
                                          graph_.data(e),
                                          e, graph_.conjugate(e));
            VERIFY(new_id == e);
            VERIFY(clone.conjugate(new_id) = graph_.conjugate(e));
        }

        return clone;
    }

    EdgeId AddEdge(EdgeData data, EdgeId id = 0) {
        return graph_.AddEdge(std::move(data), id);
    }

    EdgeId AddEdge(EdgeData data, EdgeId id, EdgeId cid) {
        return graph_.AddEdge(std::move(data), id, cid);
    }

    void LinkIncomingEdge(VertexId v, EdgeId e) {
        VERIFY(graph_.EdgeEnd(e) == VertexId());
        graph_.cvertex(v).AddOutgoingEdge(graph_.conjugate(e));
        graph_.edge(e).SetEndVertex(v);
    }

    void LinkOutgoingEdge(VertexId v, EdgeId e) {
        VERIFY(graph_.EdgeEnd(graph_.conjugate(e)) == VertexId());
        graph_.vertex(v).AddOutgoingEdge(e);
        graph_.cedge(e).SetEndVertex(graph_.conjugate(v));
    }

    void LinkEdges(EdgeId e1, EdgeId e2) {
        VertexId v = graph_.EdgeEnd(e1), w = graph_.EdgeStart(e2);
        std::vector<EdgeId> in_edges;
        std::vector<EdgeId> out_edges;
        auto in_edges_it = graph_.IncomingEdges(w);
        auto out_edges_it = graph_.OutgoingEdges(w);
        std::copy(in_edges_it.begin(), in_edges_it.end(), std::back_inserter(in_edges));
        std::copy(out_edges_it.begin(), out_edges_it.end(), std::back_inserter(out_edges));
        for (auto in_edge: in_edges) {
            DeleteLink(graph_.conjugate(w), graph_.conjugate(in_edge));
            LinkOutgoingEdge(graph_.conjugate(v), graph_.conjugate(in_edge));
        }
        for (auto out_edge: out_edges) {
            DeleteLink(w, out_edge);
            LinkOutgoingEdge(v, out_edge);
        }
    }

    void DeleteLink(VertexId v, EdgeId e) {
        bool res = graph_.vertex(v).RemoveOutgoingEdge(e);
        VERIFY(res);
        graph_.cedge(e).SetEndVertex(VertexId());
    }

    void DeleteUnlinkedEdge(EdgeId e) {
        EdgeId rc = graph_.conjugate(e);
        graph_.DestroyEdge(e, rc);
    }

    void DeleteUnlinkedVertex(VertexId v) {
        graph_.DestroyVertex(v);
    }

    VertexId CreateVertex(VertexData data, VertexId id = 0) {
        return graph_.CreateVertex(std::move(data), id);
    }

    VertexId CreateVertex(VertexData data, VertexId id1, VertexId id2) {
        return graph_.CreateVertex(std::move(data), id1, id2);
    }

};

}
