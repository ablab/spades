#pragma once
#include "utils/logger/logger.hpp"

namespace omnigraph {

template<class Graph>
void RemoveIsolatedOrCompress(Graph& g, typename Graph::VertexId v) {
    if (g.IsDeadStart(v) && g.IsDeadEnd(v)) {
        g.DeleteVertex(v);
    } else {
        g.CompressVertex(v);
    }
}

template<class Graph>
class EdgeRemover {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef std::function<void(EdgeId)> HandlerF;

    Graph& g_;
    HandlerF removal_handler_;

public:
    EdgeRemover(Graph& g, HandlerF removal_handler = nullptr)
            : g_(g),
              removal_handler_(removal_handler) {
    }

    void DeleteEdge(EdgeId e) {
        VertexId start = g_.EdgeStart(e);
        VertexId end = g_.EdgeEnd(e);
        DeleteEdgeWithNoCompression(e);
        // NOTE: e here is already dead!
        TRACE("Compressing locality");
        if (!g_.RelatedVertices(start, end)) {
            TRACE("Vertices not related");
            TRACE("Processing end");
            RemoveIsolatedOrCompress(g_, end);
            TRACE("End processed");
        }
        TRACE("Processing start");
        RemoveIsolatedOrCompress(g_, start);
        TRACE("Start processed");
    }

    void DeleteEdgeWithNoCompression(EdgeId e) {
        TRACE("Deletion of edge " << g_.str(e));
        TRACE("Start " << g_.str(g_.EdgeStart(e)));
        TRACE("End " << g_.str(g_.EdgeEnd(e)));
        if (removal_handler_) {
            TRACE("Calling handler");
            removal_handler_(e);
        }
        TRACE("Deleting edge");
        g_.DeleteEdge(e);
    }

private:
    DECL_LOGGER("EdgeRemover");
};

//todo rewrite with SmartSetIterator
template<class Graph>
class ComponentRemover {
public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef std::function<void(const std::set<EdgeId>&)> HandlerF;

private:
    Graph& g_;
    HandlerF removal_handler_;

    template<class ElemType>
    void InsertIfNotConjugate(std::set<ElemType>& elems, ElemType elem) {
        if (elems.count(g_.conjugate(elem)) == 0) {
            elems.insert(elem);
        }
    }

public:
    ComponentRemover(Graph& g, HandlerF removal_handler = 0)
            : g_(g),
              removal_handler_(removal_handler) {
    }

    template<class EdgeIt>
    void DeleteComponent(EdgeIt begin, EdgeIt end, bool alter_vertices = true) {
        using std::set;
        set<EdgeId> edges;
        set<VertexId> vertices;

        //cleaning conjugates and gathering vertices
        for (EdgeIt it = begin; it != end; ++it) {
            EdgeId e = *it;
            InsertIfNotConjugate(edges, e);
            InsertIfNotConjugate(vertices, g_.EdgeStart(e));
            InsertIfNotConjugate(vertices, g_.EdgeEnd(e));
        }

        if (removal_handler_) {
            removal_handler_(edges);
        }

        for (EdgeId e: edges) {
            g_.DeleteEdge(e);
        }

        if (alter_vertices) {
            for (VertexId v: vertices) {
                RemoveIsolatedOrCompress(g_, v);
            }
        }
    }

    template<class Container>
    void DeleteComponent(const Container& container, bool alter_vertices = true) {
        DeleteComponent(container.begin(), container.end(), alter_vertices);
    }

};
}
