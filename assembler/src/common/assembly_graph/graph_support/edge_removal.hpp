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
        DeleteEdgeNoCompress(e);
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

    void DeleteEdgeNoCompress(EdgeId e) {
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

    void DeleteEdgeOptCompress(EdgeId e, bool compress) {
        if (compress)
            DeleteEdge(e);
        else
            DeleteEdgeNoCompress(e);
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

//Removes first 'trim_len' (k+1)-mers of graph edge, disconnecting it from starting vertex
//In case edge was removed, its end will be compressed even with "compress = false" parameter
template<class Graph>
class EdgeDisconnector {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    Graph& g_;
    EdgeRemover<Graph> edge_remover_;
    const size_t trim_len_;
    typedef std::function<void(EdgeId)> HandlerF;

public:
    EdgeDisconnector(Graph& g,
                     HandlerF removal_handler = nullptr,
                     size_t trim_len = 1):
            g_(g),
            edge_remover_(g, removal_handler),
            trim_len_(trim_len) {
        VERIFY(trim_len_ > 0);
    }

    EdgeId operator()(EdgeId e, bool compress = true) {
        if (g_.length(e) <= trim_len_
                || (e == g_.conjugate(e) && g_.length(e) <= 2 * trim_len_)) {
            VertexId start = g_.EdgeStart(e);
            VertexId end = g_.EdgeEnd(e);
            edge_remover_.DeleteEdgeOptCompress(e, compress);
            if (!compress && !g_.RelatedVertices(start, end)) {
                TRACE("Processing end");
                RemoveIsolatedOrCompress(g_, end);
                TRACE("End processed");
            }
            return EdgeId();
        } else {
            pair<EdgeId, EdgeId> split_res = g_.SplitEdge(e, trim_len_);
            edge_remover_.DeleteEdgeOptCompress(split_res.first, compress);
            return split_res.second;
        }
    }
};

}
