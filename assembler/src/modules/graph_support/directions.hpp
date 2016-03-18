#pragma once

namespace omnigraph {
template<class Graph>
class AbstractDirection {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph &graph_;

protected:
    const Graph &graph() const {
        return graph_;
    }

public:
    AbstractDirection(const Graph &graph)
            : graph_(graph) { }

    virtual ~AbstractDirection() { }

    virtual const std::vector <EdgeId> OutgoingEdges(VertexId v) const = 0;

    virtual const std::vector <EdgeId> IncomingEdges(VertexId v) const = 0;

    virtual size_t OutgoingEdgeCount(VertexId v) const = 0;

    virtual size_t IncomingEdgeCount(VertexId v) const = 0;

    virtual VertexId EdgeStart(EdgeId edge) const = 0;

    virtual VertexId EdgeEnd(EdgeId edge) const = 0;

    bool CheckUniqueOutgoingEdge(VertexId v) const {
        return OutgoingEdgeCount(v) == 1;
    }

    EdgeId GetUniqueOutgoingEdge(VertexId v) const {
        return OutgoingEdges(v)[0];
    }

    bool CheckUniqueIncomingEdge(VertexId v) const {
        return IncomingEdgeCount(v) == 1;
    }

    EdgeId GetUniqueIncomingEdge(VertexId v) const {
        return IncomingEdges(v)[0];
    }

    virtual bool IsForward() const = 0;
};

template<class Graph>
class ForwardDirection : public AbstractDirection<Graph> {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
public:
    ForwardDirection(const Graph &graph)
            : AbstractDirection<Graph>(graph) {
    }

    virtual const std::vector <EdgeId> OutgoingEdges(VertexId v) const {
        return std::vector<EdgeId>(this->graph().out_begin(v), this->graph().out_end(v));
    }

    virtual const std::vector <EdgeId> IncomingEdges(VertexId v) const {
        return std::vector<EdgeId>(this->graph().in_begin(v), this->graph().in_end(v));
    }

    virtual size_t OutgoingEdgeCount(VertexId v) const {
        return this->graph().OutgoingEdgeCount(v);
    }

    virtual size_t IncomingEdgeCount(VertexId v) const {
        return this->graph().IncomingEdgeCount(v);
    }

    virtual VertexId EdgeStart(EdgeId edge) const {
        return this->graph().EdgeStart(edge);
    }

    virtual VertexId EdgeEnd(EdgeId edge) const {
        return this->graph().EdgeEnd(edge);
    }

    bool IsForward() const {
        return true;
    }
};

template<class Graph>
class BackwardDirection : public AbstractDirection<Graph> {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
public:
    BackwardDirection(const Graph &graph)
            : AbstractDirection<Graph>(graph) {
    }

    virtual const std::vector <EdgeId> OutgoingEdges(VertexId v) const {
        return std::vector<EdgeId>(this->graph().in_begin(v), this->graph().in_end(v));
    }

    virtual const std::vector <EdgeId> IncomingEdges(VertexId v) const {
        return std::vector<EdgeId>(this->graph().out_begin(v), this->graph().out_end(v));
    }

    virtual size_t OutgoingEdgeCount(VertexId v) const {
        return this->graph().IncomingEdgeCount(v);
    }

    virtual size_t IncomingEdgeCount(VertexId v) const {
        return this->graph().OutgoingEdgeCount(v);
    }

    virtual VertexId EdgeStart(EdgeId edge) const {
        return this->graph().EdgeEnd(edge);
    }

    virtual VertexId EdgeEnd(EdgeId edge) const {
        return this->graph().EdgeStart(edge);
    }

    bool IsForward() const {
        return false;
    }

};
}