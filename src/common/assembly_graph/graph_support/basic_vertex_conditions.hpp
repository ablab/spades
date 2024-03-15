#pragma once
#include "func/pred.hpp"

namespace omnigraph {

template<class Graph>
class VertexCondition : public func::AbstractPredicate<typename Graph::VertexId> {
    typedef typename Graph::VertexId VertexId;
    const Graph &g_;
protected:

    VertexCondition(const Graph &g)
            : g_(g) {
    }

    const Graph &g() const {
        return g_;
    }

};

template<class Graph>
class CompressCondition : public VertexCondition<Graph> {
    typedef typename Graph::VertexId VertexId;

public:
    CompressCondition(const Graph &g) :
            VertexCondition<Graph>(g) {
    }

    bool Check(VertexId v) const override {
        return this->g().CanCompressVertex(v);
    }
};

template<class Graph>
class IsolatedVertexCondition : public VertexCondition<Graph> {
    typedef typename Graph::VertexId VertexId;

public:
    IsolatedVertexCondition(const Graph& g) :
            VertexCondition<Graph>(g) {
    }

    bool Check(VertexId v) const override {
        return this->g().IsDeadStart(v) && this->g().IsDeadEnd(v);
    }
};

template<class Graph>
class TerminalVertexCondition : public VertexCondition<Graph> {
    typedef typename Graph::VertexId VertexId;

public:
    TerminalVertexCondition(const Graph& g) :
            VertexCondition<Graph>(g) {
    }

    bool Check(VertexId v) const override {
        return this->g().IncomingEdgeCount(v) + this->g().OutgoingEdgeCount(v) == 1;
    }

};

}