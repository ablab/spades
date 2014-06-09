#pragma once

#include "func.hpp"
#include "graph_component.hpp"
#include "coverage.hpp"

namespace omnigraph {

template<class Graph>
class ProcessingAlgorithm : private boost::noncopyable {
    Graph& g_;

 protected:
    Graph& g() {
        return g_;
    }

    const Graph& g() const {
        return g_;
    }

 public:
    ProcessingAlgorithm(Graph& g)
            : g_(g) {

    }

    virtual ~ProcessingAlgorithm() {
    }

    virtual bool Process() = 0;
};

template<class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class EdgeProcessingAlgorithm : public ProcessingAlgorithm<Graph> {
    typedef ProcessingAlgorithm<Graph> base;
    typedef typename Graph::EdgeId EdgeId;

    const Comparator comp_;
    const shared_ptr<func::Predicate<EdgeId>> proceed_condition_;

 protected:
    //todo make private
    virtual bool ProcessEdge(EdgeId e) = 0;

 public:
    EdgeProcessingAlgorithm(
            Graph& g,
            const Comparator& c = Comparator(),
            shared_ptr<func::Predicate<EdgeId>> proceed_condition = make_shared<
                    func::AlwaysTrue<EdgeId>>())
            : base(g),
              comp_(c),
              proceed_condition_(proceed_condition) {

    }

    bool Process() {
        TRACE("Start processing");
        bool triggered = false;
        for (auto it = this->g().SmartEdgeBegin(comp_); !it.IsEnd(); ++it) {
            if (!proceed_condition_->Check(*it)) {
                TRACE("Stop condition was reached.");
                break;
            }

            TRACE("Processing edge " << this->g().str(*it));
            triggered |= ProcessEdge(*it);
        }
        TRACE("Finished processing. Triggered = " << triggered);
        return triggered;
    }

 private:
    DECL_LOGGER("EdgeProcessingAlgorithm")
    ;
};

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
    typedef boost::function<void(EdgeId)> HandlerF;

    Graph& g_;
    HandlerF removal_handler_;

 public:
    EdgeRemover(Graph& g, HandlerF removal_handler = 0)
            : g_(g),
              removal_handler_(removal_handler) {
    }

    //todo how is it even compiling with const?!!!
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
    DECL_LOGGER("EdgeRemover")
    ;
};

template<class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class EdgeRemovingAlgorithm : public EdgeProcessingAlgorithm<Graph, Comparator> {
    typedef EdgeProcessingAlgorithm<Graph, Comparator> base;
    typedef typename Graph::EdgeId EdgeId;

    shared_ptr<func::Predicate<EdgeId>> remove_condition_;
    EdgeRemover<Graph> edge_remover_;

 protected:
    bool ProcessEdge(EdgeId e) {
        if (remove_condition_->Check(e)) {
            edge_remover_.DeleteEdge(e);
            return true;
        }
        return false;
    }

 public:
    EdgeRemovingAlgorithm(
            Graph& g,
            shared_ptr<func::Predicate<EdgeId>> remove_condition,
            boost::function<void(EdgeId)> removal_handler = boost::none,
            const Comparator& c = Comparator(),
            shared_ptr<func::Predicate<EdgeId>> proceed_condition = make_shared<
                    func::AlwaysTrue<EdgeId>>())
            : base(g, c, proceed_condition),
              remove_condition_(remove_condition),
              edge_remover_(g, removal_handler) {

    }

 private:
    DECL_LOGGER("EdgeRemovingAlgorithm")
    ;
};

//todo rewrite with SmartSetIterator
template<class Graph>
class ComponentRemover {
 public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef boost::function<void(const set<EdgeId>&)> HandlerF;

 private:
    Graph& g_;
    HandlerF removal_handler_;

    template<class ElemType>
    void InsertIfNotConjugate(set<ElemType>& elems, ElemType elem) {
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

        FOREACH (EdgeId e, edges) {
            g_.DeleteEdge(e);
        }

        if (alter_vertices) {
            FOREACH (VertexId v, vertices) {
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
