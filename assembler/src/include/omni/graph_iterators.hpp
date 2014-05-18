#ifndef __OMNI_GRAPH_ITERATORS_HPP__
#define __OMNI_GRAPH_ITERATORS_HPP__

#include "adt/queue_iterator.hpp"
#include "io/read_processor.hpp"
#include "omni/action_handlers.hpp"

#include <boost/iterator/iterator_facade.hpp>

namespace omnigraph {

/**
 * SmartIterator is abstract class which acts both as QueueIterator and GraphActionHandler. As QueueIterator
 * SmartIterator is able to iterate through collection content of which can be changed in process of
 * iteration. And as GraphActionHandler SmartIterator can change collection contents with respect to the
 * way graph is changed. Also one can define order of iteration by specifying Comparator.
 */
template<class Graph, typename ElementId, typename Comparator = std::less<
                                              ElementId> >
class SmartIterator : public GraphActionHandler<Graph>, public QueueIterator<
    ElementId, Comparator> {
  public:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
  private:
    bool add_new_;
  public:
    SmartIterator(const Graph &graph, const string &name, bool add_new,
                  const Comparator& comparator = Comparator())
            : GraphActionHandler<Graph>(graph, name),
              QueueIterator<ElementId, Comparator>(comparator),
              add_new_(add_new) {
    }

    virtual ~SmartIterator() {}

    virtual void HandleAdd(ElementId v) {
        if (add_new_)
            this->push(v);
    }

    virtual void HandleDelete(ElementId v) {
        this->erase(v);
    }
};

/**
 * SmartIterator is abstract class which acts both as QueueIterator and GraphActionHandler. As QueueIterator
 * SmartIterator is able to iterate through collection content of which can be changed in process of
 * iteration. And as GraphActionHandler SmartIterator can change collection contents with respect to the
 * way graph is changed. Also one can define order of iteration by specifying Comparator.
 */
template<class Graph, typename ElementId,
         typename Comparator = std::less<ElementId> >
class SmartSetIterator : public SmartIterator<Graph, ElementId, Comparator> {
  public:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
  public:
    template<class Iterator>
    SmartSetIterator(const Graph &graph, Iterator begin, Iterator end,
                     const Comparator& comparator = Comparator())
            : SmartIterator<Graph, ElementId, Comparator>(
                graph, "SmartSet " + ToString(this), false, comparator) {
        this->insert(begin, end);
    }

    virtual ~SmartSetIterator() {}
};

/*
 * ConditionedSmartSetIterator acts much like SmartSetIterator, but, unlike the above case,
 * one can (and must) provide merge handler that will decide whether to add merged edge to
 * the set being iterated or not (extending add_new_ parameter logic of SmartIterator)
 * Also has the ability to be `reset` (i.e. start from the begin-iterator with respect to
 * added and deleted values)
 * MergeHandler class/struct must provide:
 *  bool operator()(const std::vector<ElementId> &, ElementId)
 */
template<class Graph, typename ElementId, class MergeHandler>
class ConditionedSmartSetIterator : public SmartSetIterator<Graph, ElementId> {
  public:
    template <class Iterator>
    ConditionedSmartSetIterator(const Graph &graph, Iterator begin, Iterator end,
                                MergeHandler &merge_handler)
            : SmartSetIterator<Graph, ElementId>(graph, begin, end),
              merge_handler_(merge_handler),
              true_elements_() {

        for (auto it = begin; it != end; ++it) {
            true_elements_.insert(*it);
        }
    }

    virtual ~ConditionedSmartSetIterator() {}

    virtual void HandleAdd(ElementId v) {
        TRACE("handleAdd " << this->g().str(v));
        if (true_elements_.count(v)) {
            this->push(v);
        }
    }

    virtual void HandleDelete(ElementId v) {
        TRACE("handleDel " << this->g().str(v));
        super::HandleDelete(v);
        true_elements_.erase(v);
    }

    virtual void HandleMerge(const std::vector<ElementId>& old_edges, ElementId new_edge) {
        TRACE("handleMer " << this->g().str(new_edge));
        if (merge_handler_(old_edges, new_edge)) {
            true_elements_.insert(new_edge);
        }
    }

    virtual void reset() {
        TRACE("reset");
        this->operator++();
        this->insert(true_elements_.begin(), true_elements_.end());
    }

  private:
    typedef SmartSetIterator<Graph, ElementId> super;

    MergeHandler &merge_handler_;
    std::unordered_set<ElementId> true_elements_;

    DECL_LOGGER("ConditionedSmartSetIterator");
};

/**
 * SmartVertexIterator iterates through vertices of graph. It listens to AddVertex/DeleteVertex graph events
 * and correspondingly edits the set of vertices to iterate through. Note: high level event handlers are
 * triggered before low level event handlers like H>andleAdd/HandleDelete. Thus if Comparator uses certain
 * structure which is also updated with handlers make sure that all information is updated in high level
 * event handlers.
 */
template<class Graph, typename Comparator = std::less<typename Graph::VertexId> >
class SmartVertexIterator : public SmartIterator<Graph,
                                                 typename Graph::VertexId, Comparator> {
  public:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    static size_t get_id() {
        static size_t id = 0;
        return id++;
    }

  public:
    SmartVertexIterator(const Graph &graph, const Comparator& comparator =
                        Comparator())
            : SmartIterator<Graph, VertexId, Comparator>(
                graph, "SmartVertexIterator " + ToString(get_id()), true,
                comparator) {
        this->insert(graph.begin(), graph.end());
    }

    virtual ~SmartVertexIterator() {}
};

/**
 * SmartEdgeIterator iterates through edges of graph. It listens to AddEdge/DeleteEdge graph events
 * and correspondingly edits the set of edges to iterate through. Note: high level event handlers are
 * triggered before low level event handlers like HandleAdd/HandleDelete. Thus if Comparator uses certain
 * structure which is also updated with handlers make sure that all information is updated in high level
 * event handlers.
 */
template<class Graph, typename Comparator = std::less<typename Graph::EdgeId> >
class SmartEdgeIterator : public SmartIterator<Graph, typename Graph::EdgeId,
                                               Comparator> {
  public:
    typedef QueueIterator<typename Graph::EdgeId, Comparator> base;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    static size_t get_id() {
        static size_t id = 0;
        return id++;
    }
  public:
    SmartEdgeIterator(const Graph &graph, Comparator comparator = Comparator(),
                      vector<EdgeId>* edges = 0)
            : SmartIterator<Graph, EdgeId, Comparator>(
                graph, "SmartEdgeIterator " + ToString(get_id()), true,
                comparator) {
        if (edges == 0) {
            for (auto it = graph.begin(); it != graph.end(); ++it) {
                //todo: this solution doesn't work with parallel simplification
                this->base::insert(graph.out_begin(*it), graph.out_end(*it));
                //this does
                //auto out = graph.OutgoingEdges(*it);
                //this->base::insert(out.begin(), out.end());
            }
        } else {
            this->base::insert(edges->begin(), edges->end());
        }
    }
};

//todo return verifies when they can be switched off
template<class Graph>
class GraphEdgeIterator : public boost::iterator_facade<GraphEdgeIterator<Graph>,
        typename Graph::EdgeId, boost::forward_traversal_tag, typename Graph::EdgeId> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexIt const_vertex_iterator;
    typedef typename Graph::edge_const_iterator const_edge_iterator;
public:

    explicit GraphEdgeIterator(const Graph& g, const_vertex_iterator v_it)
            : g_(g),
              v_it_(v_it) {
    	if (v_it_ != g_.end()) {
    	    e_it_ = g_.out_begin(*v_it_);
            Skip();
    	}
    }

private:
    friend class boost::iterator_core_access;

    void Skip() {
        //VERIFY(v_it_ != g_.end());
        while (e_it_ == g_.out_end(*v_it_)) {
    	    v_it_++;
            if (v_it_ == g_.end()) 
                return;
            e_it_ = g_.out_begin(*v_it_);
    	}
    }

    void increment() {
        if (v_it_ == g_.end()) 
          return;
        e_it_++;
        Skip();
    }

    bool equal(const GraphEdgeIterator &other) const {
    	if (other.v_it_ != v_it_)
    		return false;
        if (v_it_ != g_.end() && other.e_it_ != e_it_)
        	return false;
        return true;
    }

    EdgeId dereference() const {
        //VERIFY(v_it_ != g_.end());
        return *e_it_;
    }

    const Graph& g_;
    const_vertex_iterator v_it_;
    const_edge_iterator e_it_;
};

template<class Graph>
class ConstEdgeIterator {
    typedef typename Graph::EdgeId EdgeId;
    GraphEdgeIterator<Graph> begin_, end_;

  public:
    ConstEdgeIterator(const Graph &g)
            : begin_(g, g.begin()), end_(g, g.end()) {
    }

    bool IsEnd() const {
        return begin_ == end_;
    }

    EdgeId operator*() const {
        return *begin_;
    }

    const ConstEdgeIterator& operator++() {
        begin_++;
        return *this;
    }
};

//template<class Graph>
//class ConstEdgeIterator :
//        public boost::iterator_facade<ConstEdgeIterator<Graph>,
//                                      typename Graph::EdgeId const,
//                                      boost::forward_traversal_tag,
//                                      typename Graph::EdgeId const> {
//  public:
//    ConstEdgeIterator(const Graph &g)
//            : graph_(g),
//              cvertex_(g.begin()), evertex_(g.end()),
//              cedge_(g.out_begin(*cvertex_)), eedge_(g.out_end(*cvertex_)) {
//        skip_empty();
//    }
//
//    bool IsEnd() const {
//        return cvertex_ == evertex_;
//    }
//
//  private:
//    friend class boost::iterator_core_access;
//
//    void skip_empty() {
//        while (cedge_ == eedge_) {
//            if (++cvertex_ == evertex_)
//                break;
//            cedge_ = graph_.out_begin(*cvertex_);
//            eedge_ = graph_.out_end(*cvertex_);
//        }
//    }
//
//    void increment() {
//        ++cedge_;
//        skip_empty();
//    }
//
//    bool equal(ConstEdgeIterator &other) const {
//        return (graph_ == other.graph_ &&
//                cvertex_ == other.cvertex_ &&
//                cedge_ == other.cedge_);
//    }
//
//    typename Graph::EdgeId const dereference() const {
//        return *cedge_;
//    }
//
//    const Graph &graph_;
//    typename Graph::VertexIt cvertex_, evertex_;
//    typename Graph::edge_const_iterator cedge_, eedge_;
//};

template<class Graph>
class ParallelEdgeProcessor {
    class ConstEdgeIteratorWrapper {
      public:
        typedef typename Graph::EdgeId ReadT;

        ConstEdgeIteratorWrapper(const Graph &g)
                : it_(g) {}

        bool eof() const { return it_.IsEnd(); }

        ConstEdgeIteratorWrapper& operator>>(typename Graph::EdgeId &val) {
            val = *it_;
            ++it_; 
            return *this;
        }

      private:
        ConstEdgeIterator<Graph> it_;
    };

  public:
    ParallelEdgeProcessor(const Graph &g, unsigned nthreads)
            : rp_(nthreads), it_(g) {}

    template <class Processor>
    bool Run(Processor &op) { return rp_.Run(it_, op); }

    bool IsEnd() const { return it_.eof(); }
    size_t processed() const { return rp_.processed(); }

  private:
    hammer::ReadProcessor rp_;
    ConstEdgeIteratorWrapper it_;
};

}

#endif
