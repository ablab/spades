//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef OMNI_UTILS_HPP_
#define OMNI_UTILS_HPP_

#include "standard_base.hpp"
#include "simple_tools.hpp"
#include "adt/queue_iterator.hpp"
#include "logger/logger.hpp"
#include "simple_tools.hpp"
#include "dijkstra.hpp"
#include "xmath.h"
#include <cmath>
#include <ostream>
#include <boost/function.hpp>
#include "perfcounter.hpp"
#include <ctime>
#include "order_and_law.hpp"

namespace omnigraph {
using std::auto_ptr;
using std::vector;
using std::string;
using std::pair;
using std::set;

/**
 * ActionHandler is base listening class for graph events. All structures and information storages
 * which are meant to synchronize with graph should use this structure. In order to make handler listen
 * to graph events one should add it to graph listeners.
 * Normally structure itself extends ActionHandler and overrides several handling methods. In
 * constructor it adds itself to graph handler list and removes itself form this list in destructor.
 * All events are divided into two levels: low level events and high level events.
 * Low level events are addition/deletion of vertices/edges. These events should be triggered only after
 * high level events when all data was already transferred and graph structure is consistent.
 * High level events should be used to keep external data synchronized with graph and keep internal data
 * consistent. Now high level events are merge, glue and split. This list can be extended in near future.
 */
template<typename VertexId, typename EdgeId>
class ActionHandler : boost::noncopyable {
    const string handler_name_;
private:
    bool attached_;
 public:
    /**
     * Create action handler with given name. With this name one can find out what tipe of handler is it.
     */
    ActionHandler(const string& name)
            : handler_name_(name), attached_(false) {
    }

    virtual ~ActionHandler() {
        TRACE("~ActionHandler " << handler_name_);
    }

    /**
     * Method returns name of this handler
     */
    const string& name() const {
        return handler_name_;
    }

    /**
     * Low level event which is triggered when vertex is added to graph.
     * @param v new vertex
     */
    virtual void HandleAdd(VertexId v) {
    }

    /**
     * Low level event which is triggered when edge is added to graph.
     * @param e new edge
     */
    virtual void HandleAdd(EdgeId e) {
    }

    /**
     * Low level event which is triggered when vertex is deleted from graph.
     * @param v vertex to delete
     */
    virtual void HandleDelete(VertexId v) {
    }

    /**
     * Low level event which is triggered when edge is deleted from graph.
     * @param e edge to delete
     */
    virtual void HandleDelete(EdgeId e) {
    }

    /**
     * High level event which is triggered when merge operation is performed on graph, which is when
     * path of edges with all inner vertices having exactly one incoming and one outgoing edge is
     * replaced with a single edge. Since this is high level operation event of creation of new edge
     * and events of deletion of old edges should not have been triggered yet when this event was triggered.
     * @param old_edges path of edges to be replaced with single edge
     * @param new_edge new edge that was added to be a replacement of path
     */
    virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
    }

    /**
     * High level event which is triggered when glue operation is performed on graph, which is when
     * edge is completely replaced with other edge. This operation is widely used in bulge removal
     * when alternative path is glued to main path. Since this is high level operation event of deletion
     * of old edge should not have been triggered yet when this event was triggered.
     * @param new_edge edge glue result
     * @param edge1 edge to be glued to edge2
     * @param edge2 edge edge1 should be glued with
     */
    virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
    }

    /**
     * High level event which is triggered when split operation is performed on graph, which is when
     * edge is split into several shorter edges. Split operation is reverse to merge operation.
     * Since this is high level operation event of deletion of old edge and events of creation of new edges
     * should not have been triggered yet when this event was triggered.
     * @param old_edge edge to be split
     * @param new_edges edges which are results of split
     */
    virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
                             EdgeId new_edge_2) {
    }

    /**
     * High level event which is triggered when vertex split operation is performed on graph, which is when
     * vertex is split into several vertices, possibly doubling edges.
     * Since this is high level operation events of creation of new edges and vertex
     * should not have been triggered yet when this event was triggered.
     * @param old_vertex vertex to be split
     * @param old_2_new_edges edges which are results of split, paired with their preimage (<preimage, new_edge>)
     * @param newVertex - resulting vertex
     */
    virtual void HandleVertexSplit(
            VertexId old_vertex, VertexId new_vertex,
            const vector<pair<EdgeId, EdgeId>>& old_2_new_edges,
            const vector<double>& split_coefficients) {
    }

    /**
     * Every thread safe descendant should override this method for correct concurrent graph processing.
     */
    virtual bool IsThreadSafe() const {
        return false;
    }

    bool IsAttached() const {
        return attached_;
    }

    void Attach() {
        VERIFY(!attached_);
//        g_.AddActionHandler(this);
        attached_ = true;
    }

    void Detach() {
        VERIFY(attached_);
//        g_.RemoveActionHandler(this);
        attached_ = false;
    }
};

template<class Graph>
class GraphActionHandler : public ActionHandler<typename Graph::VertexId,
        typename Graph::EdgeId> {
    typedef ActionHandler<typename Graph::VertexId, typename Graph::EdgeId> base;

    const Graph& g_;
 protected:
    const Graph& g() const {
        return g_;
    }

 public:
    GraphActionHandler(const Graph& g, const string& name)
            : base(name),
              g_(g) {
        TRACE("Adding new action handler: " << this->name());
        g_.AddActionHandler(this);
        this->Attach();
    }

    GraphActionHandler(const GraphActionHandler<Graph> &other)
            : base(other.name()),
              g_(other.g_) {
        TRACE("Adding new action handler: " << this->name());
        g_.AddActionHandler(this);
    }

    virtual ~GraphActionHandler() {
        TRACE("Removing action handler: " << this->name());
        this->Detach();
        g_.RemoveActionHandler(this);
    }
};

/**
 * In order to support various types of graphs and make handler structure more flexible HandlerApplier
 * structure was introduced. If certain implementation of graph requires special handler triggering scheme
 * one can store certain extension of HandlerApplier in graph and trigger HandlerApplier methods instead
 * of GraphHandler methods.
 * HandlerApplier contains one method for each of graph events which define the exact way this event
 * should be triggered.
 */
template<typename VertexId, typename EdgeId>
class HandlerApplier {
    typedef ActionHandler<VertexId, EdgeId> Handler;
 public:

    virtual void
    ApplyAdd(Handler& handler, VertexId v) const = 0;

    virtual void
    ApplyAdd(Handler& handler, EdgeId e) const = 0;

    virtual void
    ApplyDelete(Handler& handler, VertexId v) const = 0;

    virtual void
    ApplyDelete(Handler& handler, EdgeId e) const = 0;

    virtual void ApplyMerge(Handler& handler, vector<EdgeId> old_edges,
                            EdgeId new_edge) const = 0;

    virtual void ApplyGlue(Handler& handler, EdgeId new_edge, EdgeId edge1,
                           EdgeId edge2) const = 0;

    virtual void ApplySplit(Handler& handler, EdgeId old_edge,
                            EdgeId new_edge_1, EdgeId new_edge2) const = 0;

    virtual void ApplyVertexSplit(
            Handler& handler, VertexId old_vertex, VertexId new_vertex,
            const vector<pair<EdgeId, EdgeId>>& old_2_new_edges,
            const vector<double>& split_coefficients) const = 0;

    virtual ~HandlerApplier() {
    }
};

/**
 * SimpleHandlerApplier is simple implementation of handler applier with no special filtering.
 */
template<class Graph>
class SimpleHandlerApplier : public HandlerApplier<typename Graph::VertexId,
        typename Graph::EdgeId> {
 public:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef ActionHandler<VertexId, EdgeId> Handler;

    virtual void ApplyAdd(Handler& handler, VertexId v) const {
        handler.HandleAdd(v);
    }

    virtual void ApplyAdd(Handler& handler, EdgeId e) const {
        handler.HandleAdd(e);
    }

    virtual void ApplyDelete(Handler& handler, VertexId v) const {
        handler.HandleDelete(v);
    }

    virtual void ApplyDelete(Handler& handler, EdgeId e) const {
        handler.HandleDelete(e);
    }

    virtual void ApplyMerge(Handler& handler, vector<EdgeId> old_edges,
                            EdgeId new_edge) const {
        handler.HandleMerge(old_edges, new_edge);
    }

    virtual void ApplyGlue(Handler& handler, EdgeId new_edge, EdgeId edge1,
                           EdgeId edge2) const {
        handler.HandleGlue(new_edge, edge1, edge2);
    }

    virtual void ApplySplit(Handler& handler, EdgeId old_edge, EdgeId new_edge1,
                            EdgeId new_edge2) const {
        handler.HandleSplit(old_edge, new_edge1, new_edge2);
    }

    virtual void ApplyVertexSplit(
            Handler& handler, VertexId old_vertex, VertexId new_vertex,
            const vector<pair<EdgeId, EdgeId>>& old_2_new_edges,
            const vector<double>& split_coefficients) const {
        handler.HandleVertexSplit(old_vertex, new_vertex, old_2_new_edges,
                                   split_coefficients);
    }

};

/**
 * PairedHandlerApplier is implementation of HandlerApplier for graph with synchronization of actions
 * performed with vertices/edges and its reverse-complement analogues. Thus while corresponding
 * method was called only once event should be triggered twice: for the parameters with which method
 * was called and for reverse-complement parameters. Also certain assertions were added for bad cases.
 */
template<class Graph>
class PairedHandlerApplier : public HandlerApplier<typename Graph::VertexId,
        typename Graph::EdgeId> {
 private:
    Graph &graph_;
 public:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef ActionHandler<VertexId, EdgeId> Handler;

    PairedHandlerApplier(Graph &graph)
            : graph_(graph) {
    }

    virtual void ApplyAdd(Handler& handler, VertexId v) const {
        VertexId rcv = graph_.conjugate(v);
        handler.HandleAdd(v);
        if (v != rcv) {
            handler.HandleAdd(rcv);
        }
    }

    virtual void ApplyAdd(Handler& handler, EdgeId e) const {
        EdgeId rce = graph_.conjugate(e);
        handler.HandleAdd(e);
        if (e != rce) {
            handler.HandleAdd(rce);
        }
    }

    virtual void ApplyDelete(Handler& handler, VertexId v) const {
        VertexId rcv = graph_.conjugate(v);
        handler.HandleDelete(v);
        if (v != rcv) {
            handler.HandleDelete(rcv);
        }
    }

    virtual void ApplyDelete(Handler& handler, EdgeId e) const {
        EdgeId rce = graph_.conjugate(e);
        handler.HandleDelete(e);
        if (e != rce) {
            handler.HandleDelete(rce);
        }
    }

    virtual void ApplyMerge(Handler& handler, vector<EdgeId> old_edges,
                            EdgeId new_edge) const {
        EdgeId rce = graph_.conjugate(new_edge);
        handler.HandleMerge(old_edges, new_edge);
        if (new_edge != rce) {
            vector<EdgeId> rc_old_edges;
            for (int i = old_edges.size() - 1; i >= 0; i--) {
                rc_old_edges.push_back(graph_.conjugate(old_edges[i]));
            }
            handler.HandleMerge(rc_old_edges, rce);
        }
    }

    virtual void ApplyGlue(Handler& handler, EdgeId new_edge, EdgeId edge1,
                           EdgeId edge2) const {
        EdgeId rc_edge1 = graph_.conjugate(edge1);
        EdgeId rc_edge2 = graph_.conjugate(edge2);
        VERIFY(edge1 != edge2);
        VERIFY(edge2 != rc_edge2);
        handler.HandleGlue(new_edge, edge1, edge2);
        if (edge1 != rc_edge1) {
            handler.HandleGlue(graph_.conjugate(new_edge), rc_edge1, rc_edge2);
        }
    }

    virtual void ApplySplit(Handler& handler, EdgeId old_edge,
                            EdgeId new_edge_1, EdgeId new_edge2) const {
        EdgeId rce = graph_.conjugate(old_edge);
        VERIFY(old_edge != rce);
        handler.HandleSplit(old_edge, new_edge_1, new_edge2);
        if (old_edge != rce) {
            handler.HandleSplit(rce, graph_.conjugate(new_edge2),
                                 graph_.conjugate(new_edge_1));
        }
    }

    virtual void ApplyVertexSplit(
            Handler& handler, VertexId old_vertex, VertexId new_vertex,
            const vector<pair<EdgeId, EdgeId>>& old_2_new_edges,
            const vector<double>& split_coefficients) const {

        //todo add checks that there are no edges and conjugate

        VERIFY(old_2_new_edges.size() == split_coefficients.size());
        handler.HandleVertexSplit(old_vertex, new_vertex, old_2_new_edges,
                                   split_coefficients);

        vector<pair<EdgeId, EdgeId>> rc_old_2_new_edges;
        FOREACH (auto old_2_new_edge, old_2_new_edges) {
            //not reversing order
            rc_old_2_new_edges.push_back(
                    make_pair(graph_.conjugate(old_2_new_edge.first),
                              graph_.conjugate(old_2_new_edge.second)));
        }

        VERIFY(rc_old_2_new_edges.size() == split_coefficients.size());
        handler.HandleVertexSplit(graph_.conjugate(old_vertex),
                                   graph_.conjugate(new_vertex),
                                   rc_old_2_new_edges, split_coefficients);
    }

 private:
    DECL_LOGGER("PairedHandlerApplier")
};

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

    virtual ~SmartIterator() {
    }

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
template<class Graph, typename ElementId, typename Comparator = std::less<
        ElementId>>
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

    virtual ~SmartSetIterator() {
    }
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

  virtual ~ConditionedSmartSetIterator() {
  }

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

  DECL_LOGGER("ConditionedSmartSetIterator")
  ;
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

    virtual ~SmartVertexIterator() {
    }

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

/**
 * This class is a representation of how certain sequence is mapped to genome. Needs further adjustment.
 */
template<typename ElementId>
class Path {
    vector<ElementId> sequence_;
    int start_pos_;
    int end_pos_;
 public:
    typedef typename vector<ElementId>::const_iterator iterator;

    Path(const vector<ElementId>& sequence, size_t start_pos, size_t end_pos)
            : sequence_(sequence),
              start_pos_(start_pos),
              end_pos_(end_pos) {
    }

    Path()
            : sequence_(),
              start_pos_(-1),
              end_pos_(-1) {
    }

    size_t start_pos() const {
        return start_pos_;
    }

    size_t end_pos() const {
        return end_pos_;
    }

    size_t size() const {
        return sequence_.size();
    }

    const vector<ElementId>& sequence() const {
        return sequence_;
    }

    ElementId operator[](size_t index) const {
        return sequence_[index];
    }

    iterator begin() const {
        return sequence_.begin();
    }

    iterator end() const {
        return sequence_.end();
    }

};

struct Range {
    //inclusive
    size_t start_pos;
    //exclusive
    size_t end_pos;

    size_t size() const {
        VERIFY(end_pos >= start_pos);
        return end_pos - start_pos;
    }

    void shift(int shift) {
        VERIFY(shift > 0 || size_t(-shift) <= start_pos);
        start_pos += shift;
        end_pos += shift;
    }

    Range(size_t start_pos, size_t end_pos)
            : start_pos(start_pos),
              end_pos(end_pos) {
        VERIFY(end_pos >= start_pos);
    }
};

std::ostream& operator<<(std::ostream& os, const Range& range) {
    os << "[" << range.start_pos << ", " << range.end_pos << "]";
    return os;
}

struct MappingRange {
    Range initial_range;
    Range mapped_range;

    MappingRange(Range initial_range, Range mapped_range)
            : initial_range(initial_range),
              mapped_range(mapped_range) {
    }
};

std::ostream& operator<<(std::ostream& os, const MappingRange& map_range) {
    os << map_range.initial_range << " --> " << map_range.mapped_range;
    return os;
}

template<typename ElementId>
class MappingPath {
 public:

    MappingPath() {
    }

    MappingPath(const vector<ElementId>& edges,
                const vector<MappingRange> range_mappings)
            : edges_(edges),
              range_mappings_(range_mappings) {
    }

    size_t size() const {
        return edges_.size();
    }

    pair<const ElementId, const MappingRange> operator[](size_t idx) const {
        return make_pair(edges_[idx], range_mappings_[idx]);
    }

    pair<const ElementId, const MappingRange> front() const {
        return make_pair(edges_.front(), range_mappings_.front());
    }

    pair<const ElementId, const MappingRange> back() const {
        return make_pair(edges_.back(), range_mappings_.back());
    }

    size_t start_pos() const {
        return range_mappings_.front().mapped_range.start_pos;
    }

    size_t end_pos() const {
        return range_mappings_.back().mapped_range.end_pos;
    }

    Path<ElementId> simple_path() const {
        if (edges_.size() != 0)
            return Path<ElementId>(
                    edges_,
                    range_mappings_[0].mapped_range.start_pos,
                    range_mappings_[range_mappings_.size() - 1].mapped_range
                            .end_pos);
        else
            return Path<ElementId>();
    }

 private:
    vector<ElementId> edges_;
    vector<MappingRange> range_mappings_;
};

template<class Graph>
class BackwardBoundedDijkstra : public BackwardDijkstra<Graph> {
 private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef BackwardDijkstra<Graph> base;
    const size_t bound_;

 public:
    BackwardBoundedDijkstra(const Graph &g, size_t bound)
            : base(g),
              bound_(bound) {
    }

    virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
        return distance <= bound_;
    }

};

template<class Graph>
class BackwardReliableBoundedDijkstra : public BackwardDijkstra<Graph> {

    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef BackwardDijkstra<Graph> base;

 public:
    BackwardReliableBoundedDijkstra(const Graph &g, size_t bound,
                                    size_t max_vertex_number)
            : base(g),
              bound_(bound),
              max_vertex_number_(max_vertex_number),
              vertices_number_(0),
              vertex_limit_exceeded_(false) {
    }

    virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
        ++vertices_number_;

        if (vertices_number_ > max_vertex_number_)
            vertex_limit_exceeded_ = true;

        return vertices_number_ < max_vertex_number_ && distance <= bound_;
    }

    bool VertexLimitExceeded() const {
        return vertex_limit_exceeded_;
    }

 private:
    const size_t bound_;
    const size_t max_vertex_number_;
    size_t vertices_number_;
    bool vertex_limit_exceeded_;
};

template<class Graph>
class ReliableBoundedDijkstra : public Dijkstra<Graph> {

    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef Dijkstra<Graph> base;

 public:
    ReliableBoundedDijkstra(const Graph& g, size_t bound,
                            size_t max_vertex_number)
            : base(g),
              bound_(bound),
              max_vertex_number_(max_vertex_number),
              vertices_number_(0),
              vertex_limit_exceeded_(false) {
    }

    virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
        ++vertices_number_;

        if (vertices_number_ > max_vertex_number_)
            vertex_limit_exceeded_ = true;

        return (vertices_number_ < max_vertex_number_) && (distance <= bound_);
    }

    bool VertexLimitExceeded() const {
        return vertex_limit_exceeded_;
    }

 private:
    const size_t bound_;
    const size_t max_vertex_number_;
    size_t vertices_number_;
    bool vertex_limit_exceeded_;
};

template<class Graph>
struct CoverageComparator {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph& graph_;
 public:
    CoverageComparator(const Graph &graph)
            : graph_(graph) {
    }

    /**
     * Standard comparator function as used in collections.
     */
    bool operator()(EdgeId edge1, EdgeId edge2) const {
        if (math::eq(graph_.coverage(edge1), graph_.coverage(edge2))) {
            return edge1 < edge2;
        }
        return math::ls(graph_.coverage(edge1), graph_.coverage(edge2));
    }
};

/**
 * This class defines which edge is more likely to be tip. In this case we just assume shorter edges
 * are more likely tips then longer ones.
 */
template<class Graph>
struct LengthComparator {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph& graph_;
 public:
    /**
     * TipComparator should never be created with default constructor but it is necessary on order for
     * code to compile.
     */
    //  TipComparator() {
    //    VERIFY(false);
    //  }
    /**
     * Construct TipComparator for given graph
     * @param graph graph for which comparator is created
     */
    LengthComparator(const Graph &graph)
            : graph_(graph) {
    }

    /**
     * Standard comparator function as used in collections.
     */
    bool operator()(EdgeId edge1, EdgeId edge2) const {
        if (graph_.length(edge1) == graph_.length(edge2)) {
            return edge1 < edge2;
        }
        return graph_.length(edge1) < graph_.length(edge2);
    }
};

template<class Graph>
size_t CummulativeLength(const Graph& g,
                         const vector<typename Graph::EdgeId>& path) {
    size_t s = 0;
    for (auto it = path.begin(); it != path.end(); ++it) {
        s += g.length(*it);
    }
    return s;
}

template<class Graph>

class AbstractDirection {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& graph_;

 protected:
    const Graph &graph() const {
        return graph_;
    }

 public:
    AbstractDirection(const Graph& graph)
            : graph_(graph) {
    }

    virtual ~AbstractDirection() {
    }

    virtual const vector<EdgeId> OutgoingEdges(VertexId v) const = 0;

    virtual const vector<EdgeId> IncomingEdges(VertexId v) const = 0;

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

    virtual const vector<EdgeId> OutgoingEdges(VertexId v) const {
        return this->graph().OutgoingEdges(v);
    }

    virtual const vector<EdgeId> IncomingEdges(VertexId v) const {
        return this->graph().IncomingEdges(v);
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

    virtual const vector<EdgeId> OutgoingEdges(VertexId v) const {
        return this->graph().IncomingEdges(v);
    }

    virtual const vector<EdgeId> IncomingEdges(VertexId v) const {
        return this->graph().OutgoingEdges(v);
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

template<class Graph>
class UniquePathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& graph_;
 public:

    //todo use length bound if needed
    UniquePathFinder(const Graph& graph, size_t length_bound =
                             std::numeric_limits<size_t>::max())
            : graph_(graph) {

    }

    const vector<EdgeId> operator()(
            EdgeId e, const AbstractDirection<Graph> &direction) const {
        vector<EdgeId> answer;
        EdgeId curr = e;
        answer.push_back(curr);
        set<EdgeId> was;
        while (direction.CheckUniqueOutgoingEdge(direction.EdgeEnd(curr))) {
            curr = direction.GetUniqueOutgoingEdge(direction.EdgeEnd(curr));
            if (was.count(curr) > 0)
                break;
            was.insert(curr);
            answer.push_back(curr);
        }
        return answer;
    }

    const vector<EdgeId> UniquePathForward(EdgeId e) const {
        return this->operator()(e, ForwardDirection<Graph>(graph_));
    }

    const vector<EdgeId> UniquePathBackward(EdgeId e) const {
        return this->operator()(e, BackwardDirection<Graph>(graph_));
    }

//	const vector<EdgeId> UniquePathBackward(EdgeId e) const {
//		TRACE("UniquePathBackward from " << graph_.str(e));
//		vector<EdgeId> answer;
//		EdgeId curr = e;
//		answer.push_back(curr);
//		set<EdgeId> was;
//		while (graph_.CheckUniqueIncomingEdge(graph_.EdgeStart(curr))) {
//			TRACE("current " << curr);
//			curr = graph_.GetUniqueIncomingEdge(graph_.EdgeStart(curr));
//			if (was.count(curr) > 0)
//				break;
//			was.insert(curr);
//			answer.push_back(curr);
//		}
//		TRACE("UniquePathBackward from " << graph_.str(e) << " finished");
//		return vector<EdgeId>(answer.rbegin(), answer.rend());
//	}
};

template<class Graph>
class TrivialPathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

 public:

    TrivialPathFinder(const Graph&, size_t stub = 0) {

    }

    const vector<EdgeId> operator()(
            EdgeId e, const AbstractDirection<Graph> &direction) const {
        return {e};
    }

};

template<class Graph>
class PlausiblePathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    //todo remove graph_ field???
    const Graph& graph_;
    const size_t length_bound_;

    class DFS {
     private:
        const Graph &graph_;
        const AbstractDirection<Graph> &direction_;
        const size_t length_bound_;

        pair<size_t, EdgeId> find(EdgeId edge, size_t length) {
            length += graph_.length(edge);
            VertexId cross = direction_.EdgeEnd(edge);
            auto result = make_pair(length, edge);
            if (length < length_bound_
                    && direction_.CheckUniqueIncomingEdge(cross)) {
                vector<EdgeId> outgoing = direction_.OutgoingEdges(cross);
                for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
                    auto candidate = find(*it, length);
                    if (candidate.first > result.first)
                        result = candidate;
                }
            }
            return result;
        }

        vector<EdgeId> RestoreAnswer(EdgeId start, EdgeId end) {
            vector<EdgeId> result;
            while (end != start) {
                result.push_back(end);
                end = direction_.GetUniqueIncomingEdge(
                        direction_.EdgeStart(end));
            }
            result.push_back(start);
            return vector<EdgeId>(result.rbegin(), result.rend());
        }

     public:
        DFS(const Graph &graph, const AbstractDirection<Graph> &direction,
            size_t length_bound)
                : graph_(graph),
                  direction_(direction),
                  length_bound_(length_bound) {
        }

        vector<EdgeId> find(EdgeId edge) {
            vector<EdgeId> result = RestoreAnswer(edge, find(edge, 0).second);
            return result;
        }
    };

 public:
    PlausiblePathFinder(const Graph& graph, size_t length_bound)
            : graph_(graph),
              length_bound_(length_bound) {
    }

    const vector<EdgeId> operator()(
            EdgeId e, const AbstractDirection<Graph> &direction) const {
        vector<EdgeId> answer;
        return DFS(graph_, direction, length_bound_).find(e);
    }

};

template<class Graph>
class MultiplicityCounter {
 private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    size_t uniqueness_length_;
    size_t max_depth_;

    bool search(VertexId a, VertexId start, EdgeId e, size_t depth,
                set<VertexId> &was, pair<size_t, size_t> &result) const {
        if (depth > max_depth_)
            return false;
        if (was.count(a) == 1)
            return true;
        was.insert(a);
        if (graph_.OutgoingEdgeCount(a) == 0
                || graph_.IncomingEdgeCount(a) == 0)
            return false;
        for (auto I = graph_.out_begin(a), E = graph_.out_end(a); I != E; ++I) {
            if (*I == e) {
                if (a != start) {
                    return false;
                }
            } else {
                if (graph_.length(*I) >= uniqueness_length_) {
                    result.second++;
                } else {
                    if (!search(graph_.EdgeEnd(*I), start, e,
                                depth + 1 /*graph_.length(*it)*/, was, result))
                        return false;
                }
            }
        }
        vector<EdgeId> in = graph_.IncomingEdges(a);
        for (auto it = in.begin(); it != in.end(); ++it) {
            if (*it == e) {
                if (a != start) {
                    return false;
                }
            } else {
                if (graph_.length(*it) >= uniqueness_length_) {
                    result.first++;
                } else {
                    if (!search(graph_.EdgeStart(*it), start, e,
                                depth + 1 /*graph_.length(*it)*/, was, result))
                        return false;
                }
            }
        }
        return true;
    }

 public:
    MultiplicityCounter(const Graph &graph, size_t uniqueness_length,
                        size_t max_depth)
            : graph_(graph),
              uniqueness_length_(uniqueness_length),
              max_depth_(max_depth) {
    }

    size_t count(EdgeId e, VertexId start) const {
        pair<size_t, size_t> result;
        set<VertexId> was;
        bool valid = search(start, start, e, 0, was, result);
        if (!valid) {
            return (size_t) (-1);
        }
        if (graph_.EdgeStart(e) == start) {
            if (result.first < result.second) {
                return (size_t) (-1);
            }
            return result.first - result.second;
        } else {
            if (result.first > result.second) {
                return (size_t) (-1);
            }
            return -result.first + result.second;
        }
    }
};

template<class Graph>
class DominatedSetFinder {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    VertexId start_vertex_;
    size_t max_length_;
    size_t max_count_;

    size_t cnt_;
    map<VertexId, Range> dominated_;

    bool CheckCanBeProcessed(VertexId v) const {
        DEBUG( "Check if vertex " << g_.str(v) << " is dominated close neighbour");
        FOREACH (EdgeId e, g_.IncomingEdges(v)) {
            if (dominated_.count(g_.EdgeStart(e)) == 0) {
                DEBUG( "Blocked by external vertex " << g_.int_id(g_.EdgeStart(e)) << " that starts edge " << g_.int_id(e));
                DEBUG("Check fail");
                return false;
            }
        }
        DEBUG("Check ok");
        return true;
    }

    void UpdateCanBeProcessed(VertexId v,
                              std::queue<VertexId>& can_be_processed) const {
        DEBUG("Updating can be processed")
        FOREACH (EdgeId e, g_.OutgoingEdges(v)) {
            DEBUG("Considering edge " << ToString(e));
            VertexId neighbour_v = g_.EdgeEnd(e);
            if (CheckCanBeProcessed(neighbour_v)) {
                can_be_processed.push(neighbour_v);
            }
        }
    }

    Range NeighbourDistanceRange(VertexId v, bool dominated_only = true) const {
        DEBUG("Counting distance range for vertex " << g_.str(v));
        size_t min = numeric_limits<size_t>::max();
        size_t max = 0;
        VERIFY(g_.IncomingEdgeCount(v) > 0);
        VERIFY(!dominated_only || CheckCanBeProcessed(v));
        FOREACH (EdgeId e, g_.IncomingEdges(v)) {
            //in case of dominated_only == false
            if (dominated_.count(g_.EdgeStart(e)) == 0)
                continue;
            Range range = dominated_.find(g_.EdgeStart(e))->second;
            range.shift(g_.length(e));
            DEBUG("Edge " << g_.str(e) << " provide distance range " << range);
            if (range.start_pos < min)
                min = range.start_pos;
            if (range.end_pos > max)
                max = range.end_pos;
        }
        VERIFY((max > 0) && (min < numeric_limits<size_t>::max()) && (min <= max));
        Range answer(min, max);
        DEBUG("Range " << answer);
        return answer;
    }

    bool CheckNoEdgeToStart(VertexId v) {
        FOREACH (EdgeId e, g_.OutgoingEdges(v)) {
            if (g_.EdgeEnd(e) == start_vertex_) {
                return false;
            }
        }
        return true;
    }

 public:
    DominatedSetFinder(const Graph& g, VertexId v, size_t max_length = -1,
                       size_t max_count = -1)
            : g_(g),
              start_vertex_(v),
              max_length_(max_length),
              max_count_(max_count),
              cnt_(0) {

    }

    //true if no thresholds exceeded
    bool FillDominated() {
        DEBUG("Adding starting vertex " << g_.str(start_vertex_) << " to dominated set");
        dominated_.insert(make_pair(start_vertex_, Range(0, 0)));
        cnt_++;
        std::queue<VertexId> can_be_processed;
        UpdateCanBeProcessed(start_vertex_, can_be_processed);
        while (!can_be_processed.empty()) {
            if (++cnt_ > max_count_) {
                return false;
            }
            VertexId v = can_be_processed.front();
            can_be_processed.pop();
            Range r = NeighbourDistanceRange(v);
            if (r.start_pos > max_length_) {
                return false;
            }
            //Currently dominated vertices cannot have edge to start vertex
            if (CheckNoEdgeToStart(v)) {
                DEBUG("Adding vertex " << g_.str(v) << " to dominated set");
                dominated_.insert(make_pair(v, r));
                UpdateCanBeProcessed(v, can_be_processed);
            }
        }
        return true;
    }

    const map<VertexId, Range>& dominated() const {
        return dominated_;
    }

    //little meaning if FillDominated returned false
    const map<VertexId, Range> CountBorder() const {
        map<VertexId, Range> border;
        FOREACH(VertexId v, key_set(border)) {
            FOREACH(EdgeId e, g_.OutgoingEdges(v)) {
                VertexId e_end = g_.EdgeEnd(e);
                if (dominated_.count(e_end) == 0) {
                    border[e_end] = NeighbourDistanceRange(e_end, false);
                }
            }
        }
        return border;
    }

};

inline size_t PairInfoPathLengthUpperBound(size_t k, size_t insert_size,
                                           double delta) {
    double answer = 0. + insert_size + delta - k - 2;
    VERIFY(math::gr(answer, 0.));
    return std::floor(answer);
}

inline size_t PairInfoPathLengthLowerBound(size_t k, size_t l1, size_t l2,
                                           int gap, double delta) {
    double answer = 0. + gap + k + 2 - l1 - l2 - delta;
    return math::gr(answer, 0.) ? std::floor(answer) : 0;
}

}
#endif /* OMNI_UTILS_HPP_ */
