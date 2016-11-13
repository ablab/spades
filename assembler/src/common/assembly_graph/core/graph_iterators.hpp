//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/adt/queue_iterator.hpp"
#include "func/pred.hpp"
#include "action_handlers.hpp"
#include "utils/simple_tools.hpp"
#include <boost/iterator/iterator_facade.hpp>

namespace omnigraph {

/**
 * SmartIterator is able to iterate through collection content of which can be changed in process of
 * iteration. And as GraphActionHandler SmartIterator can change collection contents with respect to the
 * way graph is changed. Also one can define order of iteration by specifying Comparator.
 */
template<class Graph, typename ElementId, typename Comparator = std::less<ElementId>>
class SmartIterator : public GraphActionHandler<Graph> {
    typedef GraphActionHandler<Graph> base;
    DynamicQueueIterator<ElementId, Comparator> inner_it_;
    bool add_new_;
    bool canonical_only_;
    //todo think of checking it in HandleAdd
    func::TypedPredicate<ElementId> add_condition_;

protected:

    void push(const ElementId& el) {
        if ((!canonical_only_ || el <= this->g().conjugate(el)) &&
            add_condition_(el)) {
            inner_it_.push(el);
        }
    }

    template<typename InputIterator>
    void insert(InputIterator begin, InputIterator end) {
        for (auto it = begin; it != end; ++it) {
            push(*it);
        }
    }

    void erase(const ElementId& el) {
        if (!canonical_only_ || el <= this->g().conjugate(el)) {
            inner_it_.erase(el);
        }
    }

    void clear() {
        inner_it_.clear();
    }

    SmartIterator(const Graph &g, const std::string &name, bool add_new,
                  const Comparator& comparator, bool canonical_only,
                  func::TypedPredicate<ElementId> add_condition = func::AlwaysTrue<ElementId>())
            : base(g, name),
              inner_it_(comparator),
              add_new_(add_new),
              canonical_only_(canonical_only),
              add_condition_(add_condition) {
    }

public:

    bool canonical_only() const {
        return canonical_only_;
    }

    bool IsEnd() const {
        return inner_it_.IsEnd();
    }

    size_t size() const {
        return inner_it_.size();
    }

    ElementId operator*() {
        return *inner_it_;
    }

    void operator++() {
        ++inner_it_;
    }

    void HandleAdd(ElementId v) override {
        if (add_new_)
            push(v);
    }

    void HandleDelete(ElementId v) override {
        erase(v);
    }

    //use carefully!
    void ReleaseCurrent() {
        inner_it_.ReleaseCurrent();
    }

};

/**
 * SmartIterator is abstract class which acts both as QueueIterator and GraphActionHandler. As QueueIterator
 * SmartIterator is able to iterate through collection content of which can be changed in process of
 * iteration. And as GraphActionHandler SmartIterator can change collection contents with respect to the
 * way graph is changed. Also one can define order of iteration by specifying Comparator.
 */
template<class Graph, typename ElementId,
         typename Comparator = std::less<ElementId>>
class SmartSetIterator : public SmartIterator<Graph, ElementId, Comparator> {
    typedef SmartIterator<Graph, ElementId, Comparator> base;

public:
    SmartSetIterator(const Graph &g,
                     bool add_new = false,
                     const Comparator& comparator = Comparator(),
                     bool canonical_only = false,
                     func::TypedPredicate<ElementId> add_condition = func::AlwaysTrue<ElementId>())
            : base(g, "SmartSet " + ToString(this), add_new, comparator, canonical_only, add_condition) {
    }

    template<class Iterator>
    SmartSetIterator(const Graph &g, Iterator begin, Iterator end,
                     bool add_new = false,
                     const Comparator& comparator = Comparator(),
                     bool canonical_only = false,
                     func::TypedPredicate<ElementId> add_condition = func::AlwaysTrue<ElementId>())
            : SmartSetIterator(g, add_new, comparator, canonical_only, add_condition) {
        insert(begin, end);
    }

    template<typename InputIterator>
    void insert(InputIterator begin, InputIterator end) {
        base::insert(begin, end);
    }

    void push(const ElementId& el) {
        base::push(el);
    }

    void clear() {
        base::clear();
    }
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

    static size_t get_id() {
        static size_t id = 0;
        return id++;
    }

  public:
    SmartVertexIterator(const Graph &g, const Comparator& comparator =
                        Comparator(), bool canonical_only = false)
            : SmartIterator<Graph, VertexId, Comparator>(
                g, "SmartVertexIterator " + ToString(get_id()), true,
                comparator, canonical_only) {
        this->insert(g.begin(), g.end());
    }

};

//todo return verifies when they can be switched off
template<class Graph>
class GraphEdgeIterator : public boost::iterator_facade<GraphEdgeIterator<Graph>
                                                    , typename Graph::EdgeId, boost::forward_traversal_tag
                                                    , typename Graph::EdgeId> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexIt const_vertex_iterator;
    typedef typename Graph::edge_const_iterator const_edge_iterator;

    const Graph& g_;
    const_vertex_iterator v_it_;
    const_edge_iterator e_it_;
    bool canonical_only_;

public:

    GraphEdgeIterator(const Graph& g, const_vertex_iterator v_it, bool canonical_only = false)
            : g_(g),
              v_it_(v_it),
              canonical_only_(canonical_only) {
        if (v_it_ != g_.end()) {
            e_it_ = g_.out_begin(*v_it_);
            Skip();
        }
    }

private:
    
    bool Canonical(EdgeId e) const {
        return e <= g_.conjugate(e);
    }

    friend class boost::iterator_core_access;

    void Skip() {
        //VERIFY(v_it_ != g_.end());
        while (true) {
            if (e_it_ == g_.out_end(*v_it_)) {
                v_it_++;
                if (v_it_ == g_.end())
                    return;
                e_it_ = g_.out_begin(*v_it_);
            } else {
                if (!canonical_only_ || Canonical(*e_it_))
                    return;
                else
                    e_it_++;
            }
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
        if (other.canonical_only_ != canonical_only_)
            return false;
        return true;
    }

    EdgeId dereference() const {
        //VERIFY(v_it_ != g_.end());
        return *e_it_;
    }

};

template<class Graph>
class ConstEdgeIterator {
    typedef typename Graph::EdgeId EdgeId;
    GraphEdgeIterator<Graph> begin_, end_;

  public:
    ConstEdgeIterator(const Graph &g, bool canonical_only = false)
            : begin_(g, g.begin(), canonical_only), end_(g, g.end(), canonical_only) {
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

/**
 * SmartEdgeIterator iterates through edges of graph. It listens to AddEdge/DeleteEdge graph events
 * and correspondingly edits the set of edges to iterate through. Note: high level event handlers are
 * triggered before low level event handlers like HandleAdd/HandleDelete. Thus if Comparator uses certain
 * structure which is also updated with handlers make sure that all information is updated in high level
 * event handlers.
 */
template<class Graph, typename Comparator = std::less<typename Graph::EdgeId> >
class SmartEdgeIterator : public SmartIterator<Graph, typename Graph::EdgeId, Comparator> {
    typedef GraphEdgeIterator<Graph> EdgeIt;
  public:
    typedef typename Graph::EdgeId EdgeId;

    static size_t get_id() {
        static size_t id = 0;
        return id++;
    }

  public:
    SmartEdgeIterator(const Graph &g, Comparator comparator = Comparator(),
                      bool canonical_only = false)
            : SmartIterator<Graph, EdgeId, Comparator>(
                g, "SmartEdgeIterator " + ToString(get_id()), true,
                comparator, canonical_only) {
        this->insert(EdgeIt(g, g.begin()), EdgeIt(g, g.end()));

//        for (auto it = graph.begin(); it != graph.end(); ++it) {
//            //todo: this solution doesn't work with parallel simplification
//            this->insert(graph.out_begin(*it), graph.out_end(*it));
//            //this does
//            //auto out = graph.OutgoingEdges(*it);
//            //this->base::insert(out.begin(), out.end());
//        }
    }
};

//todo move out
template<class Graph, class ElementId>
class IterationHelper {
};

template<class Graph>
class IterationHelper<Graph, typename Graph::VertexId> {
    const Graph& g_;
public:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::VertexIt const_vertex_iterator;

    IterationHelper(const Graph& g)
            : g_(g) {
    }

    const_vertex_iterator begin() const {
        return g_.begin();
    }

    const_vertex_iterator end() const {
        return g_.end();
    }

    std::vector<const_vertex_iterator> Chunks(size_t chunk_cnt) const {
        VERIFY(chunk_cnt > 0);
        if (chunk_cnt == 1) {
            return {begin(), end()};
        }

        //trying to split vertices into equal chunks, leftovers put into first chunk
        vector<const_vertex_iterator> answer;
        size_t vertex_cnt = g_.size();
        size_t chunk_size = vertex_cnt / chunk_cnt;
        auto it = g_.begin();
        answer.push_back(it);
        for (size_t i = 0; i + chunk_cnt * chunk_size < vertex_cnt; ++i) {
            it++;
        }
        if (chunk_size > 0) {
            size_t i = 0;
            do {
                ++it;
                if (++i % chunk_size == 0)
                    answer.push_back(it);
            } while (it != g_.end());

            VERIFY(i == chunk_cnt * chunk_size);
        } else {
            VERIFY(it == g_.end());
            answer.push_back(it);
        }
        VERIFY(answer.back() == g_.end());
        return answer;
    }

};

//todo move out
template<class Graph>
class IterationHelper<Graph, typename Graph::EdgeId> {
    typedef typename Graph::VertexId VertexId;

    const Graph& g_;
public:
    typedef typename Graph::EdgeId EdgeId;
    typedef GraphEdgeIterator<Graph> const_edge_iterator;

    IterationHelper(const Graph& g)
            : g_(g) {
    }

    const_edge_iterator begin() const {
        return const_edge_iterator(g_, g_.begin());
    }

    const_edge_iterator end() const {
        return const_edge_iterator(g_, g_.end());
    }

    std::vector<omnigraph::GraphEdgeIterator<Graph>> Chunks(size_t chunk_cnt) const {
        if (chunk_cnt == 1) {
            return {begin(), end()};
        }

        vector<omnigraph::GraphEdgeIterator<Graph>> answer;

        for (auto v_it : IterationHelper<Graph, VertexId>(g_).Chunks(chunk_cnt)) {
            answer.push_back(omnigraph::GraphEdgeIterator<Graph>(g_, v_it));
        }
        return answer;
    }
};

}
