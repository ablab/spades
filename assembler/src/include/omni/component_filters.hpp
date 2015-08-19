//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
/*
 * ComponentFilters.hpp
 *
 *  Created on: Jul 10, 2013
 *      Author: anton
 */

namespace omnigraph {

template<class Element>
class AbstractFilter {
public:
    virtual ~AbstractFilter() {
    }

    virtual bool Check(const Element &element) const = 0;
};

template<class Element>
class TrueFilter : public AbstractFilter<Element> {

    bool Check(const Element &) const {
        return true;
    }
};

template<class Graph>
class GraphComponentFilter : public AbstractFilter<GraphComponent<Graph>> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph& graph_;
protected:
    GraphComponentFilter(const Graph& graph)
            : graph_(graph) {
    }

    const Graph& graph() const {
        return graph_;
    }
};

template<class Graph>
class AnyEdgeContainFilter : public GraphComponentFilter<Graph> {
    typedef GraphComponentFilter<Graph> base;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const vector<EdgeId> edges_of_interest_;
public:
    AnyEdgeContainFilter(const Graph& graph,
                         const vector<EdgeId>& edges_of_interest)
            : base(graph),
              edges_of_interest_(edges_of_interest) {

    }

    AnyEdgeContainFilter(const Graph& graph, EdgeId edge_of_interest)
            : base(graph),
              edges_of_interest_( { edge_of_interest }) {

    }

    bool ContainsEdge(const GraphComponent<Graph>& component, EdgeId e) const {
        return component.edges().find(e) != component.edges().end();
    }

    bool Check(const GraphComponent<Graph> &component) const {
        for (auto it = edges_of_interest_.begin();
                it != edges_of_interest_.end(); ++it) {
            if (ContainsEdge(component, *it)) {
                return true;
            }
        }
        return false;
    }

};

template<class Graph>
class ComponentSizeFilter : public GraphComponentFilter<Graph> {
private:
    typedef GraphComponentFilter<Graph> base;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    size_t max_length_;
    size_t min_vertex_number_;
    size_t max_vertex_number_;
public:
    ComponentSizeFilter(const Graph &graph, size_t max_length,
                        size_t min_vertex_number, size_t max_vertex_number)
            : base(graph),
              max_length_(max_length),
              min_vertex_number_(min_vertex_number),
              max_vertex_number_(max_vertex_number) {
    }

    bool Check(const GraphComponent<Graph> & component) const {
        if (component.v_size() < min_vertex_number_
                || component.v_size() > max_vertex_number_)
            return false;
        for (auto iterator = component.e_begin(); iterator != component.e_end();
                ++iterator) {
            if (this->graph().length(*iterator) <= max_length_) {
                return true;
            }
        }
        return false;
    }
};

template<class Graph>
class SmallComponentFilter : public GraphComponentFilter<Graph> {
private:
    typedef GraphComponentFilter<Graph> base;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    size_t min_vertex_number_;
public:
    SmallComponentFilter(const Graph &graph, size_t min_vertex_number)
            : base(graph), min_vertex_number_(min_vertex_number) {
    }

    bool Check(const GraphComponent<Graph> &component) const {
        return component.v_size() >= min_vertex_number_;
    }
};

}
