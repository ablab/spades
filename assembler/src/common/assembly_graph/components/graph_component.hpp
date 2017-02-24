//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/standard_base.hpp"

namespace omnigraph {

template<class Graph>
class GraphComponent {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename std::set<VertexId>::const_iterator vertex_iterator;
    typedef typename std::set<EdgeId>::const_iterator edge_iterator;
    const Graph& graph_;
    std::set<VertexId> vertices_;
    std::set<EdgeId> edges_;
    std::set<VertexId> exits_;
    std::set<VertexId> entrances_;
    std::string name_;

    template<class VertexIt>
    void FillVertices(VertexIt begin, VertexIt end, bool add_conjugate = false) {
        for (auto it = begin; it != end; ++it) {
            vertices_.insert(*it);
            if (add_conjugate)
                vertices_.insert(graph_.conjugate(*it));
        }
    }

    template<class EdgeIt>
    void FillEdges(EdgeIt begin, EdgeIt end, bool add_conjugate = false) {
        for (auto it = begin; it != end; ++it) {
            edges_.insert(*it);
            if (add_conjugate)
                edges_.insert(graph_.conjugate(*it));
        }
    }

    void FillInducedEdges() {
        for (VertexId v : vertices_) {
            for (EdgeId e : graph_.OutgoingEdges(v)) {
                if (vertices_.count(graph_.EdgeEnd(e)) > 0) {
                    edges_.insert(e);
                }
            }
        }
    }

    void FillRelevantVertices() {
        for (EdgeId e : edges_) {
            vertices_.insert(graph_.EdgeStart(e));
            vertices_.insert(graph_.EdgeEnd(e));
        }
    }

    void FindEntrancesAndExits() {
        for (auto v : vertices_) {
            for (auto e : graph_.IncomingEdges(v)) {
                if (!contains(e)) {
                    entrances_.insert(v);
                    break;
                }
            }

            for (auto e : graph_.OutgoingEdges(v)) {
                if (!contains(e)) {
                    exits_.insert(v);
                    break;
                }
            }
        }
    }

    void Swap(GraphComponent<Graph> &that) {
        VERIFY(&this->graph_ == &that.graph_);
        std::swap(this->name_, that.name_);
        std::swap(this->vertices_, that.vertices_);
        std::swap(this->edges_, that.edges_);
        std::swap(this->exits_, that.exits_);
        std::swap(this->entrances_, that.entrances_);
    }

    template<class EdgeIt>
    void FillFromEdges(EdgeIt begin, EdgeIt end,
                       bool add_conjugate) {
        FillEdges(begin, end, add_conjugate);
        FillRelevantVertices();
        FindEntrancesAndExits();
    }

    GraphComponent<Graph> &operator=(const GraphComponent<Graph> &);
    GraphComponent(const GraphComponent<Graph> &);

public:

    template<class VertexIt>
    static GraphComponent FromVertices(const Graph &g, VertexIt begin, VertexIt end,
                                       bool add_conjugate = false, const string &name = "") {
        GraphComponent answer(g, name);
        answer.FillVertices(begin, end, add_conjugate);
        answer.FillInducedEdges();
        answer.FindEntrancesAndExits();
        return answer;
    }

    template<class EdgeIt>
    static GraphComponent FromEdges(const Graph &g, EdgeIt begin, EdgeIt end,
                                    bool add_conjugate = false, const string &name = "") {
        GraphComponent answer(g, name);
        answer.FillFromEdges(begin, end, add_conjugate);
        return answer;
    }

    template<class Container>
    static GraphComponent FromVertices(const Graph &g, const Container &c,
                                       bool add_conjugate = false, const string &name = "") {
        return FromVertices(g, c.begin(), c.end(), add_conjugate, name);
    }

    template<class Container>
    static GraphComponent FromEdges(const Graph &g, const Container &c,
                                    bool add_conjugate = false, const string &name = "") {
        return FromEdges(g, c.begin(), c.end(), add_conjugate, name);
    }

    static GraphComponent WholeGraph(const Graph &g, const string &name = "") {
        return FromVertices(g, g.begin(), g.end(), false, name);
    }

    static GraphComponent Empty(const Graph &g, const string &name = "") {
        return GraphComponent(g, name);
    }

    GraphComponent(const Graph &g, const string &name = "") :
            graph_(g), name_(name) {
    }

    //may be used for conjugate closure
    GraphComponent(const GraphComponent& component,
                   bool add_conjugate,
                   const string &name = "") : graph_(component.graph_), name_(name) {
        FillFromEdges(component.e_begin(), component.e_end(), add_conjugate);
    }

    GraphComponent(GraphComponent&& that) : graph_(that.graph_) {
        Swap(that);
    }

    GraphComponent<Graph> &operator=(GraphComponent<Graph> &&that) {
        Swap(that);
        return *this;
    }

    const Graph& g() const {
        return graph_;
    }

    string name() const {
        return name_;
    }

    size_t v_size() const {
        return vertices_.size();
    }

    size_t e_size() const {
        return edges_.size();
    }

    bool contains(EdgeId e) const {
        return edges_.count(e) > 0;
    }

    bool contains(VertexId v) const {
        return vertices_.count(v) > 0;
    }

    edge_iterator e_begin() const {
        return edges_.begin();
    }

    edge_iterator e_end() const {
        return edges_.end();
    }

    const std::set<EdgeId>& edges() const {
        return edges_;
    }

    const std::set<VertexId>& vertices() const{
        return vertices_;
    }

    vertex_iterator v_begin() const {
        return vertices_.begin();
    }

    vertex_iterator v_end() const {
        return vertices_.end();
    }

    const std::set<VertexId>& exits() const {
        return exits_;
    }

    const std::set<VertexId>& entrances() const {
        return entrances_;
    }

    bool IsBorder(VertexId v) const {
        return exits_.count(v) || entrances_.count(v);
    }

    bool empty() const {
        return v_size() == 0;
    }

};

}
