//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <set>
#include <unordered_set>
#include <queue>
#include <string>

#include "assembly_graph/dijkstra/dijkstra_helper.hpp"

namespace omnigraph {

template<class Graph>
class CoverageProvider {
protected:
    typedef typename Graph::EdgeId EdgeId;
    const Graph& graph_;

public:
    CoverageProvider(const Graph &graph) :
            graph_(graph) {
    }

    virtual double Coverage(EdgeId edge) const = 0;
};

template<class Graph>
class BasicCoverageProvider : public CoverageProvider<Graph> {
    typedef typename Graph::EdgeId EdgeId;

public:
    BasicCoverageProvider(const Graph &graph) :
            CoverageProvider<Graph>(graph) {

    }

    virtual double Coverage(EdgeId edge) const {
        return this->graph_.coverage(edge);
    }
};

template<class Graph>
class CustomCoverageProvider : public CoverageProvider<Graph> {
    typedef typename Graph::EdgeId EdgeId;

    const std::unordered_map<EdgeId, double>& barcode_map_;
public:
    CustomCoverageProvider(const Graph &graph, const std::unordered_map<EdgeId, double>& barcode_map) :
            CoverageProvider<Graph>(graph), barcode_map_(barcode_map)  {
    }

    virtual double Coverage(EdgeId edge) const {
        if (barcode_map_.find(edge) != barcode_map_.end()) {
            DEBUG("coverage(" << edge << ") - " << barcode_map_.at(edge));
            return barcode_map_.at(edge);
        }
        DEBUG("coverage(" << edge << ") - 0.0" );
        return 0.0;
    }
};

template<class Graph>
class GraphComponent {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename std::unordered_set<VertexId>::const_iterator vertex_iterator;
    typedef typename std::unordered_set<EdgeId>::const_iterator edge_iterator;
    const Graph& graph_;
    std::unordered_set<VertexId> vertices_;
    std::unordered_set<EdgeId> edges_;
    std::set<VertexId> exits_;
    std::set<VertexId> entrances_;
    std::string name_;
    std::shared_ptr<CoverageProvider<Graph>> coverage_provider_;

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
        std::swap(this->coverage_provider_, that.coverage_provider_);
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

    void FillInducedEdges() {
        for (VertexId v : vertices_) {
            for (EdgeId e : graph_.OutgoingEdges(v)) {
                if (vertices_.count(graph_.EdgeEnd(e)) > 0) {
                    edges_.insert(e);
                }
            }
        }
    }


    bool HasAlternative(EdgeId initial_edge) {
        VertexId v = graph_.EdgeStart(initial_edge);
        return OutgoingEdges(v).size() > 1;
    }

    double MaxAlternativeCoverage(EdgeId initial_edge) {
        VertexId v = graph_.EdgeStart(initial_edge);
        double max_coverage = std::numeric_limits<double>::min();
        for (auto e : OutgoingEdges(v)) {
            if (e == initial_edge)
                continue;
            auto edges = GetLinearPart(e);
            max_coverage = std::max(max_coverage, Coverage(edges));
        }
        return max_coverage;
    }

    std::set<VertexId> GetSources() const {
        std::set<VertexId> result;
        for (auto v : vertices_) {
            if (VertexOutDegree(v) == 0) {
                result.insert(v);
            }
        }
        return result;
    }

    std::set<VertexId> GetSinks() const {
        std::set<VertexId> result;
        for (auto v : vertices_) {
            if (VertexInDegree(v) == 0) {
                result.insert(v);
            }
        }
        return result;
    }

    size_t Length(const std::vector<EdgeId> &edges) const  {
        size_t length = 0;
        for (auto e : edges) {
            length += graph_.length(e);
        }
        return length;
    }

    std::vector<EdgeId> GetLinearPart(EdgeId e) {
        std::vector<EdgeId> start;
        std::vector<EdgeId> end;
        std::set<EdgeId> used;
        EdgeId e2 = e;
        while (IncomingEdges(graph_.EdgeStart(e2)).size() == 1 && OutgoingEdges(graph_.EdgeStart(e2)).size() == 1) {
            e2 = IncomingEdges(graph_.EdgeStart(e2)).front();
            if (used.count(e2))
                break;
            start.push_back(e2);
            used.insert(e2);
        }
        std::reverse(start.begin(), start.end());
        start.push_back(e);
        e2 = e;
        while (OutgoingEdges(graph_.EdgeEnd(e2)).size() == 1 && OutgoingEdges(graph_.EdgeEnd(e2)).size() == 1) {
            e2 = OutgoingEdges(graph_.EdgeEnd(e2)).front();
            if (used.count(e2))
                break;
            end.push_back(e2);
            used.insert(e2);
        }
        for (auto e : end) {
            start.push_back(e);
        }
        return start;
    }

    std::set<VertexId> FilterByTheSameComponent(std::set<VertexId> &intersection, VertexId v) {
        std::set<VertexId> answer;
        std::set<VertexId> reached;


        std::queue<VertexId> q;
        q.push(v);
        while (!q.empty()) {
            VertexId current = q.front();
            q.pop();
            if (reached.count(current)) {
                continue;
            }
            reached.insert(current);
            reached.insert(graph_.conjugate(current));
            for (auto e : IncidentEdges(current)) {
                if (!reached.count(graph_.EdgeStart(e))) {
                    q.push(graph_.EdgeStart(e));
                }
                if (!reached.count(graph_.EdgeEnd(e))) {
                    q.push(graph_.EdgeEnd(e));
                }
            }
        }
        //
        for (auto int_v : intersection) {
            if (!reached.count(int_v)) {
                answer.insert(int_v);
            }
        }
        return answer;
    }


public:

    template<class VertexIt>
    static GraphComponent FromVertices(const Graph &g, VertexIt begin, VertexIt end,
                                       bool add_conjugate = false, const std::string &name = "") {
        GraphComponent answer(g, name);
        answer.FillVertices(begin, end, add_conjugate);
        answer.FillInducedEdges();
        answer.FindEntrancesAndExits();
        return answer;
    }

    template<class EdgeIt>
    static GraphComponent FromEdges(const Graph &g, EdgeIt begin, EdgeIt end,
                                    bool add_conjugate = false, const std::string &name = "") {
        GraphComponent answer(g, name);
        answer.FillFromEdges(begin, end, add_conjugate);
        return answer;
    }

    template<class Container>
    static GraphComponent FromVertices(const Graph &g, const Container &c,
                                       bool add_conjugate = false, const std::string &name = "") {
        return FromVertices(g, c.begin(), c.end(), add_conjugate, name);
    }

    template<class Container>
    static GraphComponent FromEdges(const Graph &g, const Container &c,
                                    bool add_conjugate = false, const std::string &name = "") {
        return FromEdges(g, c.begin(), c.end(), add_conjugate, name);
    }

    static GraphComponent WholeGraph(const Graph &g, const std::string &name = "") {
        return FromVertices(g, g.begin(), g.end(), false, name);
    }

    static GraphComponent Empty(const Graph &g, const std::string &name = "") {
        return GraphComponent(g, name);
    }

    GraphComponent(const Graph &g, const std::string &name = "") :
            graph_(g), name_(name) {
    }

    //may be used for conjugate closure
    GraphComponent(const GraphComponent& component,
                   bool add_conjugate,
                   const std::string &name = "") : graph_(component.graph_), name_(name) {
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

    const std::string &name() const {
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

    const std::unordered_set<EdgeId>& edges() const {
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

    const std::vector<EdgeId> IncidentEdges(VertexId v) const {
        auto all_edges = graph_.IncidentEdges(v);
        std::vector<EdgeId> result;
        for (auto e : all_edges) {
            if (edges_.count(e)) {
                result.push_back(e);
            }
        }
        return result;
    }

    const std::vector<EdgeId> OutgoingEdges(VertexId v) const {
        auto all_edges = graph_.OutgoingEdges(v);
        std::vector<EdgeId> result;
        for (auto e : all_edges) {
            if (edges_.count(e)) {
                result.push_back(e);
            }
        }
        return result;
    }

    const std::vector<EdgeId> IncomingEdges(VertexId v) const {
        auto all_edges = graph_.IncomingEdges(v);
        std::vector<EdgeId> result;
        for (auto e : all_edges) {
            if (edges_.count(e)) {
                result.push_back(e);
            }
        }
        return result;
    }

    size_t VertexDegree(VertexId v) const {
        auto all_edges = graph_.IncidentEdges(v);
        DEBUG(all_edges);
        size_t result = 0;
        for (auto e : all_edges) {
            if (edges_.count(e)) {
                DEBUG(e.int_id() << " is counted");
                result++;
            }
        }
        return result;
    }

    size_t VertexInDegree(VertexId v) const {
        auto all_edges = graph_.IncomingEdges(v);
        size_t result = 0;
        for (auto e : all_edges) {
            if (edges_.count(e)) {
                result++;
            }
        }
        return result;
    }

    size_t VertexOutDegree(VertexId v) const {
        auto all_edges = graph_.OutgoingEdges(v);
        size_t result = 0;
        for (auto e : all_edges) {
            if (edges_.count(e)) {
                result++;
            }
        }
        return result;
    }

    bool IsTip(EdgeId e) {
        return VertexOutDegree(graph_.EdgeEnd(e)) == 0;
    }

    bool IsIsolated(EdgeId e) {
        return VertexInDegree(graph_.EdgeStart(e)) == 0 && VertexOutDegree(graph_.EdgeEnd(e)) == 0;

    }

    void ChangeCoverageProvider(const std::unordered_map<EdgeId, double>& barcode_map) {
        coverage_provider_ = std::make_shared<CustomCoverageProvider<Graph>>(this->graph_, barcode_map);
    }

    void RemoveIsolated() {
        std::set<EdgeId> edges_to_delete;
        std::set<VertexId> vertices_to_delete;
        size_t max_isolated_length = 100;
        for (auto e : edges_) {

            if (IsIsolated(e) && (graph_.length(e) < max_isolated_length || coverage_provider_->Coverage(e) < 0.5)) {
                edges_to_delete.insert(e);
                edges_to_delete.insert(graph_.conjugate(e));
            }
        }

        for (auto e : edges_to_delete) {
            edges_.erase(e);
        }

        for (auto v : vertices_) {
            if (VertexDegree(v) == 0) {
                vertices_to_delete.insert(v);
                vertices_to_delete.insert(graph_.conjugate(v));
            }
        }

        for (auto v : vertices_to_delete) {
            vertices_.erase(v);
        }
    }

    void ClipTips() {
        std::set<EdgeId> edges_to_delete;
        std::set<VertexId> vertices_to_delete;
        double mult_coef = 1.0;
        size_t max_tip_length = 3000;
        for (auto e : edges_) {
            if (IsTip(e)) {
                std::vector<EdgeId> tip1 = GetLinearPart(e);
                if (Length(tip1) > max_tip_length || tip1.size() > 5 || !HasAlternative(tip1.front())) {
                    continue;
                }

                double max_coverage = MaxAlternativeCoverage(tip1.front());

                if (math::ge(max_coverage, mult_coef * Coverage(tip1))) {
                    for (auto e2 : tip1) {
                        edges_to_delete.insert(e2);
                        edges_to_delete.insert(graph_.conjugate(e2));
                        vertices_to_delete.insert(graph_.EdgeEnd(e2));
                        vertices_to_delete.insert(graph_.EdgeStart(graph_.conjugate(e2)));

                    }
                }
            }
        }
        for (auto e : edges_to_delete) {
            edges_.erase(e);
        }
        for (auto v : vertices_to_delete) {
            vertices_.erase(v);
        }
    }

    void FillGaps(size_t max_distance) {
        std::set<VertexId> sources = GetSources();
        std::set<VertexId> sinks = GetSinks();

        for (auto v : sources) {
            typename omnigraph::DijkstraHelper<Graph>::BoundedDijkstra d = omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(graph_,
                                                                                                                                   max_distance,
                                                                                                                                   1000, true);
            d.Run(v);
            auto reached = d.ReachedVertices();
            std::set<VertexId> intersection;
            std::set_intersection(sinks.begin(), sinks.end(), reached.begin(), reached.end(), std::inserter(intersection, intersection.end()));
            intersection = FilterByTheSameComponent(intersection, v);
            if (intersection.size() == 1) {
                auto path = d.GetShortestPathTo(*intersection.begin());
                for (size_t i = 0; i < path.size(); ++i) {
                    edges_.insert(path[i]);
                    edges_.insert(graph_.conjugate(path[i]));
                    vertices_.insert(graph_.EdgeStart(path[i]));
                    vertices_.insert(graph_.EdgeEnd(path[i]));

                    edges_.insert(graph_.conjugate(path[i]));
                    vertices_.insert(graph_.EdgeStart(graph_.conjugate(path[i])));
                    vertices_.insert(graph_.EdgeEnd(graph_.conjugate(path[i])));
                }
                sinks = GetSinks();
            }
        }
    }



    void RemoveLowCoveredJunctions() {
        return;
        std::set<EdgeId> to_delete;
        for (auto v : vertices_) {
            if (OutgoingEdges(v).size() > 1) {
                double min_coverage = coverage_provider_->Coverage(OutgoingEdges(v).front());
                EdgeId candidate = OutgoingEdges(v).front();
                for (auto e : OutgoingEdges(v)) {
                    if (math::ls(coverage_provider_->Coverage(e), min_coverage)) {
                        candidate = e;
                        min_coverage = coverage_provider_->Coverage(e);
                    }
                }
                to_delete.insert(candidate);
            }
            if (IncomingEdges(v).size() > 1) {
                double min_coverage = coverage_provider_->Coverage(IncomingEdges(v).front());
                EdgeId candidate = IncomingEdges(v).front();
                for (auto e : IncomingEdges(v)) {
                    if (math::ls(coverage_provider_->Coverage(e), min_coverage)) {
                        candidate = e;
                        min_coverage = coverage_provider_->Coverage(e);
                    }
                }
                to_delete.insert(candidate);
            }
        }

        for (auto e : to_delete) {
            RemoveEdge(e);
        }
    }

    double Coverage(const std::vector<EdgeId> &edges) const  {
        double cov = 0.0;

        for (size_t i = 0; i < edges.size(); ++i) {
            cov += coverage_provider_->Coverage(edges[i]) * (double) graph_.length(edges[i]);
        }
        return cov / (double) Length(edges);
    }

    void RemoveEdge(EdgeId e) {
        edges_.erase(e);
        edges_.erase(graph_.conjugate(e));
        if (VertexDegree(graph_.EdgeEnd(e)) == 0) {
            vertices_.erase(graph_.EdgeEnd(e));
            vertices_.erase(graph_.conjugate(graph_.EdgeEnd(e)));
        }
        if (VertexDegree(graph_.EdgeStart(e)) == 0) {
            vertices_.erase(graph_.EdgeStart(e));
            vertices_.erase(graph_.conjugate(graph_.EdgeStart(e)));
        }

    }

};

}
