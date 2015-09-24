//
// Created by andrey on 17.09.15.
//
#pragma once

#include "logger/logger.hpp"
#include "standard_base.hpp"
#include "debruijn_graph.hpp"
#include "paired_library.hpp"


namespace scaffold_graph {

//do NOT add using debruijn_graph due to EdgeId conflict

using path_extend::PairedInfoLibrary;


class ConnectionCondition {
public:
    virtual set<debruijn_graph::EdgeId> ConnectedWith(debruijn_graph::EdgeId e) const = 0;

    virtual double GetWeight(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const = 0;

    virtual size_t GetLibIndex() const = 0;

    virtual ~ConnectionCondition() {

    }
};

class PairedLibConnectionCondition: public ConnectionCondition {
private:

    const debruijn_graph::Graph& graph_;

    shared_ptr<PairedInfoLibrary> lib_;

    size_t lib_index_;

    size_t min_read_count_;


public:
    PairedLibConnectionCondition(const debruijn_graph::Graph& graph,
                                 shared_ptr<PairedInfoLibrary> lib,
                                 size_t lib_index,
                                 size_t min_read_count):
            graph_(graph),
            lib_(lib),
            lib_index_(lib_index),
            min_read_count_(min_read_count) {

    }

    size_t GetLibIndex() const {
        return lib_index_;
    }

    set<debruijn_graph::EdgeId> ConnectedWith(debruijn_graph::EdgeId e) const {
        set<debruijn_graph::EdgeId> all_edges;
        lib_->FindJumpEdges(e, all_edges);

        set<debruijn_graph::EdgeId> result;
        for (auto edge : all_edges) {
            if (edge != e && edge != graph_.conjugate(e) && math::ge(GetWeight(e, edge), (double) min_read_count_)) {
                result.insert(edge);
            }
        }
        return result;
    }

    double GetWeight(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const {
        return lib_->CountPairedInfo(e1, e2);
    }
};

class AssemblyGraphConnectionCondition: public ConnectionCondition {
private:
    const debruijn_graph::Graph& g_;

    size_t max_connection_length_;

public:
    AssemblyGraphConnectionCondition(const debruijn_graph::Graph& g, size_t max_connection_length):
            g_(g),
            max_connection_length_(max_connection_length) {
    }

    set<debruijn_graph::EdgeId> ConnectedWith(debruijn_graph::EdgeId e) const {
        set<debruijn_graph::EdgeId> result;

        for (auto connected: g_.OutgoingEdges(g_.EdgeEnd(e))) {
            result.insert(connected);
        }

        DijkstraHelper<debruijn_graph::Graph>::BoundedDijkstra dijkstra(DijkstraHelper<debruijn_graph::Graph>::CreateBoundedDijkstra(g_, max_connection_length_));
        dijkstra.run(g_.EdgeEnd(e));
        for (auto v: dijkstra.ReachedVertices()) {
            for (auto connected: g_.OutgoingEdges(v)) {
                result.insert(connected);
            }
        }

        return result;
    }

    double GetWeight(debruijn_graph::EdgeId, debruijn_graph::EdgeId) const {
        return 1.0;
    }

    size_t GetLibIndex() const {
        return (size_t) -1;
    }
};


class ScaffoldGraph {

public:

    typedef debruijn_graph::EdgeId ScaffoldVertex;

    typedef size_t ScaffoldEdgeIdT;

    struct ScaffoldEdge {
    private:
        ScaffoldEdgeIdT id_;
        static std::atomic<ScaffoldEdgeIdT> scaffold_edge_id_;

        ScaffoldVertex start_;
        ScaffoldVertex end_;
        size_t color_;
        double weight_;

    public:

        ScaffoldEdge(ScaffoldVertex start, ScaffoldVertex end, size_t lib_id = (size_t) -1, double weight = 0):
                id_(scaffold_edge_id_++),
                start_(start), end_(end),
                color_(lib_id),
                weight_(weight) {
        }

        ScaffoldEdgeIdT getId() const {
            return id_;
        }


        size_t getColor() const {
            return color_;
        }

        double getWeight() const {
            return weight_;
        }

        const ScaffoldVertex getStart() const {
            return start_;
        }

        const ScaffoldVertex getEnd() const {
            return end_;
        }

        bool operator==(const ScaffoldEdge& e) const {
            return color_ == e.color_ && weight_ == e.weight_ && start_ == e.start_ && end_ == e.end_;
        }

        bool operator==(const ScaffoldEdge& e)  {
            return color_ == e.color_ && weight_ == e.weight_ && start_ == e.start_ && end_ == e.end_;
        }
    };

    typedef ScaffoldVertex VertexId;
    typedef ScaffoldEdge EdgeId;
    typedef std::set<ScaffoldVertex> VertexStorage;
    typedef std::unordered_map<ScaffoldEdgeIdT, ScaffoldEdge> EdgeStotage;
    typedef std::unordered_multimap<ScaffoldVertex, ScaffoldEdgeIdT> AdjacencyStorage;

    struct ConstScaffoldEdgeIterator {
    private:
        EdgeStotage::const_iterator iter_;

    public:
        ConstScaffoldEdgeIterator(EdgeStotage::const_iterator iter): iter_(iter) {

        }

        ScaffoldEdge operator*() const {
            return iter_->second;
        }

        const ScaffoldEdge* operator->() const {
            return &iter_->second;
        }

        ConstScaffoldEdgeIterator operator++() {
            return ConstScaffoldEdgeIterator(++iter_);
        }

        ConstScaffoldEdgeIterator operator++(int) {
            return ConstScaffoldEdgeIterator(iter_++);
        }

        bool operator==(const ConstScaffoldEdgeIterator& that) const {
            return this->iter_ == that.iter_;
        }

        bool operator!=(const ConstScaffoldEdgeIterator& that) const {
            return !this->operator==(that);
        }
    };

private:
    const debruijn_graph::Graph& assembly_graph_;

    VertexStorage vertices_;

    EdgeStotage edges_;

    std::unordered_map<ScaffoldEdgeIdT, ScaffoldEdgeIdT> conjugate_;

    AdjacencyStorage outgoing_edges_;

    AdjacencyStorage incoming_edges_;

    void AddEdgeSimple(const ScaffoldEdge& e, size_t conjugate_id) {
        edges_.emplace(e.getId(), e);
        outgoing_edges_.emplace(e.getStart(), e.getId());
        incoming_edges_.emplace(e.getEnd(), e.getId());
        conjugate_[e.getId()] = conjugate_id;
    }

    void DeleteOutgoing(const ScaffoldEdge& e) {
        auto e_range = outgoing_edges_.equal_range(e.getStart());
        for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
            if (edges_.at(edge_id->second) == e) {
                outgoing_edges_.erase(edge_id);
            }
        }
    }

    void DeleteIncoming(const ScaffoldEdge& e) {
        auto e_range = incoming_edges_.equal_range(e.getEnd());
        for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
            if (edges_.at(edge_id->second) == e) {
                incoming_edges_.erase(edge_id);
            }
        }
    }

    void DeleteEdgeFromStorage(const ScaffoldEdge& e) {
        VERIFY(!Exists(e));

        size_t conjugate_id = conjugate_[e.getId()];
        edges_.erase(e.getId());
        edges_.erase(conjugate_id);
        conjugate_.erase(e.getId());
        conjugate_.erase(conjugate_id);
    }

    void DeleteAllOutgoingEdgesSimple(ScaffoldVertex v) {
        auto e_range = outgoing_edges_.equal_range(v);
        for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
            DeleteIncoming(edges_.at(edge_id->second));
        }
        outgoing_edges_.erase(v);
    }

    void DeleteAllIncomingEdgesSimple(ScaffoldVertex v) {
        auto e_range = incoming_edges_.equal_range(v);
        for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
            DeleteOutgoing(edges_.at(edge_id->second));
        }
        incoming_edges_.erase(v);
    }


public:

    ScaffoldGraph(const debruijn_graph::Graph& g): assembly_graph_(g) {
    }

    ScaffoldGraph(const debruijn_graph::Graph& g,
                  const set<debruijn_graph::EdgeId> edge_set,
                  const vector<shared_ptr<ConnectionCondition>>& conditions):
            assembly_graph_(g) {
        vertices_.insert(edge_set.begin(), edge_set.end());
        ConstructFromConditions(conditions);
    }

    bool Exists(ScaffoldVertex assembly_graph_edge) const {
        return vertices_.count(assembly_graph_edge) != 0;
    }

    bool Exists(const ScaffoldEdge &e) const {
        auto e_range = outgoing_edges_.equal_range(e.getStart());
        for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
            if (edges_.at(edge_id->second) == e) {
                return true;
            }
        }
        return false;
    }

    ScaffoldVertex conjugate(ScaffoldVertex assembly_graph_edge) const {
        return assembly_graph_.conjugate(assembly_graph_edge);
    }

    ScaffoldEdge conjugate(const ScaffoldEdge& e) const {
        auto iter = conjugate_.find(e.getId());
        if (iter != conjugate_.end()) {
            return edges_.at(iter->second);
        }
        return ScaffoldEdge(conjugate(e.getEnd()), conjugate(e.getStart()), e.getColor(), e.getWeight());
    }

    bool AddVertex(ScaffoldVertex assembly_graph_edge) {
        if (!Exists(assembly_graph_edge)) {
            VERIFY(!Exists(conjugate(assembly_graph_edge)));
            vertices_.insert(assembly_graph_edge);
            vertices_.insert(conjugate(assembly_graph_edge));
            return true;
        }
        return false;
    }

    bool AddEdge(ScaffoldVertex v1, ScaffoldVertex v2, size_t lib_id, double weight) {
        VERIFY(Exists(v1));
        VERIFY(Exists(v2));

        ScaffoldEdge e(v1, v2, lib_id, weight);
        if (Exists(e)) {
            VERIFY(Exists(conjugate(e)));
            return false;
        }

        auto conj = conjugate(e);
        AddEdgeSimple(e, conj.getId());
        AddEdgeSimple(conj, e.getId());
        return true;
    }

    bool AddEdge(const ScaffoldEdge& e) {
        return AddEdge(e.getStart(), e.getEnd(), e.getColor(), e.getWeight());
    }

    bool RemoveEdge(const ScaffoldEdge& e) {
        if (Exists(e)) {
            VERIFY(Exists(conjugate(e)));
            DeleteOutgoing(e);
            DeleteIncoming(e);
            DeleteOutgoing(conjugate(e));
            DeleteIncoming(conjugate(e));
            DeleteEdgeFromStorage(e);

            return true;
        }
        return false;
    }

    bool RemoveVertex(ScaffoldVertex assembly_graph_edge) {
        if (Exists(assembly_graph_edge)) {
            VERIFY(Exists(conjugate(assembly_graph_edge)));

            DeleteAllOutgoingEdgesSimple(assembly_graph_edge);
            DeleteAllIncomingEdgesSimple(assembly_graph_edge);
            DeleteAllOutgoingEdgesSimple(conjugate(assembly_graph_edge));
            DeleteAllIncomingEdgesSimple(conjugate(assembly_graph_edge));

            VERIFY(incoming_edges_.count(assembly_graph_edge) == 0);
            VERIFY(outgoing_edges_.count(assembly_graph_edge) == 0);
            VERIFY(incoming_edges_.count(conjugate(assembly_graph_edge)) == 0);
            VERIFY(outgoing_edges_.count(conjugate(assembly_graph_edge)) == 0);

            vertices_.erase(assembly_graph_edge);
            vertices_.erase(conjugate(assembly_graph_edge));

            return true;
        }
        return false;
    }

    bool IsVertexIsolated(ScaffoldVertex assembly_graph_edge) const {
        bool result = incoming_edges_.count(assembly_graph_edge) == 0 && outgoing_edges_.count(assembly_graph_edge) == 0;
        VERIFY((incoming_edges_.count(conjugate(assembly_graph_edge)) == 0 && incoming_edges_.count(assembly_graph_edge) == 0) == result);
        return result;
    }

    VertexStorage::const_iterator vbegin() const {
        return vertices_.cbegin();
    }

    VertexStorage::const_iterator vend() const {
        return vertices_.cend();
    }

    ConstScaffoldEdgeIterator ebegin() const {
        return ConstScaffoldEdgeIterator(edges_.cbegin());
    }

    ConstScaffoldEdgeIterator eend() const {
        return ConstScaffoldEdgeIterator(edges_.cend());
    }

    size_t int_id(ScaffoldVertex v) const {
        return assembly_graph_.int_id(v);
    }

    size_t int_id(ScaffoldEdge e) const {
        return e.getId();
    }

    ScaffoldVertex EdgeStart(ScaffoldEdge e) const {
        return e.getStart();
    }

    ScaffoldVertex EdgeEnd(ScaffoldEdge e) const {
        return e.getEnd();
    }

    size_t VertexCount() const {
        return vertices_.size();
    }

    size_t EdgeCount() const {
        return edges_.size();
    }

    const debruijn_graph::Graph& AssemblyGraph() const {
        return assembly_graph_;
    }

    vector<ScaffoldEdge> OutgoingEdges(ScaffoldVertex assembly_graph_edge) const {
        vector<ScaffoldEdge> result;
        auto e_range = outgoing_edges_.equal_range(assembly_graph_edge);
        for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
            result.push_back(edges_.at(edge_id->second));
        }
        return result;
    }

    vector<ScaffoldEdge> IncomingEdges(ScaffoldVertex assembly_graph_edge) const {
        vector<ScaffoldEdge> result;
        auto e_range = incoming_edges_.equal_range(assembly_graph_edge);
        for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
            result.push_back(edges_.at(edge_id->second));
        }
        return result;
    }

    size_t OutgoingEdgeCount(ScaffoldVertex assembly_graph_edge) const {
        return outgoing_edges_.count(assembly_graph_edge);
    }

    size_t IncomingEdgeCount(ScaffoldVertex assembly_graph_edge) const {
        return incoming_edges_.count(assembly_graph_edge);
    }

    bool HasUniqueOutgoing(ScaffoldVertex assembly_graph_edge) const {
        return OutgoingEdgeCount(assembly_graph_edge) == 1;
    }

    bool HasUniqueIncoming(ScaffoldVertex assembly_graph_edge) const {
        return IncomingEdgeCount(assembly_graph_edge) == 1;
    }

    ScaffoldEdge UniqueOutgoing(ScaffoldVertex assembly_graph_edge) const {
        VERIFY(HasUniqueOutgoing(assembly_graph_edge));
        return edges_.at(outgoing_edges_.find(assembly_graph_edge)->second);
    }

    ScaffoldEdge UniqueIncoming(ScaffoldVertex assembly_graph_edge) const {
        VERIFY(HasUniqueIncoming(assembly_graph_edge));
        return edges_.at(incoming_edges_.find(assembly_graph_edge)->second);
    }

    void ConstructFromSingleCondition(const shared_ptr<ConnectionCondition> condition) {
        for (auto v : vertices_) {
            TRACE("Vertex " << assembly_graph_.int_id(v));
            auto connected_with = condition->ConnectedWith(v);
            for (auto connected : connected_with) {
                TRACE("Connected with " << assembly_graph_.int_id(connected));
                if (vertices_.count(connected) != 0) {
                    AddEdge(v, connected, condition->GetLibIndex(), condition->GetWeight(v, connected));
                }
            }
        }
    }

    void ConstructFromConditions(const vector<shared_ptr<ConnectionCondition>>& conditions) {
        for (auto condition : conditions) {
            ConstructFromSingleCondition(condition);
        }
    }

    void Print(ostream& os) const {
        for (auto v: vertices_) {
            os << "Vertex " << int_id(v) << " ~ " << int_id(conjugate(v))
                    << ": len = " << assembly_graph_.length(v) << ", cov = " << assembly_graph_.coverage(v) << endl;
        }
        for (auto e_iter = ebegin(); e_iter != eend(); ++e_iter) {
            os << "Edge " << e_iter->getId() << " ~ " << conjugate(*e_iter).getId() <<
                    ": " << int_id(e_iter->getStart()) << " -> " << int_id(e_iter->getEnd()) <<
                    ", lib index = " << e_iter->getColor() << ", weight " << e_iter->getWeight() << endl;
        }
    }

};


};

