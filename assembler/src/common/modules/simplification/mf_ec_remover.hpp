//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/components/splitters.hpp"
#include "cleaner.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"

#include <map>
#include <queue>

namespace omnigraph {

template<class Graph>
class FlowGraph {
public:
    typedef size_t FlowVertexId;
    typedef std::pair<FlowVertexId, FlowVertexId> FlowEdgeId;

private:
    typedef typename Graph::VertexId OuterVertexId;
    std::map<OuterVertexId, FlowVertexId> vertex_mapping_;
    std::map<FlowVertexId, std::map<FlowVertexId, int>> capacities_;
    std::set<FlowVertexId> vertices_;
    size_t vertex_number_;
    FlowVertexId source_;
    FlowVertexId sink_;

    FlowVertexId AddVertex() {
        vertices_.insert(vertex_number_);
        capacities_[vertex_number_];
        vertex_number_++;
        return vertex_number_ - 1;
    }

    void PushFlow(FlowEdgeId edge, int capacity) {
        VERIFY(capacities_[EdgeStart(edge)][EdgeEnd(edge)] >= capacity);
        capacities_[EdgeStart(edge)][EdgeEnd(edge)] -= capacity;
        capacities_[EdgeEnd(edge)][EdgeStart(edge)] += capacity;
    }

    void AddEdge(FlowVertexId first, FlowVertexId second, int capacity = 10000) {
        capacities_[first][second] += capacity; // operator [] creates entry with default values in case argument is not in keyset
        capacities_[second][first] += 0;
    }

public:
    FlowGraph() :
            vertex_number_(0), source_(AddVertex()), sink_(AddVertex()) {
    }

    FlowVertexId GetCorrespondingVertex(OuterVertexId v) const {
        return vertex_mapping_.find(v)->second;
    }

    bool HasCorrespondingVertex(OuterVertexId v) const {
        return vertex_mapping_.find(v) == vertex_mapping_.end();
    }

    FlowVertexId AddVertex(OuterVertexId vertex) {
        FlowVertexId new_vertex = AddVertex();
        vertex_mapping_[vertex] = new_vertex;
        return new_vertex;
    }

    void AddEdge(OuterVertexId outer_first, OuterVertexId outer_second,
            int capacity = 10000) {
        VERIFY(
                vertex_mapping_.find(outer_first) != vertex_mapping_.end()
                        && vertex_mapping_.find(outer_second)
                                != vertex_mapping_.end());
        FlowVertexId first = vertex_mapping_[outer_first];
        FlowVertexId second = vertex_mapping_[outer_second];
        AddEdge(first, second, capacity);
    }

    void AddSource(OuterVertexId vertex, int capacity) {
        AddEdge(source_, GetCorrespondingVertex(vertex), capacity);
    }

    void AddSink(OuterVertexId vertex, int capacity) {
        AddEdge(GetCorrespondingVertex(vertex), sink_, capacity);
    }

    FlowVertexId Source() const {
        return source_;
    }

    FlowVertexId Sink() const {
        return sink_;
    }

    bool Connected(FlowVertexId start, FlowVertexId end) const {
        return capacities_.find(start) != capacities_.end()
                && capacities_.find(start)->second.find(end)
                        != capacities_.find(start)->second.end()
                && capacities_.find(start)->second.find(end)->second > 0;
    }

    std::vector<FlowEdgeId> OutgoingEdges(FlowVertexId v) const {
        std::vector<FlowEdgeId> result;
        const auto &outgoing = capacities_.find(v)->second;
        for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
            if (it->second > 0) {
                result.emplace_back(v, it->first);
            }
        }
        return result;
    }

    std::vector<FlowEdgeId> IncomingEdges(FlowVertexId v) const {
        std::vector<FlowEdgeId> result;
        const auto &outgoing = capacities_.find(v)->second;
        for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
            if (Connected(it->first, v)) {
                result.emplace_back(it->first, v);
            }
        }
        return result;
    }

    size_t OutgoingEdgesCount(FlowVertexId v) const {
        return OutgoingEdges(v).size();
    }

    size_t IncomingEdgesCount(FlowVertexId v) const {
        return IncomingEdges(v).size();
    }

    FlowVertexId EdgeStart(FlowEdgeId edge) const {
        return edge.first;
    }

    FlowVertexId EdgeEnd(FlowEdgeId edge) const {
        return edge.second;
    }

    std::set<FlowVertexId>::iterator begin() const {
        return vertices_.begin();
    }

    std::set<FlowVertexId>::iterator end() const {
        return vertices_.end();
    }

    int GetCapacity(FlowVertexId first, FlowVertexId second) const {
        auto it1 = capacities_.find(first);
        if (it1 == capacities_.end())
            return 0;
        auto it2 = it1->second.find(second);
        if (it2 == it1->second.end())
            return 0;
        return it2->second;
    }

    void PushFlow(std::vector<FlowVertexId> path, int capacity) {
        size_t n = path.size();
        VERIFY(path[0] == source_ && path[n - 1] == sink_);
        for (size_t i = 0; i + 1 < n; i++) {
            PushFlow(std::make_pair(path[i], path[i + 1]), capacity);
        }
    }

//    void Print() const {
//        for(auto it = vertex_mapping_.begin(); it != vertex_mapping_.end(); ++it) {
//            TRACE(it->first << " " << it->second);
//        }
//        for(auto it = vertices_.begin(); it != vertices_.end();) {
//            auto out = OutgoingEdges(*it);
//            for(auto it1 = out.begin(); it1 != out.end(); ++it1) {
//                TRACE("edge " << (*it1) << " " << GetCapacity(*it, it1->second));
//            }
//            ++it;
//            if(it == vertices_.end())
//                break;
//        }
//    }
};

template<class Graph>
class BFS {
private:
    const Graph &graph_;
    typedef typename Graph::FlowVertexId FlowVertexId;
    typedef typename Graph::FlowEdgeId FlowEdgeId;

    std::vector<FlowVertexId> RestoreAnswer(FlowVertexId start, FlowVertexId end,
                                            const std::map<FlowVertexId, FlowVertexId> &prev) {
        std::vector<FlowVertexId> result = {end};
        FlowVertexId current = end;
        while (current != start) {
            current = prev.find(current)->second;
            result.push_back(current);
        }
        return std::vector<FlowVertexId>(result.rbegin(), result.rend());
    }

public:
    BFS(const Graph &graph) :
            graph_(graph) {
    }

    std::vector<FlowVertexId> Go(FlowVertexId start, FlowVertexId finish) {
        std::queue<FlowVertexId> q;
        q.push(start);
        std::map<FlowVertexId, FlowVertexId> prev;
        prev[start] = start;
        while (!q.empty()) {
            FlowVertexId current = q.front();
            q.pop();
            auto outgoing = graph_.OutgoingEdges(current);
            for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
                if (prev.find(it->second) == prev.end()) {
                    q.push(it->second);
                    prev[it->second] = current;
                }
                if (it->second == finish) {
                    return RestoreAnswer(start, finish, prev);
                }
            }
        }
        return {};
    }
};

template<class Graph>
class MaxFlowFinder {
private:
    FlowGraph<Graph> &graph_;
    typedef typename FlowGraph<Graph>::FlowVertexId FlowVertexId;
    typedef typename FlowGraph<Graph>::FlowEdgeId FlowEdgeId;

    int MinCapacity(const std::vector<FlowVertexId> &path) {
        VERIFY(path.size() >= 2);
        int result = graph_.GetCapacity(path[0], path[1]);
        for (size_t i = 1; i + 1 < path.size(); i++) {
            result = std::min(result, graph_.GetCapacity(path[i], path[i + 1]));
        }
        return result;
    }

public:
    MaxFlowFinder(FlowGraph<Graph> &graph) :
            graph_(graph) {
    }

    void Find() {
        BFS<FlowGraph<Graph> > bfs(graph_);
        while (true) {
            auto path = bfs.Go(graph_.Source(), graph_.Sink());
            if (path.size() == 0)
                break;
            int capacity = MinCapacity(path);
            VERIFY(capacity > 0);
            graph_.PushFlow(path, capacity);
//            graph_.Print();
        }
    }
};

template<class Graph>
class TopSorter {
private:
    typedef typename Graph::FlowVertexId FlowVertexId;
    typedef typename Graph::FlowEdgeId FlowEdgeId;
    const Graph &graph_;

    void Find(FlowVertexId v, std::vector<FlowVertexId> &result, std::set<FlowVertexId> &visited) {
        visited.insert(v);
        auto outgoing = graph_.OutgoingEdges(v);
        for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
            FlowVertexId next = graph_.EdgeEnd(*it);
            if (visited.count(next) == 0) {
                Find(next, result, visited);
            }
        }
        result.push_back(v);
    }

public:
    TopSorter(const Graph &graph) :
            graph_(graph) {
    }

    std::vector<FlowVertexId> Sort() {
        std::vector<FlowVertexId> result;
        std::set<FlowVertexId> visited;
        for (auto it = graph_.begin(); it != graph_.end(); ++it) {
            if (visited.count(*it) == 0) {
                Find(*it, result, visited);
            }
        }
        return result;
    }
};

template<class Graph>
class ReverseDFSComponentFinder {
private:
    typedef typename Graph::FlowVertexId FlowVertexId;
    typedef typename Graph::FlowEdgeId FlowEdgeId;

    const Graph &graph_;

    void Find(FlowVertexId v, std::map<FlowVertexId, size_t> &result, size_t cc) {
        result[v] = cc;
        auto incoming = graph_.IncomingEdges(v);
        for (auto it = incoming.begin(); it != incoming.end(); ++it) {
            FlowVertexId next = graph_.EdgeStart(*it);
            if (result.count(next) == 0) {
                Find(next, result, cc);
            }
        }
    }
public:
    ReverseDFSComponentFinder(const Graph &graph) :
            graph_(graph) {
    }

    std::map<FlowVertexId, size_t> Find(const std::vector<FlowVertexId> &order) {
        size_t cc = 0;
        std::map<FlowVertexId, size_t> result;
        for (auto it = order.rbegin(); it != order.rend(); ++it) {
            if (result.count(*it) == 0) {
                Find(*it, result, cc);
                cc++;
            }
        }
        return result;
    }
};

template<class Graph>
class StroglyConnectedComponentFinder {
private:
    typedef typename Graph::FlowVertexId FlowVertexId;
    const Graph &graph_;
    bool ready_;
public:
    StroglyConnectedComponentFinder(const Graph &graph) :
            graph_(graph), ready_(false) {
    }

    std::map<FlowVertexId, size_t> ColourComponents() {
        std::map<FlowVertexId, size_t> result;
        auto order = TopSorter<Graph>(graph_).Sort();
        return ReverseDFSComponentFinder<Graph>(graph_).Find(order);
    }
};

template<class Graph>
class MaxFlowECRemover {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    Graph& g_;
    size_t max_length_;
    size_t uniqueness_length_;
    size_t plausibility_length_;
    ComponentRemover<Graph> component_remover_;

    bool IsTerminal(VertexId vertex) {
        return g_.OutgoingEdgeCount(vertex)
                + g_.IncomingEdgeCount(vertex) == 1;
    }

    bool IsTip(EdgeId edge) {
        VertexId start = g_.EdgeStart(edge);
        VertexId end = g_.EdgeEnd(edge);
        return IsTerminal(start) || IsTerminal(end);
    }


    bool IsSuspicious(EdgeId edge) {
        return g_.length(edge) <= max_length_ && !IsTip(edge);
    }

    std::set<EdgeId> CollectUnusedEdges(const std::set<VertexId> &component, const FlowGraph<Graph> &fg,
                                        const std::map<typename FlowGraph<Graph>::FlowVertexId, size_t> &colouring) {
        std::set<EdgeId> result;
        for (auto it_start = component.begin(); it_start != component.end();
                ++it_start) {
            VertexId start = *it_start;
            auto outgoing = g_.OutgoingEdges(start);
            for (auto it_edge = outgoing.begin(); it_edge != outgoing.end();
                    ++it_edge) {
                EdgeId edge = *it_edge;
                VertexId end = g_.EdgeEnd(edge);
                if (component.count(end) == 1 && IsSuspicious(edge)
                        && colouring.find(fg.GetCorrespondingVertex(start))->second
                                != colouring.find(
                                        fg.GetCorrespondingVertex(end))->second) {
                    result.insert(edge);
                }
            }
        }
        return result;
    }

    bool CheckCompleteFlow(FlowGraph<Graph> &fg) {
        return fg.OutgoingEdges(fg.Source()).size() == 0
                && fg.IncomingEdges(fg.Sink()).size() == 0;
    }

    bool IsPlausible(EdgeId edge) {
        return g_.length(edge) >= plausibility_length_ && !IsTip(edge);
    }

    bool IsUnique(EdgeId edge) {
        return g_.length(edge) >= uniqueness_length_;
    }

    bool IsInnerShortEdge(const std::set<VertexId> &component, EdgeId edge) {
        return !IsUnique(edge) && component.count(g_.EdgeStart(edge)) == 1
                && component.count(g_.EdgeEnd(edge)) == 1;
    }

    void ProcessShortEdge(FlowGraph<Graph> &fg, const std::set<VertexId> &component, EdgeId edge) {
        if (IsInnerShortEdge(component, edge)) {
            fg.AddEdge(g_.EdgeStart(edge), g_.EdgeEnd(edge));
        }
    }

    void ProcessSource(FlowGraph<Graph> &fg, const std::set<VertexId> &/*component*/, EdgeId edge) {
        if (IsPlausible(edge) || IsUnique(edge)) {
            fg.AddSource(g_.EdgeEnd(edge), 1);
        }
    }

    void ProcessSink(FlowGraph<Graph> &fg, const std::set<VertexId> &/*component*/, EdgeId edge) {
        if (IsPlausible(edge) || IsUnique(edge)) {
            fg.AddSink(g_.EdgeStart(edge), 1);
        }
    }

    void ConstructFlowGraph(FlowGraph<Graph> &fg, const std::set<VertexId> &component) {
        for (auto it = component.begin(); it != component.end(); ++it) {
            fg.AddVertex(*it);
        }
        for (auto it = component.begin(); it != component.end(); ++it) {
            VertexId vertex = *it;
            auto outgoing = g_.OutgoingEdges(vertex);
            for (auto it_edge = outgoing.begin(); it_edge != outgoing.end();
                    ++it_edge) {
                EdgeId edge = *it_edge;
                ProcessShortEdge(fg, component, edge);
                ProcessSink(fg, component, edge);
            }
            auto incoming = g_.IncomingEdges(vertex);
            for (auto it_edge = incoming.begin(); it_edge != incoming.end();
                    ++it_edge) {
                EdgeId edge = *it_edge;
                ProcessSource(fg, component, edge);
            }
        }
    }

public:
    MaxFlowECRemover(Graph& g, size_t max_length, size_t uniqueness_length,
            size_t plausibility_length, std::function<void (EdgeId)>
    /*fixme ignored, fix after merge with relative coverage branch!!! removal_handler*/) :
            g_(g), max_length_(max_length), uniqueness_length_(
                    uniqueness_length), plausibility_length_(
                    plausibility_length), component_remover_(g, typename ComponentRemover<Graph>::HandlerF(nullptr)) {
        VERIFY(uniqueness_length >= plausibility_length);
        VERIFY(plausibility_length > max_length);
    }

    bool Process() {
        for (std::shared_ptr<GraphSplitter<Graph>> splitter_ptr =
                     LongEdgesExclusiveSplitter<Graph>(g_, uniqueness_length_); splitter_ptr->HasNext();) {
            auto component = splitter_ptr->Next().vertices();
            FlowGraph<Graph> fg;
            ConstructFlowGraph(fg, component);
//            fg.Print();
            MaxFlowFinder<Graph> mf_finder(fg);
            mf_finder.Find();
            if (!CheckCompleteFlow(fg)) {
                TRACE("Suspicious component! No edge delition!");
                continue;
            }
            StroglyConnectedComponentFinder<FlowGraph<Graph>> component_finder(
                    fg);
            auto colouring = component_finder.ColourComponents();
            auto to_remove = CollectUnusedEdges(component, fg, colouring);
            component_remover_.DeleteComponent(to_remove.begin(), to_remove.end(), false);
        }
        CompressAllVertices(g_);
        Cleaner<Graph>(g_).Run();

        return false;
    }
private:
    DECL_LOGGER("MaxFlowECRemover");
};
}
