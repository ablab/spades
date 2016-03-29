//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <map>
#include <queue>

#include "assembly_graph/components/splitters.hpp"
#include "cleaner.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"

namespace omnigraph {

using std::set;
using std::map;
using std::vector;
using std::pair;
using std::queue;
using std::make_pair;

template<class Graph>
class FlowGraph {
public:
    typedef size_t FlowVertexId;
    typedef pair<FlowVertexId, FlowVertexId> FlowEdgeId;

private:
    typedef typename Graph::VertexId OuterVertexId;
    map<OuterVertexId, FlowVertexId> vertex_mapping_;
    map<FlowVertexId, map<FlowVertexId, int>> capacities_;
    set<FlowVertexId> vertices_;
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

    vector<FlowEdgeId> OutgoingEdges(FlowVertexId v) const {
        vector<FlowEdgeId> result;
        const map<FlowVertexId, int> &outgoing = capacities_.find(v)->second;
        for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
            if (it->second > 0) {
                result.push_back(make_pair(v, it->first));
            }
        }
        return result;
    }

    vector<FlowEdgeId> IncomingEdges(FlowVertexId v) const {
        vector<FlowEdgeId> result;
        const map<FlowVertexId, int> &outgoing = capacities_.find(v)->second;
        for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
            if (Connected(it->first, v)) {
                result.push_back(make_pair(it->first, v));
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

    set<FlowVertexId>::iterator begin() const {
        return vertices_.begin();
    }

    set<FlowVertexId>::iterator end() const {
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

    void PushFlow(vector<FlowVertexId> path, int capacity) {
        size_t n = path.size();
        VERIFY(path[0] == source_ && path[n - 1] == sink_);
        for (size_t i = 0; i + 1 < n; i++) {
            PushFlow(make_pair(path[i], path[i + 1]), capacity);
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

    vector<FlowVertexId> RestoreAnswer(FlowVertexId start, FlowVertexId end,
            const map<FlowVertexId, FlowVertexId> &prev) {
        vector<FlowVertexId> result;
        result.push_back(end);
        FlowVertexId current = end;
        while (current != start) {
            current = prev.find(current)->second;
            result.push_back(current);
        }
        return vector<FlowVertexId>(result.rbegin(), result.rend());
    }

public:
    BFS(const Graph &graph) :
            graph_(graph) {
    }

    vector<FlowVertexId> Go(FlowVertexId start, FlowVertexId finish) {
        queue<FlowVertexId> q;
        q.push(start);
        map<FlowVertexId, FlowVertexId> prev;
        prev[start] = start;
        while (!q.empty()) {
            FlowVertexId current = q.front();
            q.pop();
            vector<FlowEdgeId> outgoing = graph_.OutgoingEdges(current);
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
        return vector<FlowVertexId>();
    }
};

template<class Graph>
class MaxFlowFinder {
private:
    FlowGraph<Graph> &graph_;
    typedef typename FlowGraph<Graph>::FlowVertexId FlowVertexId;
    typedef typename FlowGraph<Graph>::FlowEdgeId FlowEdgeId;

    int MinCapacity(vector<FlowVertexId> path) {
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
            vector<FlowVertexId> path = bfs.Go(graph_.Source(), graph_.Sink());
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

    void Find(FlowVertexId v, vector<FlowVertexId> &result, set<FlowVertexId> &visited) {
        visited.insert(v);
        vector<FlowEdgeId> outgoing = graph_.OutgoingEdges(v);
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

    vector<FlowVertexId> Sort() {
        vector<FlowVertexId> result;
        set<FlowVertexId> visited;
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

    void Find(FlowVertexId v, map<FlowVertexId, size_t> &result, size_t cc) {
        result[v] = cc;
        vector<FlowEdgeId> incoming = graph_.IncomingEdges(v);
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

    map<FlowVertexId, size_t> Find(const vector<FlowVertexId> &order) {
        size_t cc = 0;
        map<FlowVertexId, size_t> result;
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

    map<FlowVertexId, size_t> ColourComponents() {
        map<FlowVertexId, size_t> result;
        vector<FlowVertexId> order = TopSorter<Graph>(graph_).Sort();
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

    set<EdgeId> CollectUnusedEdges(set<VertexId> component, FlowGraph<Graph> fg,
            const map<typename FlowGraph<Graph>::FlowVertexId, size_t> &colouring) {
        set<EdgeId> result;
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

    bool IsInnerShortEdge(set<VertexId> component, EdgeId edge) {
        return !IsUnique(edge) && component.count(g_.EdgeStart(edge)) == 1
                && component.count(g_.EdgeEnd(edge)) == 1;
    }

    void ProcessShortEdge(FlowGraph<Graph> &fg, set<VertexId> component,
            EdgeId edge) {
        if (IsInnerShortEdge(component, edge)) {
            fg.AddEdge(g_.EdgeStart(edge), g_.EdgeEnd(edge));
        }
    }

    void ProcessSource(FlowGraph<Graph> &fg, set<VertexId> /*component*/,
            EdgeId edge) {
        if (IsPlausible(edge) || IsUnique(edge)) {
            fg.AddSource(g_.EdgeEnd(edge), 1);
        }
    }

    void ProcessSink(FlowGraph<Graph> &fg, set<VertexId> /*component*/,
            EdgeId edge) {
        if (IsPlausible(edge) || IsUnique(edge)) {
            fg.AddSink(g_.EdgeStart(edge), 1);
        }
    }

    void ConstructFlowGraph(FlowGraph<Graph> &fg, set<VertexId> component) {
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
                    plausibility_length), component_remover_(g, (std::function<void (set<EdgeId>)>) 0) {
        VERIFY(uniqueness_length >= plausibility_length);
        VERIFY(plausibility_length > max_length);
    }

    bool Process() {
        for (shared_ptr<GraphSplitter<Graph>> splitter_ptr = LongEdgesExclusiveSplitter<Graph>(g_,
                uniqueness_length_); splitter_ptr->HasNext();) {
            set<VertexId> component = splitter_ptr->Next().vertices();
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
            map<typename FlowGraph<Graph>::FlowVertexId, size_t> colouring =
                    component_finder.ColourComponents();
            set<EdgeId> to_remove = CollectUnusedEdges(component, fg,
                    colouring);
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
