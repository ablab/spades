//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "path_scaffolder.hpp"

#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"

namespace path_extend {

void SimplePathScaffolder::CondenseSimplePaths(const std::vector<ScaffoldEdge> &scaffold_edges) const {
    std::unordered_map<ScaffoldVertex, ScaffoldVertex> merge_connections;
    for (const auto &edge: scaffold_edges) {
        ScaffoldVertex start = edge.getStart();
        ScaffoldVertex end = edge.getEnd();
        DEBUG(start.int_id() << " -> " << end.int_id());
        DEBUG("Weight: " << edge.getWeight());
        VERIFY(merge_connections.find(start) == merge_connections.end());
        merge_connections.insert({start, end});
    }

    for (const auto &connection: merge_connections) {
        auto start = connection.first;
        auto end = connection.second;
        auto start_conjugate = start.GetConjugateFromGraph(gp_.g);
        auto end_conjugate = end.GetConjugateFromGraph(gp_.g);
        if (merge_connections.find(end_conjugate) == merge_connections.end() or
            merge_connections.at(end_conjugate) != start_conjugate) {
            WARN("Conjugate connection does not correspond to direct connection")
            merge_connections.at(end_conjugate) = start_conjugate;
        } else {
            merge_connections.insert({end_conjugate, start_conjugate});
        }
    }

    StartFinder start_finder(gp_.g);
    auto starts = start_finder.GetStarts(merge_connections);
    std::unordered_map<ScaffoldVertex, size_t> start_to_distance;
    for (const auto &edge: scaffold_edges) {
        start_to_distance.insert({edge.getStart(), edge.getLength()});
        start_to_distance.insert({edge.getEnd().GetConjugateFromGraph(gp_.g), edge.getLength()});
    }
    for (const auto &connection: merge_connections) {
        DEBUG(connection.first.int_id() << " -> " << connection.second.int_id());
    }
    INFO(starts.size() << " starts.");
    for (const auto &start: starts) {
        ScaffoldVertex current = start;
        bool next_found = merge_connections.find(current) != merge_connections.end();
        DEBUG("Start: " << current.int_id());
        while (next_found and merge_connections.at(current) != start) {
            current = merge_connections.at(current);
            next_found = merge_connections.find(current) != merge_connections.end();
            DEBUG(current.int_id());
        }
    }
    for (const auto &start: starts) {
        if (not start.ToPath(gp_.g)->Empty()) {
            ExtendPathAlongConnections(start, merge_connections, start_to_distance);
        }
    }
}

void SimplePathScaffolder::MergePaths(const ScaffoldGraph &scaffold_graph) const {
    INFO(scaffold_graph.VertexCount() << " vertices and " << scaffold_graph.EdgeCount()
                                           << " edges in path scaffold graph");
    for (const ScaffoldVertex &vertex: scaffold_graph.vertices()) {
        VERIFY_DEV(vertex.GetType() == scaffold_graph::ScaffoldVertexT::Path);
    }
    read_cloud::ScaffoldGraphExtractor graph_extractor;
    auto reliable_edges = graph_extractor.ExtractReliableEdges(scaffold_graph);
    INFO("Found " << reliable_edges.size() << " reliable edges");
    CondenseSimplePaths(reliable_edges);
}

std::unordered_set<StartFinder::ScaffoldVertex> StartFinder::GetStarts(const TransitionMap &transition_map) const {
    std::unordered_set<ScaffoldVertex> starts;
    std::unordered_set<ScaffoldVertex> used;
    for (const auto &connection: transition_map) {
        auto start = connection.first;
        auto current = start;
        auto current_conjugate = current.GetConjugateFromGraph(g_);
        if (used.find(current) != used.end()) {
            continue;
        }
        bool prev_found = transition_map.find(current_conjugate) != transition_map.end();
        bool prev_used = false;
        while (prev_found) {
            used.insert(current);
            used.insert(current_conjugate);
            auto prev_conjugate = transition_map.at(current_conjugate);
            if (used.find(prev_conjugate) != used.end()) {
                prev_used = true;
                break;
            }
            current = prev_conjugate.GetConjugateFromGraph(g_);
            current_conjugate = current.GetConjugateFromGraph(g_);
            prev_found = transition_map.find(current_conjugate) != transition_map.end();
        }
        starts.insert(current);
        if (not prev_used) {
            bool next_found = transition_map.find(current) != transition_map.end();
            while (next_found) {
                current = transition_map.at(current);
                used.insert(current);
                used.insert(current.GetConjugateFromGraph(g_));
                next_found = transition_map.find(current) != transition_map.end();
            }
        } else {
            VERIFY_DEV(used.find(start) != used.end());
        }
    }
    return starts;
}

void SimplePathScaffolder::ExtendPathAlongConnections(
        const ScaffoldVertex &start,
        const std::unordered_map<ScaffoldVertex,ScaffoldVertex> &merge_connections,
        const std::unordered_map<ScaffoldVertex, size_t> &start_to_distance) const {
    auto current = start;
    bool next_found = merge_connections.find(current) != merge_connections.end();
    auto start_path = start.ToPath(gp_.g);
    while (next_found) {
        auto next = merge_connections.at(current);
        auto next_path = next.ToPath(gp_.g);
        if (start_path->GetId() == next_path->GetId()) {
            break;
        }
        DEBUG("First path: " << start_path->GetId() << ", length : " << start_path->Length());
        DEBUG("Second path: " << next_path->GetId() << ", length: " << next_path->Length());
        DEBUG("First conj: " << start_path->GetConjPath()->GetId() << ", length : "
                             << start_path->GetConjPath()->Length());
        DEBUG(
            "Second conj: " << next_path->GetConjPath()->GetId() << ", length: " << next_path->GetConjPath()->Length());
        DEBUG("Got paths")
        int gap_length = static_cast<int>(start_to_distance.at(current));
        if (gap_length == 0) {
            gap_length = default_gap_;
        }
        Gap path_distance_gap(gap_length);
        DEBUG("Push back")
        start_path->PushBack(*next_path, path_distance_gap);
        DEBUG("Clear");
        next_path->Clear();
        DEBUG("Second path: " << next_path->GetId() << ", length: " << next_path->Length());
        DEBUG(next_path->Empty());
        DEBUG("Conjugate: " << next_path->GetConjPath()->GetId() << ", length: " << next_path->GetConjPath()->Length());
        DEBUG("Conjugate empty: " << next_path->GetConjPath()->Empty());
        current = next;
        next_found = merge_connections.find(current) != merge_connections.end();
    }
}
SimplePathScaffolder::SimplePathScaffolder(const conj_graph_pack &gp, int default_gap) :
    gp_(gp), default_gap_(default_gap) {}

StartFinder::StartFinder(const Graph &g): g_(g) {}
}