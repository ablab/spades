//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_graph_storage.hpp"

namespace path_extend {
namespace read_cloud {

ScaffoldGraphStorage::ScaffoldGraphStorage(ScaffoldGraphStorage::ScaffoldGraph &&large_scaffold_graph,
                                           ScaffoldGraphStorage::ScaffoldGraph &&small_scaffold_graph,
                                           size_t large_length_threshold, size_t small_length_threshold)
    : large_scaffold_graph_(large_scaffold_graph), small_scaffold_graph_(small_scaffold_graph),
      large_length_threshold_(large_length_threshold), small_length_threshold_(small_length_threshold) {}

const ScaffoldGraphStorage::ScaffoldGraph &ScaffoldGraphStorage::GetLargeScaffoldGraph() const {
    return large_scaffold_graph_;
}

const ScaffoldGraphStorage::ScaffoldGraph &ScaffoldGraphStorage::GetSmallScaffoldGraph() const {
    return small_scaffold_graph_;
}

void ScaffoldGraphStorage::SetSmallScaffoldGraph(const ScaffoldGraph &small_scaffold_graph) {
    ReplaceScaffoldGraph(small_scaffold_graph, small_scaffold_graph_);
}

void ScaffoldGraphStorage::SetLargeScaffoldGraph(const ScaffoldGraph &large_scaffold_graph) {
    ReplaceScaffoldGraph(large_scaffold_graph, large_scaffold_graph_);
}

ScaffoldGraphStorage::ScaffoldGraphStorage(const debruijn_graph::Graph &g) :
    large_scaffold_graph_(g), small_scaffold_graph_(g), large_length_threshold_(0), small_length_threshold_(0) {}

void ScaffoldGraphStorage::ReplaceScaffoldGraph(const ScaffoldGraphStorage::ScaffoldGraph &from, ScaffoldGraph &to) {
    to = from;
    VERIFY(to.EdgeCount() == from.EdgeCount());
    VERIFY(to.VertexCount() == from.VertexCount());
    INFO("Finished replacing");
}

void ScaffoldGraphStorage::Load(const std::string &path, const std::map<size_t, debruijn_graph::EdgeId> &edge_map) {
    std::ifstream fin(path);
    ScaffoldGraphSerializer loader;
    fin >> small_length_threshold_ >> large_length_threshold_;
    loader.LoadScaffoldGraph(fin, large_scaffold_graph_, edge_map);
    loader.LoadScaffoldGraph(fin, small_scaffold_graph_, edge_map);
}
void ScaffoldGraphStorage::Save(const std::string &path) const {
    std::ofstream fout(path);
    ScaffoldGraphSerializer saver;
    fout << small_length_threshold_ << " " << large_length_threshold_ << std::endl;
    saver.SaveScaffoldGraph(fout, large_scaffold_graph_);
    saver.SaveScaffoldGraph(fout, small_scaffold_graph_);
}
ScaffoldGraphStorage &ScaffoldGraphStorage::operator=(const ScaffoldGraphStorage &other) {
    SetSmallScaffoldGraph(other.small_scaffold_graph_);
    small_length_threshold_ = other.small_length_threshold_;
    SetLargeScaffoldGraph(other.large_scaffold_graph_);
    large_length_threshold_ = other.large_length_threshold_;
    return *this;
}
size_t ScaffoldGraphStorage::GetLargeLengthThreshold() const {
    return large_length_threshold_;
}
size_t ScaffoldGraphStorage::GetSmallLengthThreshold() const {
    return small_length_threshold_;
}
void ScaffoldGraphSerializer::SaveScaffoldGraph(std::ofstream &fout,
                                                const ScaffoldGraphSerializer::ScaffoldGraph &graph) const {
    fout << graph.VertexCount() << std::endl;
    for (const ScaffoldGraph::ScaffoldGraphVertex &vertex: graph.vertices()) {
        fout << vertex.int_id() << std::endl;
    }
    fout << graph.EdgeCount() << std::endl;
    for (const ScaffoldGraph::ScaffoldEdge &edge: graph.edges()) {
        fout << edge.getStart().int_id() << " " << edge.getEnd().int_id() << " " << edge.getColor() << " "
             << edge.getWeight() << " " << edge.getLength() << std::endl;
    }
}
void ScaffoldGraphSerializer::LoadScaffoldGraph(std::ifstream &fin,
                                                ScaffoldGraphSerializer::ScaffoldGraph &graph,
                                                const std::map<size_t, debruijn_graph::EdgeId> &edge_map) const {
    size_t number_of_vertices = 0;
    fin >> number_of_vertices;
    for (size_t i = 0; i < number_of_vertices; ++i) {
        size_t edge_id = 0;
        fin >> edge_id;
        EdgeId e = edge_map.at(edge_id);
        graph.AddVertex(e);
    }
    size_t number_of_edges = 0;
    fin >> number_of_edges;
    DEBUG(number_of_edges << " edges");
    size_t block_size = number_of_edges / 20;
    std::vector<ScaffoldGraph::ScaffoldEdge> edges;
    for (size_t i = 0; i < number_of_edges; ++i) {
        TRACE("Start")
        size_t start_id = 0;
        size_t end_id = 0;
        size_t lib = -1;
        double weight = 0;
        size_t length = 0;
        fin >> start_id >> end_id >> lib >> weight >> length;
        if (block_size != 0 and i % block_size == 0) {
            DEBUG("Loaded " << i << " edges out of " << number_of_edges);
        }
        TRACE("Read values")
        edges.emplace_back(edge_map.at(start_id), edge_map.at(end_id), lib, weight, length);
    }
    DEBUG(edges.size() << " edges loaded");
    DEBUG("Adding edges");
    size_t current = 0;
    block_size = edges.size() / 200;
    for (const auto &edge: edges) {
        graph.AddEdgeSimple(edge);
        ++current;
        if (block_size != 0 and current % block_size == 0) {
            DEBUG("Added " << current << " edges out of " << edges.size());
        }
    }
}
}
}