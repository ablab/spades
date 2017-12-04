#include "scaffold_graph_storage.hpp"

namespace path_extend {
ScaffoldGraphStorage::ScaffoldGraphStorage(ScaffoldGraphStorage::ScaffoldGraph&& large_scaffold_graph,
                                           ScaffoldGraphStorage::ScaffoldGraph&& small_scaffold_graph)
    : large_scaffold_graph_(large_scaffold_graph), small_scaffold_graph_(small_scaffold_graph) {}

const ScaffoldGraphStorage::ScaffoldGraph& ScaffoldGraphStorage::GetLargeScaffoldGraph() const {
    return large_scaffold_graph_;
}

const ScaffoldGraphStorage::ScaffoldGraph& ScaffoldGraphStorage::GetSmallScaffoldGraph() const {
    return small_scaffold_graph_;
}

void ScaffoldGraphStorage::SetSmallScaffoldGraph(const ScaffoldGraph& small_scaffold_graph) {
    ReplaceScaffoldGraph(small_scaffold_graph, small_scaffold_graph_);
}

void ScaffoldGraphStorage::SetLargeScaffoldGraph(const ScaffoldGraph& large_scaffold_graph) {
    ReplaceScaffoldGraph(large_scaffold_graph, large_scaffold_graph_);
}

ScaffoldGraphStorage::ScaffoldGraphStorage(const debruijn_graph::Graph& g) : large_scaffold_graph_(g), small_scaffold_graph_(g) {}

void ScaffoldGraphStorage::ReplaceScaffoldGraph(const ScaffoldGraphStorage::ScaffoldGraph &from, ScaffoldGraph &to) {
    to = from;
    VERIFY(to.EdgeCount() == from.EdgeCount());
    VERIFY(to.VertexCount() == from.VertexCount());
    INFO("Finished replacing");
}

void ScaffoldGraphStorage::Load(const string& path, const std::map<size_t, debruijn_graph::EdgeId>& edge_map) {
    ifstream fin(path);
    ScaffoldGraphSerializer loader;
    loader.LoadScaffoldGraph(fin, large_scaffold_graph_, edge_map);
    loader.LoadScaffoldGraph(fin, small_scaffold_graph_, edge_map);
}
void ScaffoldGraphStorage::Save(const string& path) const {
    ofstream fout(path);
    ScaffoldGraphSerializer saver;
    saver.SaveScaffoldGraph(fout, large_scaffold_graph_);
    saver.SaveScaffoldGraph(fout, small_scaffold_graph_);
}
void ScaffoldGraphSerializer::SaveScaffoldGraph(ofstream &fout, const ScaffoldGraphSerializer::ScaffoldGraph &graph) const {
    fout << graph.VertexCount() << std::endl;
    for (const ScaffoldGraph::ScaffoldGraphVertex& vertex: graph.vertices()) {
        fout << vertex.int_id() << std::endl;
    }
    fout << graph.EdgeCount() << std::endl;
    for (const ScaffoldGraph::ScaffoldEdge& edge: graph.edges()) {
        fout << edge.getStart().int_id() << " " << edge.getEnd().int_id() << " " << edge.getColor() << " "
             << edge.getWeight() << " " << edge.getLength() << std::endl;
    }
}
void ScaffoldGraphSerializer::LoadScaffoldGraph(ifstream &fin,
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
    for (size_t i = 0; i < number_of_edges; ++i) {
        size_t start_id = 0;
        size_t end_id = 0;
        size_t lib = -1;
        double weight = 0;
        size_t length = 0;
        fin >> start_id >> end_id >> lib >> weight >> length;
        ScaffoldGraph::ScaffoldEdge new_edge(edge_map.at(start_id), edge_map.at(end_id), lib, weight, length);
        graph.AddEdge(new_edge);
    }
}
}