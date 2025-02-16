//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2023-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/construction_helper.hpp"
#include "assembly_graph/graph_support/v_overlaps_support.hpp"
#include "io/binary/graph.hpp"
#include "io/graph/gfa_reader.hpp"
#include "io/reads/file_reader.hpp"

#include <filesystem>
#include <gtest/gtest.h>

using namespace debruijn_graph;

using IdMapper = io::IdMapper<std::string>;

void CheckSegmentLen(std::ifstream &graph_stream, const Graph &graph, const IdMapper &id_mapper) {
    int num_edges;
    graph_stream >> num_edges;
    for (int i = 0; i < num_edges; ++i) {
        std::string edge_name;
        int edge_len;
        graph_stream >> edge_name >> edge_len;
        EdgeId edge = id_mapper[edge_name];
        EXPECT_EQ(graph.length(edge) + graph.k(), edge_len);

    }
}

std::unordered_map<VertexId, std::string> CheckStructure(std::ifstream &graph_stream,
                                                         const Graph &graph,
                                                         const IdMapper &id_mapper) {
    std::unordered_map<std::string, int> vertex_to_len;
    std::unordered_map<std::string, std::string> edge_to_end;
    std::unordered_map<std::string, std::string> edge_to_start;
    std::unordered_map<VertexId, std::string> vid_to_name;

    int num_vertices;
    graph_stream >> num_vertices;
    for (int i = 0; i < num_vertices; ++i) {
        std::string vertex_name;
        int vertex_len;
        graph_stream >> vertex_name >> vertex_len;
        vertex_to_len[vertex_name] = vertex_len;
    }

    int num_incidence_links;
    graph_stream >> num_incidence_links;
    for (int i = 0; i < num_incidence_links; ++i) {
        std::string first_name, second_name;
        graph_stream >> first_name >> second_name;
        if (vertex_to_len.find(first_name) != vertex_to_len.end()) {
            edge_to_start[second_name] = first_name;
            vid_to_name[graph.EdgeStart(id_mapper[second_name])] = first_name;
        } else {
            edge_to_end[first_name] = second_name;
            vid_to_name[graph.EdgeEnd(id_mapper[first_name])] = second_name;
        }
    }
    int checked_incidence_links = 0;
    for (const auto &entry: vid_to_name) {
        VertexId vertex = entry.first;
        std::string name = entry.second;
        for (const EdgeId &in_edge: graph.IncomingEdges(vertex)) {
            const auto &in_edge_name = id_mapper[in_edge.int_id()];
            EXPECT_EQ(edge_to_end.at(in_edge_name), name);
            ++checked_incidence_links;
        }
        for (const EdgeId &out_edge: graph.OutgoingEdges(vertex)) {
            const auto &out_edge_name = id_mapper[out_edge.int_id()];
            EXPECT_EQ(edge_to_start.at(out_edge_name), name);
            ++checked_incidence_links;
        }
    }
    EXPECT_EQ(checked_incidence_links, num_incidence_links);
    return vid_to_name;
}

void CheckLinks(std::ifstream &graph_stream, const Graph &graph, const IdMapper &id_mapper) {
    int num_links;
    graph_stream >> num_links;
    bool is_complex = false;
    for (int i = 0; i < num_links; ++i) {
        std::string first_edge_name, second_edge_name;
        int overlap;
        graph_stream >> first_edge_name >> second_edge_name >> overlap;
        EdgeId in_edge = id_mapper[first_edge_name];
        EdgeId out_edge = id_mapper[second_edge_name];
        EdgeId in_conjugate = graph.conjugate(in_edge);
        EdgeId out_conjugate = graph.conjugate(out_edge);
        VertexId vertex = graph.EdgeEnd(in_edge);
        EXPECT_EQ(vertex, graph.EdgeStart(out_edge));
        EXPECT_EQ(graph.conjugate(vertex), graph.EdgeEnd(out_conjugate));
        if (graph.is_complex(vertex)) {
            bool link_found = false;
            bool conj_link_found = false;
            for (const auto &link: graph.links(vertex)) {
                auto link_in_edge = graph.link(link).link.first;
                auto link_out_edge = graph.link(link).link.second;
                if (link_in_edge == in_edge and link_out_edge == out_edge) {
                    link_found = true;
                }
            }
            for (const auto &link: graph.links(graph.conjugate(vertex))) {
                auto link_in_edge = graph.link(link).link.first;
                auto link_out_edge = graph.link(link).link.second;
                if (link_in_edge == out_conjugate and link_out_edge == in_conjugate) {
                    conj_link_found = true;
                }
            }
            EXPECT_TRUE(link_found && conj_link_found);
            is_complex = true;
        }
    }
    if (is_complex) {
        EXPECT_EQ(num_links * 2, graph.link_size());
    }
}

void PerformSplits(debruijn_graph::Graph &graph, std::ifstream &ops_stream, const IdMapper &id_mapper) {
    int num_splits = 0;
    ops_stream >> num_splits;
    auto helper = graph.GetConstructionHelper();
    for (int i = 0; i < num_splits; ++i) {
        std::string op_name, in_edge_name, out_edge_name;
        ops_stream >> op_name >> in_edge_name >> out_edge_name;
        EdgeId in_edge = id_mapper[in_edge_name];
        EdgeId out_edge = id_mapper[out_edge_name];
        VertexId vertex = graph.EdgeEnd(in_edge);
        VERIFY_MSG(graph.EdgeStart(out_edge) == vertex, "Edges " << in_edge_name << " and " << out_edge_name <<
                                                                 " are not incident, operations are not consistent with the graph");
        VertexId new_vertex;
        if (graph.is_complex(vertex)) {
            LinkId split_link, conj_split_link;
            bool link_found = false;
            bool conj_link_found = false;
            for (auto &link_id: graph.links(vertex)) {
                if (graph.link(link_id).link.first == in_edge && graph.link(link_id).link.second == out_edge) {
                    split_link = link_id;
                    link_found = true;
                }
            }
            for (auto &link_id: graph.links(graph.conjugate(vertex))) {
                if (graph.link(link_id).link.first == graph.conjugate(out_edge) &&
                    graph.link(link_id).link.second == graph.conjugate(in_edge)) {
                    conj_split_link = link_id;
                    conj_link_found = true;
                }
            }

            EXPECT_TRUE(link_found && conj_link_found);
            std::vector<debruijn_graph::LinkId> empty;
            new_vertex = helper.CreateVertex(debruijn_graph::DeBruijnVertexData(empty));
            graph.add_link(new_vertex, split_link);
            graph.add_link(graph.conjugate(new_vertex), conj_split_link);

            graph.erase_links_with_outedge(vertex, out_edge);
            graph.erase_links_with_inedge(vertex, in_edge);
        } else {
            auto overlap = static_cast<uint>(graph.link_length(vertex, in_edge, out_edge));
            new_vertex = helper.CreateVertex(debruijn_graph::DeBruijnVertexData(overlap));
        }
        helper.DeleteLink(vertex, out_edge);
        helper.DeleteLink(graph.conjugate(vertex), graph.conjugate(in_edge));

        helper.LinkIncomingEdge(new_vertex, in_edge);
        helper.LinkOutgoingEdge(new_vertex, out_edge);
    }
}

void CheckPathOperations(debruijn_graph::Graph &graph,
                         const std::string &operations_path,
                         const std::string &fasta_path,
                         const IdMapper &id_mapper) {
    std::ifstream ops_stream(operations_path);
    INFO("Performing splits");
    PerformSplits(graph, ops_stream, id_mapper);
    int num_merges = 0;
    ops_stream >> num_merges;
    INFO("Performed splits");

    std::unordered_map<std::string, EdgeId> merged_id_map;
    for (int i = 0; i < num_merges; ++i) {
        std::string ctg_name;
        int num_edges;
        ops_stream >> ctg_name >> num_edges;
        std::vector<EdgeId> edges;
        for (int j = 0; j < num_edges; ++j) {
            std::string edge_name;
            ops_stream >> edge_name;
            edges.push_back(id_mapper[edge_name]);
        }
        std::vector<uint32_t> overlaps;
        for (size_t j = 1; j < edges.size(); ++j) {
            size_t overlap = graph.link_length(graph.EdgeEnd(edges[j - 1]), edges[j - 1], edges[j]);
            overlaps.push_back(static_cast<uint32_t>(overlap));
        }
        EdgeId resulting_edge = graph.MergePath(edges, true, overlaps);
        merged_id_map[ctg_name] = resulting_edge;
    }

    std::unordered_map<std::string, Sequence> ctg_name_to_seq;
    io::FileReadStream path_contigs(fasta_path);
    io::SingleRead contig;
    while (!path_contigs.eof()) {
        path_contigs >> contig;
        auto merge_result = merged_id_map.find(contig.name());
        EXPECT_NE(merge_result, merged_id_map.end());
        EXPECT_EQ(graph.EdgeNucls(merge_result->second), contig.sequence());
    }
}

void CheckGraphWithPaths(const std::filesystem::path &graph_basename) {
    std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());
    auto gfa_path = graph_basename;
    auto graph_path = graph_basename;
    auto fasta_path = graph_basename;
    auto operations_path = graph_basename;
    gfa_path += ".gfa";
    graph_path += ".graph";
    fasta_path += ".fasta";
    operations_path += ".ops";

    gfa::GFAReader gfa_reader(gfa_path);
    Graph graph(0);
    gfa_reader.to_graph(graph, id_mapper.get());

    omnigraph::OverlapHandler<Graph> overlap_handler(graph);

    std::ifstream graph_stream(graph_path.c_str());
    std::string graph_name;
    graph_stream >> graph_name;
    CheckSegmentLen(graph_stream, graph, *id_mapper);
    auto vid_to_name = CheckStructure(graph_stream, graph, *id_mapper);
    CheckLinks(graph_stream, graph, *id_mapper);

    CheckPathOperations(graph, operations_path, fasta_path, *id_mapper);
}

void CheckBinaryIO(const std::string &path_to_save,
                   const std::filesystem::path &graph_basename) {
    auto gfa_path = graph_basename;
    auto graph_path = graph_basename;
    gfa_path += ".gfa";
    graph_path += ".graph";

    size_t K = 55;
    std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());
    gfa::GFAReader gfa_reader(gfa_path);
    Graph graph(0);
    gfa_reader.to_graph(graph, id_mapper.get());
    size_t local_k = gfa_reader.to_graph(graph, id_mapper.get());

    io::binary::Save(path_to_save, graph);
    io::binary::Load(path_to_save, graph);

    std::ifstream graph_stream(graph_path.c_str());
    std::string graph_name;
    graph_stream >> graph_name;
    CheckSegmentLen(graph_stream, graph, *id_mapper);
    CheckStructure(graph_stream, graph, *id_mapper);
}

TEST(VariableOverlaps, BasicOperations) {
    CheckGraphWithPaths("src/test/debruijn/graph_fragments/v_overlaps/bone");
    CheckGraphWithPaths("src/test/debruijn/graph_fragments/v_overlaps/conjugate_bone");
    CheckGraphWithPaths("src/test/debruijn/graph_fragments/v_overlaps/conjugate_triple");
    CheckGraphWithPaths("src/test/debruijn/graph_fragments/v_overlaps/triple_repeat");
}

TEST(VariableOverlaps, BinaryIO) {
    std::string save_path = "src/test/debruijn/graph_fragments/saves/test_save";
    CheckBinaryIO(save_path, "src/test/debruijn/graph_fragments/v_overlaps/bone");
    CheckBinaryIO(save_path, "src/test/debruijn/graph_fragments/v_overlaps/conjugate_triple");
    CheckBinaryIO(save_path, "src/test/debruijn/graph_fragments/v_overlaps/triple_repeat");
}