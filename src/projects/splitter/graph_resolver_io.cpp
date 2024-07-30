//***************************************************************************
//* Copyright (c) 2021-2023 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graph_resolver_io.hpp"

#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "io/graph/gfa_writer.hpp"

namespace cont_index {

void TransformedGraphIO::PrintGraph(const debruijn_graph::Graph &graph,
                                    const GraphResolver::GraphResolverInfo &resolver_info,
                                    const std::filesystem::path &output_base) const {
    auto resolved_graph_out = std::ofstream(output_base / ("resolved_graph.gfa"));
    path_extend::GFAPathWriter resolved_graph_writer(graph, resolved_graph_out);
    resolved_graph_writer.WriteSegmentsAndLinks();

    auto new_name_generator = std::make_shared<path_extend::DefaultContigNameGenerator>();
    path_extend::ContigWriter resolved_writer(graph, new_name_generator);
    path_extend::PathContainer resolved_edges;
    for (const auto &edge: graph.canonical_edges()) {
        resolved_edges.Create(graph, edge);
    }
    std::vector<path_extend::PathsWriterT> edge_writers;
    edge_writers.push_back([&](const path_extend::ScaffoldStorage &scaffold_storage) {
      auto fn = output_base / ("resolved_edges.fasta");
      INFO("Outputting edges to " << fn);
      path_extend::ContigWriter::WriteScaffolds(scaffold_storage, fn);
    });
    resolved_writer.OutputPaths(resolved_edges, edge_writers);

    auto edge_out = output_base / "edge_transform.tsv";
    auto edge_out_stream = std::ofstream(edge_out);
    edge_out_stream << "Original edge id\tResolved graph edge id\n";
    for (const auto &entry: resolver_info.original_edge_to_transformed) {
        edge_out_stream << (*id_mapper_)[entry.first.int_id()] << "\t" << entry.second.int_id() << std::endl;
    }
}
}
