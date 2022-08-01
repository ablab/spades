//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "cloud_only_links.hpp"

namespace cont_index {

void CloudOnlyLinksConstructor::CompareLinks(const scaffold_graph::ScaffoldGraph &hifi_graph,
                                             const scaffold_graph::ScaffoldGraph &tellseq_graph,
                                             cont_index::LinkIndexGraphConstructor::BarcodeScoreFunctionPtr score_function,
                                             const double &score_threshold,
                                             const string &output_path) {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    std::set<std::pair<ScaffoldVertex, ScaffoldVertex>> hifi_links;
    std::map<std::pair<ScaffoldVertex, ScaffoldVertex>, double> hifi_score_map;
    INFO(hifi_graph.EdgeCount() << " edges in hifi graph");
    INFO(tellseq_graph.EdgeCount() << " edges in tellseq graph");
    for (const auto &edge: hifi_graph.edges()) {
        std::pair<ScaffoldVertex, ScaffoldVertex> vertex_pair(edge.getStart(), edge.getEnd());
        hifi_links.emplace(edge.getStart(), edge.getEnd());
        hifi_score_map.emplace(vertex_pair, edge.getWeight());
    }
    std::vector<ScoreEntry> score_entries;
    if (tellseq_graph.VertexCount() != 0) {
        size_t new_tellseq_links = 0;
        for (const auto &edge: tellseq_graph.edges()) {
            if (math::ge(edge.getWeight(), score_threshold)) {
                std::pair<ScaffoldVertex, ScaffoldVertex> link(edge.getStart(), edge.getEnd());
                if (hifi_links.find(link) == hifi_links.end()) {
                    new_tellseq_links++;
                } else {
                    score_entries.emplace_back(edge.getStart(), edge.getEnd(), hifi_score_map.at(link), edge.getWeight());
                }
            }
        }
        INFO("Found " << new_tellseq_links << " new tellseq links");
    } else {
        for (const auto &edge: hifi_graph.edges()) {
            double tellseq_score = score_function->GetScore(edge);
            if (math::ge(tellseq_score, score_threshold)) {
                score_entries.emplace_back(edge.getStart(), edge.getEnd(), edge.getWeight(), tellseq_score);
            }
        }
    }
    INFO(score_entries.size() << " common links in hifi and tellseq graphs");
    std::ofstream os(output_path);
    os << "First\tSecond\tHiFi links\tTellSeq links" << "\n";
    for (const auto &entry: score_entries) {
        os << (*id_mapper_)[entry.first_.int_id()] << "\t" << (*id_mapper_)[entry.second_.int_id()]
           << "\t" << entry.hifi_score_ << "\t" << entry.tellseq_score_ << "\n";
    }
}


void CloudOnlyLinksConstructor::NormalizeTellseqLinks(const scaffold_graph::ScaffoldGraph &tellseq_graph,
                                                      const std::filesystem::path &output_path) const {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    const auto &assembly_graph = tellseq_graph.AssemblyGraph();
    std::map<std::pair<ScaffoldVertex, ScaffoldVertex>, double> score_map;
    double score_threshold = 2.0;
    if (tellseq_graph.VertexCount() != 0) {
        for (const auto &edge: tellseq_graph.edges()) {
            const auto &first_vertex = edge.getStart();
            const auto &first_conj = first_vertex.GetConjugateFromGraph(assembly_graph);
            const auto &second_vertex = edge.getEnd();
            const auto &second_conj = second_vertex.GetConjugateFromGraph(assembly_graph);
            size_t first_len = first_vertex.GetLengthFromGraph(assembly_graph);
            size_t second_len = second_vertex.GetLengthFromGraph(assembly_graph);
            double num_barcodes = edge.getWeight();
            bool is_link = first_len >= length_threshold_
                           and second_len >= length_threshold_
                           and math::ge(num_barcodes, score_threshold);
            if (is_link) {
                std::pair<ScaffoldVertex, ScaffoldVertex> link(first_vertex, second_vertex);
                std::pair<ScaffoldVertex, ScaffoldVertex> rc_rc_link(second_conj, first_conj);
                if (first_vertex != second_vertex and first_vertex != second_conj) {
                    score_map[link] += num_barcodes / 2;
                    score_map[rc_rc_link] += num_barcodes / 2;
                }
            }
        }
    }
    INFO(score_map.size() << " filtered tellseq links");
    std::unordered_map<ScaffoldVertex, size_t> vertex_to_head_barcodes;
    for (const auto &vertex: tellseq_graph.vertices()) {
        vertex_to_head_barcodes[vertex] = barcode_extractor_ptr_->GetBarcodesFromHead(vertex.GetFirstEdge(),
                                                                                      count_threshold_,
                                                                                      tail_threshold_).size();
    }
    std::ofstream os(output_path / "normalized_tellseq_links.tsv");
    os << "First\tSecond\tTotal links\tJaccard Index" << "\n";
    for (const auto &entry: score_map) {
        const auto &first_vertex = entry.first.first;
        const auto &first_conj = first_vertex.GetConjugateFromGraph(assembly_graph);
        const auto &second_vertex = entry.first.second;
        double shared_barcodes = entry.second;
        auto first_tail_barcodes = vertex_to_head_barcodes[first_conj];
        auto second_head_barcodes = vertex_to_head_barcodes[second_vertex];
        auto union_size = static_cast<double>(first_tail_barcodes + second_head_barcodes) - shared_barcodes;
        double jaccard_index = shared_barcodes / union_size;
        os << (*id_mapper_)[first_vertex.int_id()] << "\t" << (*id_mapper_)[second_vertex.int_id()]
           << "\t" << shared_barcodes << "\t" << jaccard_index << "\n";
    }
}
void CloudOnlyLinksConstructor::ConstructCloudOnlyLinks(const debruijn_graph::Graph &graph,
                                                        bool bin_load,
                                                        bool debug,
                                                        const std::filesystem::path &output_dir) const {
    auto tellseq_graph = cont_index::GetTellSeqScaffoldGraph(graph, barcode_extractor_ptr_, graph_score_threshold_,
                                                             length_threshold_, tail_threshold_, count_threshold_,
                                                             max_threads_, bin_load, debug, output_dir, id_mapper_);

    LinkIndexGraphConstructor link_index_constructor(graph, barcode_extractor_ptr_, graph_score_threshold_,
                                                     tail_threshold_, length_threshold_, count_threshold_, max_threads_);
    auto score_function = link_index_constructor.ConstructScoreFunction();
    NormalizeTellseqLinks(tellseq_graph, output_dir);
    auto compare_output_path = output_dir / "hifi_tellseq_scores.tsv";
//    CompareLinks(hifi_graph, tellseq_graph, score_function, id_mapper.get(), compare_output_path);
}
}
