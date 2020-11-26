//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graph_pack.hpp"
#include "assembly_graph/components/connected_component.hpp"
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"
#include "assembly_graph/graph_support/genomic_quality.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "assembly_graph/paths/bidirectional_path_container.hpp"
#include "common/modules/alignment/rna/ss_coverage.hpp"
#include "modules/alignment/edge_index.hpp"
#include "modules/alignment/kmer_mapper.hpp"
#include "modules/alignment/long_read_storage.hpp"
#include "paired_info/paired_info.hpp"
#include "sequence/genome_storage.hpp"
#include "visualization/position_filler.hpp"
#include "genomic_info.hpp"

namespace debruijn_graph {

GraphPack::GraphPack(size_t k, const std::string &workdir, size_t lib_count,
                     const std::vector<std::string> &genome,
                     size_t flanking_range, size_t max_mapping_gap, size_t max_gap_diff,
                     bool detach_indices) : k_(k), workdir_(workdir) {
    using namespace omnigraph::de;
    Graph &g = emplace<Graph>(k);
    emplace<EdgeIndex<Graph>>(g, workdir);
    emplace<KmerMapper<Graph>>(g);
    emplace<FlankingCoverage<Graph>>(g, flanking_range);
    emplace<UnclusteredPairedInfoIndicesT<Graph>>(g, lib_count);
    emplace_with_key<PairedInfoIndicesT<Graph>>("clustered_indices", g, lib_count);
    emplace_with_key<PairedInfoIndicesT<Graph>>("scaffolding_indices", g, lib_count);
    emplace<LongReadContainer<Graph>>(g, lib_count);
    emplace<path_extend::TrustedPathsContainer>(lib_count);
    emplace<SSCoverageContainer>(g, lib_count);
    emplace<GenomicInfo>();
    emplace<GenomeStorage>(genome);
    emplace<EdgeQuality<Graph>>(g);
    emplace<EdgesPositionHandler<Graph>>(g, max_mapping_gap + k, max_gap_diff);
    emplace<ConnectedComponentCounter>(g);
    emplace_with_key<path_extend::PathContainer>("exSPAnder paths");
    if (detach_indices)
        DetachAll();
}

void GraphPack::FillQuality() {
    const auto &index = get<EdgeIndex<Graph>>();
    const auto &kmer_mapper = get<KmerMapper<Graph>>();
    const auto &genome = get<GenomeStorage>();
    get_mutable<EdgeQuality<Graph>>().Fill(index, kmer_mapper, genome.GetSequence());
}

//todo remove with usages after checking
void GraphPack::ClearQuality() {
    get_mutable<EdgeQuality<Graph>>().clear();
}

void GraphPack::EnsureIndex() {
    auto &index = get_mutable<EdgeIndex<Graph>>();
    if (index.IsAttached())
        return;

    INFO("Index refill");
    index.Refill();
    index.Attach();
}

void GraphPack::EnsureBasicMapping() {
    auto &kmer_mapper = get_mutable<KmerMapper<Graph>>();

    VERIFY(kmer_mapper.IsAttached());
    EnsureIndex();
    INFO("Normalizing k-mer map. Total " << kmer_mapper.size() << " kmers to process");
    kmer_mapper.Normalize();
    INFO("Normalizing done");
}

void GraphPack::EnsureQuality() {
    auto &edge_qual = get_mutable<EdgeQuality<Graph>>();

    if (edge_qual.IsAttached())
        return;

    ClearQuality();
    FillQuality();
    edge_qual.Attach();
}

void GraphPack::EnsurePos() {
    auto &edge_pos = get_mutable<EdgesPositionHandler<Graph>>();

    if (!edge_pos.IsAttached())
        edge_pos.Attach();

    // Positions are refilled every time
    edge_pos.clear();

    const auto &genome = get<GenomeStorage>();
    visualization::position_filler::FillPos(*this, genome.str(), "ref0");
    visualization::position_filler::FillPos(*this, ReverseComplement(genome.str()), "ref1");
}

void GraphPack::EnsureDebugInfo() {
    EnsureBasicMapping();
    EnsureQuality();
    EnsurePos();
}

void GraphPack::InitRRIndices() {
    using Indices = omnigraph::de::PairedInfoIndicesT<Graph>;

    get_mutable<Indices>("clustered_indices").Init();
    get_mutable<Indices>("scaffolding_indices").Init();
}

void GraphPack::ClearRRIndices() {
    using UnclusteredIndices = omnigraph::de::UnclusteredPairedInfoIndicesT<Graph>;
    using Indices = omnigraph::de::PairedInfoIndicesT<Graph>;

    get_mutable<UnclusteredIndices>().Clear();
    get_mutable<Indices>("clustered_indices").Clear();
    get_mutable<Indices>("scaffolding_indices").Clear();

    get_mutable<LongReadContainer<Graph>>().Clear();
}

void GraphPack::ClearPaths() {
    get_mutable<path_extend::PathContainer>("exSPAnder paths").clear();
}

void GraphPack::DetachAll() {
    get_mutable<EdgeIndex<Graph>>().Detach();
    get_mutable<KmerMapper<Graph>>().Detach();
    get_mutable<EdgesPositionHandler<Graph>>().Detach();
    get_mutable<EdgeQuality<Graph>>().Detach();
}

void GraphPack::PrepareForStage(const char*) {
    get_mutable<Graph>().clear_state();
}


void GraphPack::DetachEdgeIndex() {
    get_mutable<EdgeIndex<Graph>>().Detach();
}

} // namespace debruijn_graph
