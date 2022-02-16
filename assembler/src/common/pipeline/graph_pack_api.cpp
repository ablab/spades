//***************************************************************************
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graph_pack.hpp"
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/graph_support/genomic_quality.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "assembly_graph/paths/bidirectional_path_container.hpp"
#include "modules/alignment/edge_index.hpp"
#include "modules/alignment/kmer_mapper.hpp"
#include "modules/alignment/long_read_storage.hpp"
#include "paired_info/paired_info.hpp"
#include "sequence/genome_storage.hpp"
#include "visualization/position_filler.hpp"
#include "genomic_info.hpp"

namespace graph_pack {

using namespace debruijn_graph;

void FillQuality(GraphPack& gp) {
    const auto &index = gp.get<EdgeIndex<Graph>>();
    const auto &kmer_mapper = gp.get<KmerMapper<Graph>>();
    const auto &genome = gp.get<GenomeStorage>();
    gp.get_mutable<EdgeQuality<Graph>>().Fill(index, kmer_mapper, genome.GetSequence());
}

//todo remove with usages after checking
void ClearQuality(GraphPack& gp) {
    gp.get_mutable<EdgeQuality<Graph>>().clear();
}

void EnsureIndex(GraphPack& gp) {
    auto &index = gp.get_mutable<EdgeIndex<Graph>>();
    if (index.IsAttached())
        return;

    INFO("Index refill");
    index.Refill();
    index.Attach();
}

void EnsureBasicMapping(GraphPack& gp) {
    auto &kmer_mapper = gp.get_mutable<KmerMapper<Graph>>();

    VERIFY(kmer_mapper.IsAttached());
    EnsureIndex(gp);
    INFO("Normalizing k-mer map. Total " << kmer_mapper.size() << " kmers to process");
    kmer_mapper.Normalize();
    INFO("Normalizing done");
}

void EnsureQuality(GraphPack& gp) {
    auto &edge_qual = gp.get_mutable<EdgeQuality<Graph>>();

    if (edge_qual.IsAttached())
        return;

    ClearQuality(gp);
    FillQuality(gp);
    edge_qual.Attach();
}

void EnsurePos(GraphPack& gp) {
    auto &edge_pos = gp.get_mutable<EdgesPositionHandler<Graph>>();

    if (!edge_pos.IsAttached())
        edge_pos.Attach();

    // Positions are refilled every time
    edge_pos.clear();

    const auto &genome = gp.get<GenomeStorage>();
    visualization::position_filler::FillPos(gp, genome.str(), "ref0");
    visualization::position_filler::FillPos(gp, ReverseComplement(genome.str()), "ref1");
}

void EnsureDebugInfo(GraphPack& gp) {
    EnsureBasicMapping(gp);
    EnsureQuality(gp);
    EnsurePos(gp);
}

void InitRRIndices(GraphPack& gp) {
    using Indices = omnigraph::de::PairedInfoIndicesT<Graph>;

    gp.get_mutable<Indices>("clustered_indices").Init();
    gp.get_mutable<Indices>("scaffolding_indices").Init();
}

void ClearRRIndices(GraphPack& gp) {
    using UnclusteredIndices = omnigraph::de::UnclusteredPairedInfoIndicesT<Graph>;
    using Indices = omnigraph::de::PairedInfoIndicesT<Graph>;

    gp.get_mutable<UnclusteredIndices>().Clear();
    gp.get_mutable<Indices>("clustered_indices").Clear();
    gp.get_mutable<Indices>("scaffolding_indices").Clear();

    gp.get_mutable<LongReadContainer<Graph>>().Clear();
}

void ClearPaths(GraphPack& gp) {
    gp.get_mutable<path_extend::PathContainer>("exSPAnder paths").clear();
}

void DetachAll(GraphPack& gp) {
    gp.DetachAll();
}

void PrepareForStage(GraphPack& gp, const char*) {
    gp.get_mutable<Graph>().clear_state();
}


void DetachEdgeIndex(GraphPack& gp) {
    gp.get_mutable<EdgeIndex<Graph>>().Detach();
}

} // namespace debruijn_graph
