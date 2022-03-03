//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graph_pack.hpp"

#include "genomic_info.hpp"

#include "alignment/edge_index.hpp"
#include "alignment/kmer_mapper.hpp"
#include "alignment/long_read_storage.hpp"
#include "alignment/rna/ss_coverage.hpp"
#include "assembly_graph/components/connected_component.hpp"
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"
#include "assembly_graph/graph_support/genomic_quality.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "assembly_graph/paths/bidirectional_path_container.hpp"
#include "paired_info/paired_info.hpp"
#include "sequence/genome_storage.hpp"
#include "visualization/position_filler.hpp"

namespace graph_pack {

using namespace debruijn_graph;

GraphPack::GraphPack(size_t k, const std::filesystem::path &workdir, size_t lib_count,
                     const std::vector<std::string> &genome,
                     size_t flanking_range, size_t max_mapping_gap, size_t max_gap_diff,
                     bool detach_indices) : k_(k), workdir_(workdir) {
    using namespace omnigraph::de;
    Graph &g = emplace<Graph>(k);
    emplace<EdgeIndex<Graph>>(g, workdir);
    emplace<KmerMapper<Graph>>(g);
    emplace<omnigraph::FlankingCoverage<Graph>>(g, flanking_range);
    emplace<UnclusteredPairedInfoIndicesT<Graph>>(g, lib_count);
    emplace_with_key<PairedInfoIndicesT<Graph>>("clustered_indices", g, lib_count);
    emplace_with_key<PairedInfoIndicesT<Graph>>("scaffolding_indices", g, lib_count);
    emplace<LongReadContainer<Graph>>(g, lib_count);
    emplace<path_extend::TrustedPathsContainer>(lib_count);
    emplace<SSCoverageContainer>(g, lib_count);
    emplace<GenomicInfo>();
    emplace<GenomeStorage>(genome);
    emplace<EdgeQuality<Graph>>(g);
    emplace<omnigraph::EdgesPositionHandler<Graph>>(g, max_mapping_gap + k, max_gap_diff);
    emplace<ConnectedComponentCounter>(g);
    emplace_with_key<path_extend::PathContainer>("exSPAnder paths");
    if (detach_indices)
        DetachAll();
}

void GraphPack::DetachAll() {
    get_mutable<EdgeIndex<Graph>>().Detach();
    get_mutable<KmerMapper<Graph>>().Detach();
    get_mutable<omnigraph::EdgesPositionHandler<Graph>>().Detach();
    get_mutable<EdgeQuality<Graph>>().Detach();
}

} // namespace debruijn_graph
