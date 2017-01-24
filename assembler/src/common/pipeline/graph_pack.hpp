//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/indices/edge_position_index.hpp"
#include "utils/indices/storing_traits.hpp"
#include "sequence/genome_storage.hpp"
#include "assembly_graph/handlers/id_track_handler.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "assembly_graph/core/graph.hpp"
#include "paired_info/paired_info.hpp"
#include "pipeline/config_struct.hpp"
#include "modules/alignment/edge_index.hpp"
#include "assembly_graph/graph_support/genomic_quality.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "genomic_info.hpp"
#include "modules/alignment/long_read_storage.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"
#include "assembly_graph/components/connected_component.hpp"
#include "modules/alignment/kmer_mapper.hpp"
#include "common/visualization/position_filler.hpp"
#include "common/assembly_graph/paths/bidirectional_path.hpp"

namespace debruijn_graph {

template<class Graph>
struct graph_pack: private boost::noncopyable {
    typedef Graph graph_t;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef RtSeq seq_t;
    typedef EdgeIndex<graph_t> index_t;
    using PairedInfoIndicesT = omnigraph::de::PairedInfoIndicesT<Graph>;
    //typedef omnigraph::de::PairedInfoIndicesT<Graph> PairedInfoIndicesT;
    typedef omnigraph::de::UnclusteredPairedInfoIndicesT<Graph> UnclusteredPairedInfoIndicesT;
    typedef LongReadContainer<Graph> LongReadContainerT;

    size_t k_value;

    graph_t g;
    index_t index;
    KmerMapper<graph_t> kmer_mapper;
    FlankingCoverage<graph_t> flanking_cov;
    UnclusteredPairedInfoIndicesT paired_indices;
    PairedInfoIndicesT clustered_indices;
    PairedInfoIndicesT scaffolding_indices;
    LongReadContainerT single_long_reads;
    GenomicInfo ginfo;

    GenomeStorage genome;
    EdgeQuality<Graph> edge_qual;
    mutable EdgesPositionHandler<graph_t> edge_pos;
    ConnectedComponentCounter components;
    path_extend::PathContainer contig_paths;

    graph_pack(size_t k, const std::string &workdir, size_t lib_count,
                        const std::string &genome = "",
                        size_t flanking_range = 50,
                        size_t max_mapping_gap = 0,
                        size_t max_gap_diff = 0,
                        bool detach_indices = true)
            : k_value(k), g(k), index(g, workdir),
              kmer_mapper(g),
              flanking_cov(g, flanking_range),
              paired_indices(g, lib_count),
              clustered_indices(g, lib_count),
              scaffolding_indices(g, lib_count),
              single_long_reads(g, lib_count),
              genome(genome),
              edge_qual(g),
              edge_pos(g, max_mapping_gap + k, max_gap_diff),
              components(g),
              contig_paths()
    { 
        if (detach_indices) {
            DetachAll();
        }
    }

    void FillQuality() {
        edge_qual.Fill(index, kmer_mapper, genome.GetSequence());
    }

    //todo remove with usages after checking
    void ClearQuality() {
        edge_qual.clear();
    }

    void EnsureIndex() {
        if (!index.IsAttached()) {
            INFO("Index refill");
            index.Refill();
            index.Attach();
        }
    }

    void EnsureBasicMapping() {
        VERIFY(kmer_mapper.IsAttached());
        EnsureIndex();
        INFO("Normalizing k-mer map. Total " << kmer_mapper.size() << " kmers to process");
        kmer_mapper.Normalize();
        INFO("Normalizing done");
    }

    void EnsureQuality() {
        if (!edge_qual.IsAttached()) {
            ClearQuality();
            FillQuality();
            edge_qual.Attach();
        }
    }

    //positions are refilled every time
    void EnsurePos() {
        if (!edge_pos.IsAttached()) {
            edge_pos.Attach();
        }
        edge_pos.clear();
        visualization::position_filler::FillPos(*this, genome.GetSequence(), "ref0");
        visualization::position_filler::FillPos(*this, !genome.GetSequence(), "ref1");
    }
    
    void EnsureDebugInfo() {
        EnsureBasicMapping();
        EnsureQuality();
        EnsurePos();
    }

    void InitRRIndices() {
        clustered_indices.Init();
        scaffolding_indices.Init();
    }

    void ClearRRIndices() {
        for (auto& pi : paired_indices) {
            pi.clear();
        }
        clustered_indices.Clear();
        scaffolding_indices.Clear();
        single_long_reads.Clear();
    }

    void ClearPaths() {
        contig_paths.DeleteAllPaths();
    }

    void DetachAll() {
        index.Detach();
        kmer_mapper.Detach();
        edge_pos.Detach();
        edge_qual.Detach();
    }

};

typedef graph_pack<ConjugateDeBruijnGraph> conj_graph_pack;
typedef conj_graph_pack::index_t Index;

typedef conj_graph_pack::PairedInfoIndicesT PairedIndicesT;
typedef conj_graph_pack::UnclusteredPairedInfoIndicesT UnclusteredPairedIndicesT;
typedef conj_graph_pack::LongReadContainerT LongReadContainerT;
typedef omnigraph::de::PairedInfoIndexT<ConjugateDeBruijnGraph> PairedIndexT;
typedef omnigraph::de::UnclusteredPairedInfoIndexT<ConjugateDeBruijnGraph> UnclusteredPairedIndexT;

} // namespace debruijn_graph
