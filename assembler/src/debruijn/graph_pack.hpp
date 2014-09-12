//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "omni/id_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "de/paired_info.hpp"
#include "debruijn_graph.hpp"
#include "config_struct.hpp"
#include "edge_index.hpp"
#include "genomic_quality.hpp"
#include "sequence_mapper.hpp"
#include "genomic_info.hpp"
#include "long_read_storage.hpp"
#include "detail_coverage.hpp"

namespace debruijn_graph {

/*KmerFree*//*KmerStoring*/
template<class Graph, class SeqType, class KmerEdgeIndex = KmerStoringEdgeIndex<Graph, SeqType, kmer_index_traits<SeqType>, DefaultStoring>>
struct graph_pack: private boost::noncopyable {
    typedef Graph graph_t;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef SeqType seq_t;
    typedef EdgeIndex<graph_t, seq_t, KmerEdgeIndex> index_t;
    typedef omnigraph::de::PairedInfoIndicesT<Graph> PairedInfoIndicesT;
    typedef omnigraph::de::UnclusteredPairedInfoIndicesT<Graph> UnclusteredPairedInfoIndicesT;
    typedef LongReadContainer<Graph> LongReadContainerT;

    size_t k_value;

    graph_t g;
    index_t index;
    KmerMapper<graph_t, seq_t> kmer_mapper;
    FlankingCoverage<graph_t> flanking_cov;
    UnclusteredPairedInfoIndicesT paired_indices;
    PairedInfoIndicesT clustered_indices;
    PairedInfoIndicesT scaffolding_indices;
    LongReadContainerT single_long_reads;
    GenomicInfo ginfo;

    Sequence genome;
	EdgeQuality<Graph> edge_qual;
    EdgesPositionHandler<graph_t> edge_pos;
 
    graph_pack(size_t k, const std::string &workdir, size_t lib_count,
                        Sequence genome = Sequence(),
                        size_t flanking_range = 50,
                        size_t max_mapping_gap = 0,
                        size_t max_gap_diff = 0,
                        bool detach_indices = true)
            : k_value(k), g(k), index(g, workdir),
              kmer_mapper(g),
              flanking_cov(g, flanking_range),
              paired_indices(lib_count),
              clustered_indices(g, lib_count),
              scaffolding_indices(g, lib_count),
              single_long_reads(g, lib_count),
              genome(genome),
              edge_qual(g),
              edge_pos(g, max_mapping_gap + k, max_gap_diff)
    { 
        if (detach_indices) {
            DetachAll();
        }
    }

    void FillQuality() {
        edge_qual.Fill(index, kmer_mapper, genome);
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
        FillPos(*this, genome, "ref0");
        FillPos(*this, !genome, "ref1");
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

    void DetachAll() {
        index.Detach();
        kmer_mapper.Detach();
        edge_pos.Detach();
        edge_qual.Detach();
    }

};

typedef graph_pack<ConjugateDeBruijnGraph, runtime_k::RtSeq, KmerFreeEdgeIndex<Graph, runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq>, DefaultStoring>> conj_graph_pack;
typedef conj_graph_pack::index_t Index;

typedef conj_graph_pack::PairedInfoIndicesT PairedIndicesT;
typedef conj_graph_pack::UnclusteredPairedInfoIndicesT UnclusteredPairedIndicesT;
typedef conj_graph_pack::LongReadContainerT LongReadContainerT;
typedef omnigraph::de::PairedInfoIndexT<ConjugateDeBruijnGraph> PairedIndexT;
typedef omnigraph::de::UnclusteredPairedInfoIndexT<ConjugateDeBruijnGraph> UnclusteredPairedIndexT;


} // namespace debruijn_graph
