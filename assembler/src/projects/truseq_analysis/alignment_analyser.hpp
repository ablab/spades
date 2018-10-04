//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/graph_pack.hpp"
#include "consistent_mapping.h"

namespace alignment_analysis {

    class AlignmentAnalyser {
    private:
        typedef debruijn_graph::conj_graph_pack::graph_t Graph;
        typedef Graph::EdgeId EdgeId;
        typedef Graph::VertexId VertexId;
        typedef debruijn_graph::BasicSequenceMapper<Graph, debruijn_graph::conj_graph_pack::index_t> Mapper;
        std::stringstream log_;
        const Graph &graph_;
        const Mapper &mapper_;
        const std::vector<io::SingleRead> &scaffolds_;
        const std::vector<io::SingleRead> &genome_;

        std::string str(const EdgeRange &er) const;

    public:
        AlignmentAnalyser(const std::vector<io::SingleRead> &scaffolds, const std::vector<io::SingleRead> &genome,
                          const Graph &graph, const Mapper &mapper);

        std::string Analyse(const io::SingleRead &genome_part);

    private:
        std::vector<ConsistentMapping> ExtractConsistentMappings(const MappingPath<EdgeId> &path);

        std::vector<ConsistentMapping> DetectAndMaskShortMutations(const std::vector<ConsistentMapping> &vector);
    };
}
