//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/alignment/sequence_mapper.hpp"
#include "pipeline/graph_pack.hpp"
#include "consistent_mapping.h"

namespace alignment_analysis {

    class AlignmentAnalyser {
    private:
        using Mapper = debruijn_graph::BasicSequenceMapper<Graph, debruijn_graph::EdgeIndex<Graph>>;
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
        std::vector<ConsistentMapping> ExtractConsistentMappings(const debruijn_graph::MappingPath<EdgeId> &path);

        std::vector<ConsistentMapping> DetectAndMaskShortMutations(const std::vector<ConsistentMapping> &vector);
    };
}
