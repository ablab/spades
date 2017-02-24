//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/standard_base.hpp"
#include "pipeline/graph_pack.hpp"
#include "consistent_mapping.h"

namespace alignment_analysis {

    class AlignmentAnalyser {
    private:
        typedef debruijn_graph::conj_graph_pack::graph_t Graph;
        typedef Graph::EdgeId EdgeId;
        typedef Graph::VertexId VertexId;
        typedef debruijn_graph::BasicSequenceMapper<Graph, debruijn_graph::conj_graph_pack::index_t> Mapper;
        stringstream log_;
        const Graph &graph_;
        const Mapper &mapper_;
        const vector <io::SingleRead> &scaffolds_;
        const vector <io::SingleRead> &genome_;

        string str(const EdgeRange &er) const;

    public:
        AlignmentAnalyser(const vector <io::SingleRead> &scaffolds, const vector <io::SingleRead> &genome,
                          const Graph &graph, const Mapper &mapper);

        string Analyse(const io::SingleRead &genome_part);

    private:
        vector <ConsistentMapping> ExtractConsistentMappings(const MappingPath<EdgeId> &path);

        vector <ConsistentMapping> DetectAndMaskShortMutations(const vector <ConsistentMapping> &vector);
    };
}
