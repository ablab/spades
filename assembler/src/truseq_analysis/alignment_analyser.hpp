#pragma once

#include <io/single_read.hpp>
#include <graph_pack.hpp>
#include "standard_base.hpp"
#include "consistent_mapping.h"

namespace alignment_analysis {

    class AlignmentAnalyser {
    private:
        typedef debruijn_graph::conj_graph_pack::graph_t Graph;
        typedef Graph::EdgeId EdgeId;
        typedef Graph::VertexId VertexId;
        typedef debruijn_graph::NewExtendedSequenceMapper<Graph, debruijn_graph::conj_graph_pack::index_t> Mapper;
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
