//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#pragma once
#include "visualization/graph_labeler.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "sequence/sequence.hpp"
#include "pipeline/graph_pack.hpp"
#include "visualization/position_filler.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"

namespace debruijn_graph {


using path_extend::BidirectionalPath;
using path_extend::ScaffoldingUniqueEdgeStorage;

struct PathScore{
    size_t misassemblies;
    size_t wrong_gap_size;
    size_t mapped_length;
    PathScore(size_t m, size_t w, size_t ml): misassemblies(m), wrong_gap_size(w), mapped_length(ml) {}
};
class GenomeConsistenceChecker {

private:
    const conj_graph_pack &gp_;
    //EdgesPositionHandler<Graph> &position_handler_;
    Sequence genome_;
    const ScaffoldingUniqueEdgeStorage &storage_;
    size_t absolute_max_gap_;
    double relative_max_gap_;
    set<EdgeId> excluded_unique_;
    EdgeId circular_edge_;
//map from unique edges to their order in genome spelling;
    mutable map<EdgeId, size_t> genome_spelled_;
    bool consequent(const Range &mr1, const Range &mr2) const;
    bool consequent(const MappingRange &mr1, const MappingRange &mr2) const ;

    PathScore CountMisassembliesWithStrand(const BidirectionalPath &path, const string strand) const;
//constructs longest sequence of consequetive ranges, stores result in used_mappings
    void FindBestRangeSequence(const set<MappingRange>& old_mappings, vector<MappingRange>& used_mappings) const;
//Refills genomic positions uniting alingments separated with small gaps
    void RefillPos();
    void RefillPos(const string &strand);
    void RefillPos(const string &strand, const EdgeId &e);
DECL_LOGGER("GenomeConsistenceChecker");


public:
    GenomeConsistenceChecker(const conj_graph_pack &gp, const ScaffoldingUniqueEdgeStorage &storage, size_t max_gap, double relative_max_gap /*= 0.2*/) : gp_(gp),
            genome_(gp.genome.GetSequence()), storage_(storage),
        absolute_max_gap_(max_gap), relative_max_gap_(relative_max_gap), excluded_unique_(), circular_edge_() {
        if (!gp.edge_pos.IsAttached()) {
            gp.edge_pos.Attach();
        }
        gp.edge_pos.clear();
        visualization::position_filler::FillPos(gp_, gp_.genome.GetSequence(), "0");
        visualization::position_filler::FillPos(gp_, !gp_.genome.GetSequence(), "1");
        RefillPos();
    }
    PathScore CountMisassemblies(const BidirectionalPath &path) const;
    vector<pair<EdgeId, MappingRange> > ConstructEdgeOrder() const;

//spells genome in language of long unique edges from storage;
    void SpellGenome();


};


}
