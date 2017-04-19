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
#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "modules/path_extend/pe_utils.hpp"

namespace debruijn_graph {

using path_extend::BidirectionalPath;
using path_extend::ScaffoldingUniqueEdgeStorage;
using omnigraph::MappingPath;

struct PathScore{
    size_t misassemblies;
    size_t wrong_gap_size;
    size_t mapped_length;
    PathScore(size_t m, size_t w, size_t ml): misassemblies(m), wrong_gap_size(w), mapped_length(ml) {}
};

class GenomeConsistenceChecker {
    const conj_graph_pack &gp_;
    //EdgesPositionHandler<Graph> &position_handler_;
    Sequence genome_;
    const ScaffoldingUniqueEdgeStorage &storage_;
    size_t absolute_max_gap_;
    double relative_max_gap_;
    set<EdgeId> excluded_unique_;
//Edges containing zero point for each reference
//TODO: do we need circular/linear chromosomes support?
    set<EdgeId> circular_edges_;
    size_t unresolvable_len_;
    const vector<shared_ptr<path_extend::GraphCoverageMap>> &long_reads_cov_map_;
//FIXME:: rename
//map from unique edges to their order in genome spelling;
    map<EdgeId, pair<std::string, size_t>> genome_spelled_;
    map<std::string, vector<EdgeId> > chromosomes_spelled_;
    bool consequent(const Range &mr1, const Range &mr2) const;
    bool consequent(const MappingRange &mr1, const MappingRange &mr2) const ;
    void PrintMisassemblyInfo(EdgeId e1, EdgeId e2) const;

    PathScore InternalCountMisassemblies(const BidirectionalPath &path) const;
//map from unique edges to their order in genome spelling;
// for each chromosome provides its spelling

//constructs longest sequence of consequetive ranges, stores result in used_mappings
    void FindBestRangeSequence(const set<MappingRange>& old_mappings, vector<MappingRange>& used_mappings) const;
//Refills genomic positions uniting alingments separated with small gaps
    void RefillPos();
    void RefillPos(const string &strand);
    void RefillPos(const string &strand, const EdgeId &e);
    void ReportEdge(EdgeId e, double w) const;
    void ReportVariants(std::vector<pair<double, EdgeId>> &sorted_w) const;
    size_t GetSupportingPathCount(EdgeId e1, EdgeId e2, size_t lib_index) const;
DECL_LOGGER("GenomeConsistenceChecker");


public:
    GenomeConsistenceChecker(const conj_graph_pack &gp, const ScaffoldingUniqueEdgeStorage &storage, size_t max_gap,
                             double relative_max_gap /*= 0.2*/, size_t unresolvable_len,
                             const vector<shared_ptr<path_extend::GraphCoverageMap>> &long_reads_cov_map) :
            gp_(gp),
            storage_(storage), absolute_max_gap_(max_gap),
            relative_max_gap_(relative_max_gap),
            unresolvable_len_(unresolvable_len), long_reads_cov_map_(long_reads_cov_map) {
        if (!gp.edge_pos.IsAttached()) {
            gp.edge_pos.Attach();
        }
        gp.edge_pos.clear();
        auto chromosomes = gp_.genome.GetChromosomes();
        for (auto chr: chromosomes) {
            auto seq = Sequence(chr.sequence);
            visualization::position_filler::FillPos(gp_, seq, "_0_" + chr.name);
            visualization::position_filler::FillPos(gp_, !seq, "_1_" + chr.name);
        }
        RefillPos();
    }

    PathScore CountMisassemblies(const BidirectionalPath &path) const;
    void CheckPathEnd(const BidirectionalPath &path) const;
    map<std::string, MappingPath<EdgeId>> ConstructEdgeOrder() const;
    MappingPath<EdgeId> ConstructEdgeOrder(const std::string &chr_name) const;

//spells genome in language of long unique edges from storage;
    void SpellGenome();
    void SpellChromosome(const std::string &chr, const MappingPath<EdgeId> &mapping_path, vector<size_t>& lengths);

};


}
