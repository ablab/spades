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

class ChromosomeInfo {
    std::string name_;
    std::vector<EdgeId> path_;
    std::multimap<EdgeId, size_t> edge_idxs_;

public:
    ChromosomeInfo() {}

    explicit ChromosomeInfo(const string &name, const MappingPath<EdgeId> &mapping_path) :
            name_(name),
            path_(mapping_path.simple_path()) {
        for (size_t i = 0; i < path_.size(); ++i) {
            edge_idxs_.insert(std::make_pair(path_[i], i));
        }
    }

    size_t Multiplicity(EdgeId e) const {
        return edge_idxs_.count(e);
    }

    size_t IsUnique(EdgeId e) const {
        return Multiplicity(e) == 1;
    }

    EdgeId EdgeAt(size_t idx) const {
        VERIFY(idx < path_.size());
        return path_[idx];
    }

    vector<size_t> EdgeIdxs(EdgeId e) const {
        return get_all(edge_idxs_, e);
    }

    size_t UniqueEdgeIdx(EdgeId e) const {
        vector<size_t> idxs = EdgeIdxs(e);
        VERIFY(idxs.size() == 1);
        return idxs.front();
    }

    const std::string& name() const {
        return name_;
    }

    size_t size() const {
        return path_.size();
    }
};

class GenomeInfo {
    std::map<string, ChromosomeInfo> chr_infos_;
public:
    void AddInfo(ChromosomeInfo info) {
        VERIFY(!chr_infos_.count(info.name()));
        chr_infos_[info.name()] = std::move(info);
    }

    const ChromosomeInfo& ChrInfo(const string &name) const {
        return get(chr_infos_, name);
    }

    vector<string> ChromosomesByEdge(EdgeId e) const {
        vector<string> answer;
        for (const auto& chr_info: chr_infos_)
            if (chr_info.second.Multiplicity(e))
                answer.push_back(chr_info.first);
        return answer;
    }

    size_t Multiplicity(EdgeId e) const {
        size_t ans = 0;
        for (const auto& chr_info: chr_infos_)
            ans += chr_info.second.Multiplicity(e);
        return ans;
    }

    bool IsUnique(EdgeId e) const {
        return Multiplicity(e) == 1;
    }

    bool InUniqueChromosome(EdgeId e) const {
        return ChromosomesByEdge(e).size() == 1;
    }

    const ChromosomeInfo& UniqueChromosomeInfo(EdgeId e) const {
        auto chr_names = ChromosomesByEdge(e);
        VERIFY(chr_names.size() == 1);
        return ChrInfo(chr_names.front());
    }

    pair<string, size_t> UniqueChromosomeIdx(EdgeId e) const {
        VERIFY(IsUnique(e));
        auto chrs = ChromosomesByEdge(e);
        VERIFY(chrs.size() == 1);
        return std::make_pair(chrs.front(), ChrInfo(chrs.front()).UniqueEdgeIdx(e));
    }
};

class GenomeConsistenceChecker {
    const conj_graph_pack &gp_;
    //EdgesPositionHandler<Graph> &position_handler_;
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
    GenomeInfo genome_info_;
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
    vector<size_t> SpellChromosome(const std::string &chr, const MappingPath<EdgeId> &mapping_path);

};


}
