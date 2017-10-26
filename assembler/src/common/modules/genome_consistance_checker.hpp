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
        return utils::get_all(edge_idxs_, e);
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
        return utils::get(chr_infos_, name);
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

    vector<string> Chromosomes() const {
        vector<string> answer;
        utils::push_back_all(answer, utils::key_set(chr_infos_));
        return answer;
    }
};

class GenomeConsistenceChecker {
    typedef omnigraph::MappingPath<EdgeId> MappingPathT;
    const conj_graph_pack &gp_;
    const size_t absolute_max_gap_;
    const double relative_max_gap_;
    const size_t unresolvable_len_;

    const ScaffoldingUniqueEdgeStorage &storage_;
    const vector<path_extend::GraphCoverageMap> &long_reads_cov_map_;
    static const size_t SIGNIFICANT_LENGTH_LOWER_LIMIT = 10000;
    GenomeInfo genome_info_;
    //Edges containing zero point for each reference
    //TODO: do we need circular/linear chromosomes support?
    set<EdgeId> circular_edges_;

    io::DataSet<config::LibraryData> reads_;
    bool Consequent(const MappingRange &mr1, const MappingRange &mr2) const ;

    void PrintMisassemblyInfo(EdgeId e1, EdgeId e2) const;

    void ClassifyPosition(size_t prev_pos, size_t cur_pos, const BidirectionalPath & path, PathScore &res) const;

    PathScore InternalCountMisassemblies(const BidirectionalPath &path) const;

//constructs longest sequence of consequetive ranges, stores result in used_mappings
    vector<MappingRange> FindBestRangeSequence(const set<MappingRange>& mappings) const;

    string ChromosomeByUniqueEdge(const EdgeId &e,
                                  const EdgesPositionHandler<Graph> &tmp_edge_pos,
                                  size_t &total) const;

    pair<MappingRange, size_t> Merge(const vector<MappingRange>& mappings) const;

    void FillPos(EdgeId e, const EdgesPositionHandler<Graph> &tmp_edge_pos);

    void ReportPathEndByPairedLib(const shared_ptr<path_extend::PairedInfoLibrary> paired_lib, EdgeId current_edge) const;

    void ReportPathEndByLongLib(const path_extend::BidirectionalPathSet &covering_paths, EdgeId current_edge) const;

    void ReportEdge(EdgeId e, double w) const;

    void ReportVariants(std::vector<pair<double, EdgeId>> &sorted_w) const;

    size_t GetSupportingPathCount(EdgeId e1, EdgeId e2, size_t lib_index) const;

    void TheoreticLenStats(vector<size_t> theoretic_lens) const;

    map<string, size_t> TotalAlignedLengths(const EdgesPositionHandler<Graph> &tmp_edge_pos, EdgeId e) const;;

    MappingPathT ConstructEdgeOrder(const std::string &chr_name) const;

    //returns lengths of mapped regions, divided by "unresolvable_len_"
    vector<size_t> MappedRegions(const MappingPathT &mapping_path) const;

    bool IsCloseToEnd(MappingRange range, const ChromosomeInfo &chr_info) const {
        auto last_range = gp_.edge_pos.GetUniqueEdgePosition(chr_info.EdgeAt(chr_info.size() - 1), chr_info.name());
        return range.initial_range.end_pos + SIGNIFICANT_LENGTH_LOWER_LIMIT > last_range.initial_range.end_pos;
    }

    bool IsCloseToStart(MappingRange range, const ChromosomeInfo &) const {
        return range.initial_range.start_pos <= SIGNIFICANT_LENGTH_LOWER_LIMIT;
    }

    DECL_LOGGER("GenomeConsistenceChecker");
public:
    GenomeConsistenceChecker(const conj_graph_pack &gp,
                             size_t max_gap,
                             double relative_max_gap /*= 0.2*/,
                             size_t unresolvable_len,
                             const ScaffoldingUniqueEdgeStorage &storage,
                             const vector<path_extend::GraphCoverageMap> &long_reads_cov_map,
                             const io::DataSet<config::LibraryData> reads) :
            gp_(gp),
            absolute_max_gap_(max_gap),
            relative_max_gap_(relative_max_gap),
            unresolvable_len_(unresolvable_len),
            storage_(storage),
            long_reads_cov_map_(long_reads_cov_map),
            reads_(reads) {
        //Fixme call outside
        Fill();
    }

    void Fill();

    PathScore CountMisassemblies(const BidirectionalPath &path) const;

    void CheckPathEnd(const BidirectionalPath &path) const;

    map<EdgeId, string> EdgeLabels() const;

};

}
