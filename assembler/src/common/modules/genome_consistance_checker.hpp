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

    vector<string> Chromosomes() const {
        vector<string> answer;
        push_back_all(answer, key_set(chr_infos_));
        return answer;
    }
};

class GenomeConsistenceChecker {
    typedef omnigraph::MappingPath<EdgeId> MappingPathT;
    const conj_graph_pack &gp_;
    //EdgesPositionHandler<Graph> &position_handler_;
    const ScaffoldingUniqueEdgeStorage &storage_;
    const size_t absolute_max_gap_;
    const double relative_max_gap_;
//Edges containing zero point for each reference
//TODO: do we need circular/linear chromosomes support?
    set<EdgeId> circular_edges_;
    const size_t unresolvable_len_;
    const vector<shared_ptr<path_extend::GraphCoverageMap>> &long_reads_cov_map_;
    GenomeInfo genome_info_;

    bool Consequent(const MappingRange &mr1, const MappingRange &mr2) const ;
    void PrintMisassemblyInfo(EdgeId e1, EdgeId e2) const;

    PathScore InternalCountMisassemblies(const BidirectionalPath &path) const;
//map from unique edges to their order in genome spelling;
// for each chromosome provides its spelling

//constructs longest sequence of consequetive ranges, stores result in used_mappings
    vector<MappingRange> FindBestRangeSequence(const set<MappingRange>& mappings) const;

    string ChromosomeByUniqueEdge(const EdgeId &e,
                                  const EdgesPositionHandler<Graph> &tmp_edge_pos,
                                  size_t &total) const {
        DEBUG("Positioning edge " << gp_.g.str(e));
        map<string, size_t> total_al_lens = TotalAlignedLengths(tmp_edge_pos, e);
        total = 0;
        for (size_t c : value_set(total_al_lens))
            total += c;

        if (total > size_t(math::round((double) gp_.g.length(e) * 1.5))) {
            DEBUG("Not unique, excluding");
            return "";
        }

        string chr = "";
        size_t max_l = 0;
        for (const auto &p : total_al_lens) {
            if (p.second > max_l) {
                max_l = p.second;
                chr = p.first;
            }
        }

        DEBUG("Most likely chromosome " << chr << ". Mapped bp: " << max_l);
        //TODO: support non-unique edges;
        if (max_l < size_t(math::round((double) gp_.g.length(e) * 0.5))) {
            DEBUG("Too small a portion mapped. Edge not used");
            return "";
        }

        return chr;
    }

    pair<MappingRange, size_t> Merge(const vector<MappingRange>& mappings) const {
        VERIFY(mappings.size() > 0);

        MappingRange mr = mappings.front();
        size_t total_mapped = mr.initial_range.size();
        for (size_t i = 1; i < mappings.size(); ++i) {
            total_mapped += mappings[i].initial_range.size();
            //FIXME why do we need merge?
            mr = mr.Merge(mappings[i]);
        }
        return make_pair(mr, total_mapped);
    }

    void FillPos(EdgeId e, const EdgesPositionHandler<Graph> &tmp_edge_pos) {
        size_t total_mapped;
        string chr = ChromosomeByUniqueEdge(e, tmp_edge_pos, total_mapped);
        if (chr.empty())
            return;

        auto mapping_info = Merge(FindBestRangeSequence(tmp_edge_pos.GetEdgePositions(e, chr)));

        //FIXME what is the logic here?
        //used less that 0.9 of aligned length
        VERIFY(total_mapped > mapping_info.second);
        if ((total_mapped - mapping_info.second) * 10 >=  gp_.g.length(e)) {
            INFO ("Edge " << gp_.g.int_id(e) << " length "<< gp_.g.length(e)  << "is potentially misassembled! mappings: ");
            for (auto mp : tmp_edge_pos.GetEdgePositions(e, chr)) {
                INFO("mp_range "<< mp.mapped_range.start_pos << " - " << mp.mapped_range.end_pos << " init_range " << mp.initial_range.start_pos << " - " << mp.initial_range.end_pos );
                if (mp.initial_range.start_pos < absolute_max_gap_) {
                    INFO ("Fake(linear order) misassembly on edge "<< e.int_id());
                    circular_edges_.insert(e);
                }
            }
        }

        gp_.edge_pos.AddEdgePosition(e, chr, mapping_info.first);
    }

    void ReportEdge(EdgeId e, double w) const;
    void ReportVariants(std::vector<pair<double, EdgeId>> &sorted_w) const;
    size_t GetSupportingPathCount(EdgeId e1, EdgeId e2, size_t lib_index) const;

//spells genome in language of long unique edges from storage;
    void TheoreticLenStats(vector<size_t> theoretic_lens) const {
        size_t total_len = std::accumulate(theoretic_lens.begin(), theoretic_lens.end(),
                                           0, std::plus<size_t>());

        std::sort(theoretic_lens.begin(), theoretic_lens.end());
        std::reverse(theoretic_lens.begin(), theoretic_lens.end());
        size_t cur = 0;
        size_t i = 0;
        while (cur < total_len / 2) {
            cur += theoretic_lens[i];
            i++;
        }
        INFO("Assuming gaps of length > " << storage_.GetMinLength() << " unresolvable..");
        if (theoretic_lens.size() > 0)
            INFO("Rough estimates on N50/L50:" << theoretic_lens[i - 1] << " / " << i - 1 << " with len " << total_len);
    }
    DECL_LOGGER("GenomeConsistenceChecker");

    map<string, size_t> TotalAlignedLengths(const EdgesPositionHandler<Graph> &tmp_edge_pos, EdgeId e) const {
        map<string, size_t> chr2len;
        for (const auto &edge_pos: tmp_edge_pos.GetEdgePositions(e)) {
            chr2len[edge_pos.contigId] += edge_pos.mr.initial_range.size();
        }
        return chr2len;
    };

public:
    GenomeConsistenceChecker(const conj_graph_pack &gp, const ScaffoldingUniqueEdgeStorage &storage, size_t max_gap,
                             double relative_max_gap /*= 0.2*/, size_t unresolvable_len,
                             const vector<shared_ptr<path_extend::GraphCoverageMap>> &long_reads_cov_map) :
            gp_(gp),
            storage_(storage), absolute_max_gap_(max_gap),
            relative_max_gap_(relative_max_gap),
            unresolvable_len_(unresolvable_len), long_reads_cov_map_(long_reads_cov_map) {
        //Fixme call outside
        Fill();
    }

    //returns lengths of mapped regions, divided by "unresolvable_len_"
    vector<size_t> MappedRegions(const MappingPathT &mapping_path) const {
        vector<size_t> mapped_regions;
        size_t pos = mapping_path.front().second.initial_range.start_pos;
        for (size_t i = 0; i < mapping_path.size(); i++) {
            auto current_range = mapping_path[i].second;
            INFO("Pos: " << i << " init_range " << current_range.initial_range
                         << " mapped to edge " << gp_.g.str(mapping_path[i].first)
                         << " range " << current_range.mapped_range);

            size_t curr_start = current_range.initial_range.start_pos;
            if (i > 0) {
                auto prev_range = mapping_path[i - 1].second;
                size_t prev_end = prev_range.initial_range.end_pos;
                if (curr_start - prev_end > unresolvable_len_) {
                    INFO ("Large gap " << current_range.initial_range.start_pos -
                                          prev_range.initial_range.end_pos);
                    mapped_regions.push_back(prev_end - pos);
                    pos = curr_start;
                }
            }
        }
        mapped_regions.push_back(mapping_path.back().second.initial_range.end_pos - pos);
        return mapped_regions;
    }

    void Fill() {
        gp_.edge_pos.clear();
        if (!gp_.edge_pos.IsAttached()) {
            gp_.edge_pos.Attach();
        }

        //FIXME set the parameters to something more reasonable
        EdgesPositionHandler<Graph> tmp_edge_pos(gp_.g, 0, 0);
        visualization::position_filler::PosFiller<Graph> pos_filler(gp_.g, MapperInstance(gp_), tmp_edge_pos);

        for (const auto &chr: gp_.genome.GetChromosomes()) {
            pos_filler.Process(chr.sequence, "0_" + chr.name);
            pos_filler.Process(ReverseComplement(chr.sequence), "1_" + chr.name);
        }

        for (auto it = gp_.g.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            FillPos(*it, tmp_edge_pos);
        }

        vector<size_t> theoretic_lens;
        for (const auto &prefix: vector<std::string>{"0_", "1_"}) {
            for (const auto &chr: gp_.genome.GetChromosomes()) {
                string label = prefix + chr.name;
                INFO("Spelling label " << label);
                auto mapping_path = ConstructEdgeOrder(label);
                genome_info_.AddInfo(ChromosomeInfo(label, mapping_path));

                push_back_all(theoretic_lens, MappedRegions(mapping_path));
            }
        }

        TheoreticLenStats(theoretic_lens);
    }

    PathScore CountMisassemblies(const BidirectionalPath &path) const;
    void CheckPathEnd(const BidirectionalPath &path) const;
    MappingPathT ConstructEdgeOrder(const std::string &chr_name) const;

    map<EdgeId, string> EdgeLabels() const {
        INFO("Constructing reference labels");
        map<EdgeId, string> answer;
        size_t count = 0;
        for (const auto &chr: genome_info_.Chromosomes()) {
            const auto &chr_info = genome_info_.ChrInfo(chr);
            for (size_t pos = 0; pos < chr_info.size(); ++pos) {
                EdgeId e = chr_info.EdgeAt(pos);
                auto mr = gp_.edge_pos.GetUniqueEdgePosition(e, chr);
                VERIFY(!answer.count(e));
                answer[e] += chr +
                             "order: " + ToString(count) +
                             "\n mapped range: " +
                             ToString(mr.mapped_range.start_pos) + " : "
                             + ToString(mr.mapped_range.end_pos) +
                             "\n init range: " +
                             ToString(mr.initial_range.start_pos) + " : "
                             + ToString(mr.initial_range.end_pos) + "\n";
            }
        }
        return answer;
    }

};

}
