#pragma once
#include "assembly_graph/index/edge_multi_index.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"

#include "modules/alignment/edge_index_refiller.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "modules/alignment/gap_info.hpp"

#include "modules/alignment/pacbio/pacbio_read_structures.hpp"
#include "modules/alignment/pacbio/gap_filler.hpp"
#include "modules/alignment/pacbio/pac_index.hpp"

namespace sensitive_aligner {

struct OneReadMapping {
    std::vector<vector<debruijn_graph::EdgeId>> main_storage;
    std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId>> bwa_paths;
    std::vector<GapDescription> gaps;
    std::vector<PathRange> read_ranges;
    OneReadMapping(const std::vector<vector<debruijn_graph::EdgeId>> &main_storage_,
                   const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId>> &bwa_paths_,
                   const std::vector<GapDescription>& gaps_,
                   const std::vector<PathRange> &read_ranges_) :
            main_storage(main_storage_), bwa_paths(bwa_paths_), gaps(gaps_), read_ranges(read_ranges_){}
};

typedef std::pair<QualityRange, int> ColoredRange;

class GAligner {
    PacBioMappingIndex pac_index_;
    const debruijn_graph::Graph &g_;
    debruijn_graph::config::pacbio_processor pb_config_;
//TODO:: shouldn't it be somewhere in debruijn_graph::config?
    GapClosingConfig gap_cfg_;
    GapFiller gap_filler_;

    void ProcessCluster(const Sequence &s,
                             std::vector<QualityRange> &cur_cluster,
                             std::vector<QualityRange> &start_clusters,
                             std::vector<QualityRange> &end_clusters,
                             std::vector<vector<debruijn_graph::EdgeId> > &sorted_edges,
                             std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
                             std::vector<bool> &block_gap_closer) const;
    
    void FillGapsInCluster(const vector<QualityRange> &cur_cluster,
                      const Sequence &s,
                      std::vector<vector<debruijn_graph::EdgeId> > &edges,
                      std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &bwa_hits) const;

    std::pair<int, int> GetPathLimits(const QualityRange &a,
                                      const QualityRange &b,
                                      int s_add_len, int e_add_len) const;

    bool TopologyGap(debruijn_graph::EdgeId first, debruijn_graph::EdgeId second, bool oriented) const;

    OneReadMapping AddGapDescriptions(const std::vector<QualityRange> &start_clusters,
                              const std::vector<QualityRange> &end_clusters,
                              const std::vector<vector<debruijn_graph::EdgeId> > &sorted_edges,
                              const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
                              const std::vector<PathRange> &read_ranges,
                              const Sequence &s,
                              const std::vector<bool> &block_gap_closer) const;
    void RestoreEnds(const Sequence &s,
                     const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
                     std::vector<vector<debruijn_graph::EdgeId> > &sorted_edges,
                     PathRange &cur_range) const;
public:
    OneReadMapping GetReadAlignment(const io::SingleRead &read) const;
    GAligner(const debruijn_graph::Graph &g,
             debruijn_graph::config::pacbio_processor pb_config,
             alignment::BWAIndex::AlignmentMode mode,
             GapClosingConfig gap_cfg = GapClosingConfig())
            : pac_index_(g, pb_config, mode), g_(g), pb_config_(pb_config), gap_cfg_(gap_cfg), gap_filler_(g, pb_config, gap_cfg){}



};
}