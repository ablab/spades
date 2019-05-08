//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "assembly_graph/index/edge_multi_index.hpp"

#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "modules/alignment/gap_info.hpp"

#include "modules/alignment/pacbio/pacbio_read_structures.hpp"
#include "modules/alignment/pacbio/gap_filler.hpp"
#include "modules/alignment/pacbio/pac_index.hpp"

namespace sensitive_aligner {

struct OneReadMapping {
  std::vector<std::vector<debruijn_graph::EdgeId>> edge_paths;
  std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId>> bwa_paths;
  std::vector<GapDescription> gaps;
  std::vector<PathRange> read_ranges;
  std::vector<int> scores;
  std::vector<std::string> aln_type;

  OneReadMapping(const std::vector<std::vector<debruijn_graph::EdgeId>> &edge_paths_,
                 const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId>> &bwa_paths_,
                 const std::vector<GapDescription>& gaps_,
                 const std::vector<PathRange> &read_ranges_,
                 const std::vector<int> &scores_ = std::vector<int>(),
                 const std::vector<std::string> &aln_type_ = std::vector<std::string>()) :
    edge_paths(edge_paths_), bwa_paths(bwa_paths_), gaps(gaps_), read_ranges(read_ranges_), scores(scores_), aln_type(aln_type_) {}
};

typedef std::pair<QualityRange, int> ColoredRange;

class GAligner {
 public:
  OneReadMapping GetAlignment(const io::SingleRead &read) const;
  OneReadMapping GetProteinAlignment(const io::SingleRead &read) const;
  OneReadMapping GetReadAlignment(const io::SingleRead &read) const;
  
  GAligner(const debruijn_graph::Graph &g,
           const GAlignerConfig &cfg)
    : pac_index_(g, cfg.pb, cfg.data_type), g_(g), cfg_(cfg), mode_(cfg.data_type), restore_ends_(cfg.restore_ends), gap_filler_(g, cfg) {
      if (alignment::BWAIndex::AlignmentMode::Protein == cfg.data_type) {
        gap_filler_.ProteinModeOn();
      }
    }

  GAligner(const debruijn_graph::Graph &g,
           const debruijn_graph::config::pacbio_processor &pb_config,
           const alignment::BWAIndex::AlignmentMode &mode)
    : pac_index_(g, pb_config, mode), g_(g), cfg_(GAlignerConfig(pb_config, mode)), mode_(mode), restore_ends_(false), gap_filler_(g, GAlignerConfig(pb_config, mode)) {
      if (alignment::BWAIndex::AlignmentMode::Protein == mode) {
        gap_filler_.ProteinModeOn();
      }
    }


 private:
  PacBioMappingIndex pac_index_;
  const debruijn_graph::Graph &g_;
  const GAlignerConfig cfg_;
  const alignment::BWAIndex::AlignmentMode mode_;
  bool restore_ends_;
  GapFiller gap_filler_;

  void ProcessCluster(const Sequence &s,
                      std::vector<QualityRange> &cur_cluster,
                      std::vector<QualityRange> &start_clusters,
                      std::vector<QualityRange> &end_clusters,
                      std::vector<std::vector<debruijn_graph::EdgeId> > &sorted_edges,
                      std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
                      std::vector<bool> &block_gap_closer) const;

  void FillGapsInCluster(const std::vector<QualityRange> &cur_cluster,
                         const Sequence &s,
                         std::vector<std::vector<debruijn_graph::EdgeId> > &edges,
                         std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &bwa_hits) const;

  std::pair<int, int> GetPathLimits(const QualityRange &a,
                                    const QualityRange &b,
                                    int s_add_len, int e_add_len) const;

  bool TopologyGap(debruijn_graph::EdgeId first, debruijn_graph::EdgeId second, bool oriented) const;

  OneReadMapping AddGapDescriptions(const std::vector<QualityRange> &start_clusters,
                                    const std::vector<QualityRange> &end_clusters,
                                    const std::vector<std::vector<debruijn_graph::EdgeId> > &sorted_edges,
                                    const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &sorted_bwa_hits,
                                    const std::vector<PathRange> &read_ranges,
                                    const Sequence &s,
                                    const std::vector<bool> &block_gap_closer) const;
  int RestoreEndsF(const Sequence &s,
                   int end,
                   std::vector<debruijn_graph::EdgeId> &sorted_edges,
                   PathRange &cur_range) const;

  int RestoreEndsB(const Sequence &s,
                   int start,
                   std::vector<debruijn_graph::EdgeId> &sorted_edges,
                   PathRange &cur_range) const;

  std::pair<int, std::string> FindScore(const PathRange &range, const std::vector<debruijn_graph::EdgeId> &edges, const Sequence &s) const;
};
}