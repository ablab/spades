//***************************************************************************
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/assembly_graph/core/graph.hpp"

namespace gfa_tools {

struct StatInfo;

void CountLinksNumber(const debruijn_graph::Graph& g, StatInfo& info);

void CountDeadEndsInfo(const debruijn_graph::Graph& g, StatInfo& info);

std::map<double, size_t> GetPercentile(const std::vector<size_t>& v, const std::vector<double>& percentiles);

void CountSequencesLengthInfo(const debruijn_graph::Graph& g, const std::vector<double>& percentiles,
                                                                StatInfo& info);

void CountConnectedComponentsInfo(const debruijn_graph::Graph& g, StatInfo& info);

size_t GetNx(const debruijn_graph::Graph& g, double percent);

void CountPercentileCoverage(const debruijn_graph::Graph& g, const std::vector<double>& percentiles,
                                                                StatInfo& info);

void CountEstimatedGenomeLength(const debruijn_graph::Graph& g, StatInfo& info);

void CountSequencesDegreeInfo(const debruijn_graph::Graph& g, StatInfo& info);

void CalculateStat(const debruijn_graph::Graph& g, const std::vector<double>& n_percentiles,
               const std::vector<double>& length_percentiles,
               const std::vector<double>& sequences_coverage_percentiles,
               const std::filesystem::path& yaml_output_path);
}
