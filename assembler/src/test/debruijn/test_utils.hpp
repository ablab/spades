//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/io/reads/single_read.hpp"
#include "common/io/reads/paired_read.hpp"

#include <unordered_set>
#include <string>
#include <map>

typedef std::string MyRead;
typedef std::pair<MyRead, MyRead> MyPairedRead;
typedef std::string MyEdge;
typedef std::pair<MyEdge, MyEdge> MyEdgePair;
typedef std::multimap<MyEdgePair, std::pair<int, double>> EdgePairInfo;
typedef std::map<MyEdge, double> CoverageInfo;
typedef std::unordered_set<MyEdge> Edges;

namespace debruijn_graph {
class DeBruijnGraph;
typedef DeBruijnGraph ConjugateDeBruijnGraph;
typedef ConjugateDeBruijnGraph Graph;
}

namespace test_utils {
const Edges AddComplement(const Edges& edges);
const CoverageInfo AddComplement(const CoverageInfo& coverage_info);
const EdgePairInfo AddBackward(const EdgePairInfo& pair_info);
const EdgePairInfo AddComplement(const EdgePairInfo& pair_info);
void EdgesEqual(const Edges& s1, const Edges& s2);
const io::SingleRead MakeRead(const MyRead& read);
const std::vector<io::SingleRead> MakeReads(const std::vector<MyRead>& reads);
const std::vector<io::PairedRead> MakePairedReads(const std::vector<MyPairedRead>& paired_reads, size_t insert_size);
void AssertEdges(debruijn_graph::Graph& g, const Edges& etalon_edges);
void AssertGraph(size_t k, const std::vector<std::string>& reads, const std::vector<std::string>& etalon_edges);
bool EqualDouble(double d1, double d2);
void AssertCoverage(debruijn_graph::Graph& g, const CoverageInfo& etalon_coverage);
void AssertGraph(size_t k, const std::vector<MyPairedRead> &paired_reads, size_t /*rl*/, size_t insert_size,
                 const std::vector<MyEdge> &etalon_edges, const CoverageInfo &etalon_coverage,
                 const EdgePairInfo &etalon_pair_info);
};
