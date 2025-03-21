//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "test_utils.hpp"

#include "alignment/sequence_mapper.hpp"
#include "alignment/sequence_mapper_notifier.hpp"
#include "io/reads/converting_reader_wrapper.hpp"
#include "io/reads/rc_reader_wrapper.hpp"
#include "io/reads/read_stream_vector.hpp"
#include "io/reads/vector_reader.hpp"
#include "modules/graph_construction.hpp"
#include "paired_info/pair_info_filler.hpp"
#include "paired_info/paired_info_helpers.hpp"
#include "paired_info/weights.hpp"
#include "pipeline/graph_pack.hpp"
#include "pipeline/graph_pack_helpers.h"
#include "pipeline/sequence_mapper_gp_api.hpp"
#include "utils/stl_utils.hpp"

#include <unordered_set>

#include <gtest/gtest.h>

using namespace debruijn_graph;


//void ConstructGraphFromGenome(size_t k, Graph& g, EdgeIndex<Graph>& index/*, CoverageHandler<DeBruijnGraph>& coverage_handler*/
//, PairedInfoIndexT<Graph>& paired_index, const string& genome
//        , size_t read_size) {
//    typedef read_generator::ReadGenerator<read_generator::SmoothPositionChooser> Stream;
//    size_t coverage = 2 * read_size;
//    size_t gap = 0;
//    Stream raw_stream(2, read_size, genome, coverage, gap);
//    typedef io::RCReaderWrapper<io::PairedRead> RCStream;
//    RCStream read_stream(raw_stream);
//    ConstructGraphWithPairedInfo (k, g, index/*, coverage_handler*/,
//            paired_index, read_stream);
//}

using io::SingleRead;
using io::PairedRead;

typedef std::string MyRead;
typedef std::pair<MyRead, MyRead> MyPairedRead;
typedef std::string MyEdge;
typedef std::pair<MyEdge, MyEdge> MyEdgePair;
typedef std::multimap<MyEdgePair, std::pair<int, double>> EdgePairInfo;
typedef std::map<MyEdge, double> CoverageInfo;
typedef std::unordered_set<MyEdge> Edges;

namespace test_utils {

const Edges AddComplement(const Edges& edges) {
    Edges ans;
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        ans.insert(*it);
        ans.insert(ReverseComplement(*it));
    }
    return ans;
}

const CoverageInfo AddComplement(const CoverageInfo& coverage_info) {
    CoverageInfo ans;
    for (auto it = coverage_info.begin(); it != coverage_info.end(); ++it) {
        ans.insert(*it);
        ans.emplace(ReverseComplement((*it).first), (*it).second);
    }
    return ans;
}

const EdgePairInfo AddBackward(const EdgePairInfo& pair_info) {
    EdgePairInfo ans;
    for (auto it = pair_info.begin(); it != pair_info.end(); ++it) {
        ans.insert(*it);
        ans.insert({{(*it).first.second, (*it).first.first}, {-(*it).second.first, (*it).second.second}});
    }
    return ans;
}

const EdgePairInfo AddComplement(const EdgePairInfo& pair_info) {
    EdgePairInfo ans;
    for (auto it = pair_info.begin(); it != pair_info.end(); ++it) {
        ans.insert(*it);
        ans.insert({{ReverseComplement((*it).first.second), ReverseComplement((*it).first.first)}, (*it).second});
    }
    return ans;
}

void EdgesEqual(const Edges& s1, const Edges& s2) {
    ASSERT_EQ(s2.size(), s1.size());
    for (auto it = s1.begin(); it != s1.end(); ++it) {
        ASSERT_TRUE(s2.count(*it) > 0);
    }
}

const std::vector<io::SingleRead> MakeReads(const std::vector<MyRead>& reads) {
    std::vector<io::SingleRead> ans;
    for (size_t i = 0; i < reads.size(); ++i) {
        ans.push_back(reads[i]);
    }
    return ans;
}

const std::vector<PairedRead> MakePairedReads(const std::vector<MyPairedRead>& paired_reads, size_t insert_size) {
    DEBUG("Making paired reads");
    std::vector<PairedRead> ans;
    for (size_t i = 0; i < paired_reads.size(); ++i) {
        ans.push_back(PairedRead(paired_reads[i].first, paired_reads[i].second, insert_size));
    }
    DEBUG("Made paired reads");
    return ans;
}

void AssertEdges(Graph& g, const Edges& etalon_edges) {
    DEBUG("Asserting edges");
    Edges edges;
    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        edges.insert(g.EdgeNucls(*it).str());
    }
    EdgesEqual(edges, etalon_edges);
}

void AssertGraph(size_t k, const std::vector<std::string>& reads, const std::vector<std::string>& etalon_edges) {
    DEBUG("Asserting graph");
    typedef io::VectorReadStream<io::SingleRead> RawStream;
    Graph g(k);
    auto workdir = fs::tmp::make_temp_dir("tmp", "tests");

    io::ReadStreamList<io::SingleRead> streams(io::RCWrap<io::SingleRead>(RawStream(MakeReads(reads))));
    ConstructGraph(config::debruijn_config::construction(), workdir, streams, g);

    AssertEdges(g, AddComplement(Edges(etalon_edges.begin(), etalon_edges.end())));
}

bool EqualDouble(double d1, double d2) {
    return std::abs(d1 - d2) < 1e-5;
}

void AssertCoverage(Graph& g, const CoverageInfo& etalon_coverage) {
    DEBUG("Asserting coverage");
    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        auto etalon_cov_it = etalon_coverage.find(g.EdgeNucls(*it).str());
        ASSERT_NE(etalon_cov_it, etalon_coverage.end()) << "Etalon didn't contain edge '" << g.EdgeNucls(*it) << "'";
        ASSERT_TRUE(EqualDouble(g.coverage(*it), (*etalon_cov_it).second)) << "Coverage for edge '" << g.EdgeNucls(*it) << "' was " << g.coverage(*it) << " but etalon is " << (*etalon_cov_it).second;
    }
}

typedef omnigraph::de::PairedInfoIndexT<Graph> PairedIndex;
typedef omnigraph::de::PairInfo<EdgeId> PairInfo;
typedef std::vector<PairInfo> PairInfos;

template<class PairedIndex>
void AssertPairInfo(const Graph& g, /*todo const */PairedIndex& paired_index, const EdgePairInfo& etalon_pair_info) {
    for (auto i = omnigraph::de::pair_begin(paired_index); i != omnigraph::de::pair_end(paired_index); ++i) {
      for (auto j : *i) {
        PairInfo pair_info(i.first(), i.second(), j);
        if (pair_info.first == pair_info.second && rounded_d(pair_info) == 0)
          continue;

        auto my_edge_pair = std::make_pair(g.EdgeNucls(pair_info.first).str(), g.EdgeNucls(pair_info.second).str());
        auto equal_range = etalon_pair_info.equal_range(my_edge_pair);

        std::string my_edge_pair_str = "[" + my_edge_pair.first + ", " + my_edge_pair.second + "]";
        ASSERT_NE(equal_range.first, equal_range.second) << "Pair of edges " << my_edge_pair_str << " wasn't found in etalon";

        double etalon_weight = -1.0;

        for (auto range_it = equal_range.first; range_it != equal_range.second; ++range_it) {
          if ((*range_it).second.first == rounded_d(pair_info)) {
            etalon_weight = (*range_it).second.second;
          }
        }
        ASSERT_GT(etalon_weight, 0) << "Etalon didn't contain distance=" << rounded_d(pair_info) << " for edge pair " << my_edge_pair_str;
        ASSERT_TRUE(EqualDouble(etalon_weight, pair_info.weight())) << "Actual weight for edge pair " << my_edge_pair_str << " on distance " << rounded_d(pair_info) << " was " << pair_info.weight() << " but etalon is " <<  etalon_weight;
      }
    }
}

void AssertGraph(size_t k, const std::vector<MyPairedRead> &paired_reads, size_t /*rl*/, size_t insert_size,
                 const std::vector<MyEdge> &etalon_edges, const CoverageInfo &etalon_coverage,
                 const EdgePairInfo &etalon_pair_info) {
    typedef io::VectorReadStream<io::PairedRead> RawStream;

    DEBUG("Asserting graph with etalon data");

    DEBUG("Streams initialized");

    graph_pack::GraphPack gp(k, "tmp", 1);
    auto workdir = fs::tmp::make_temp_dir(gp.workdir(), "tests");
    DEBUG("Graph pack created");

    auto &graph = gp.get_mutable<Graph>();
    auto &index = gp.get_mutable<EdgeIndex<Graph>>();
    auto &flanking_cov = gp.get_mutable<omnigraph::FlankingCoverage<Graph>>();
    auto &kmer_mapper = gp.get_mutable<KmerMapper<Graph>>();
    auto &paired_indices = gp.get_mutable<omnigraph::de::UnclusteredPairedInfoIndicesT<Graph>>();

    RawStream paired_stream = MakePairedReads(paired_reads, insert_size);
    using SquashingWrapper = io::SquashingWrapper<io::PairedRead>;
    io::ReadStreamList<io::SingleRead> single_stream_vector(SquashingWrapper(std::move(paired_stream)));

    ConstructGraphWithCoverage(config::debruijn_config::construction(), workdir,
                               single_stream_vector, graph, index, flanking_cov);

    AssertEdges(graph, AddComplement(Edges(etalon_edges.begin(), etalon_edges.end())));
    AssertCoverage(graph, AddComplement(etalon_coverage));
    
    InitRRIndices(gp);
    kmer_mapper.Attach();
    EnsureBasicMapping(gp);

    io::ReadStreamList<io::PairedRead> paired_streams(std::move(single_stream_vector[0].recover<SquashingWrapper>().recover<RawStream>()));

    SequenceMapperNotifier notifier;
    LatePairedIndexFiller pif(graph, PairedReadCountWeight, 0, paired_indices[0]);
    notifier.Subscribe(&pif);
    notifier.ProcessLibrary(paired_streams, *MapperInstance(gp));
    
    AssertPairInfo(graph, paired_indices[0], AddComplement(AddBackward(etalon_pair_info)));
}

}
