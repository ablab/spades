//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graphio.hpp"
#include "test_utils.hpp"
#include "random_graph.hpp"
#include "assembly_graph/handlers/id_track_handler.hpp"
#include "io/binary/graph.hpp"
#include "io/binary/kmer_mapper.hpp"
#include "io/binary/paired_index.hpp"
#include "io/graph/gfa_reader.hpp"
#include "io/graph/gfa_writer.hpp"

#include <filesystem>
#include <gtest/gtest.h>

using namespace debruijn_graph;

template<typename T>
void CompareContainers(const T &lhs, const T &rhs) {
    auto i = lhs.begin();
    auto j = rhs.begin();
    for (; i != lhs.end() && j != rhs.end(); ++i, ++j) {
        EXPECT_EQ(*i, *j);
    }
    EXPECT_EQ(i, lhs.end());
    EXPECT_EQ(j, lhs.end());
}

template<typename I>
void CompareGraphIterators(I it1, I it2) {
    while(!it1.IsEnd() && !it2.IsEnd()) {
        EXPECT_EQ((*it1).int_id(), (*it2).int_id());
        ++it1;
        ++it2;
    }
    EXPECT_TRUE(it1.IsEnd());
    EXPECT_TRUE(it2.IsEnd());
}

const Graph &CommonGraph() {
    static Graph graph(55);
    if (!graph.size())
        RandomGraph<Graph>(graph, /*max_size*/100).Generate(/*iterations*/1000);
    //Add more self-conjugates
    auto it = graph.SmartVertexBegin(true);
    for (size_t i = 0; i < 5; ++i, ++it)
        graph.AddEdge(*it, graph.conjugate(*it), RandomSequence(60));
    return graph;
}

using namespace io::binary;
const char *file_name = "src/test/debruijn/graph_fragments/saves/test_save";

void CheckGFAInOut(const std::string &path_to_save,
                   const std::string &graph_name,
                   const std::filesystem::path &gfa_basename) {
    size_t K = 55;
    Graph graph(K);
    graphio::ScanBasicGraph(path_to_save, graph);
    std::filesystem::path gfa_path = gfa_basename / std::filesystem::path(graph_name + ".gfa");
    std::ofstream gfa_out(gfa_path);
    gfa::GFAWriter gfa_writer(graph, gfa_out);
    gfa_writer.WriteSegmentsAndLinks();
    gfa_out.close();

    std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());
    gfa::GFAReader gfa_reader(gfa_path);
    Graph gfa_graph(K);
    size_t local_k = gfa_reader.to_graph(gfa_graph, id_mapper.get());

    auto other_gfa_path = gfa_basename / std::filesystem::path(graph_name + "_other.gfa");
    std::ofstream other_out(other_gfa_path);
    gfa::GFAWriter other_writer(gfa_graph, other_out, io::MapNamingF<debruijn_graph::ConjugateDeBruijnGraph>(*id_mapper));
    other_writer.WriteSegmentsAndLinks();

    EXPECT_EQ(K, local_k);
    EXPECT_EQ(graph.size(), gfa_graph.size());
    EXPECT_EQ(graph.e_size(), gfa_graph.e_size());
}

TEST(Io, Order) {
    const auto &graph = CommonGraph();

    Save(file_name, graph);

    Graph new_graph(graph.k());
    Load(file_name, new_graph);

    CompareGraphIterators(graph.SmartVertexBegin(), new_graph.SmartVertexBegin());
    CompareGraphIterators(graph.SmartEdgeBegin(), new_graph.SmartEdgeBegin());
}

TEST(Io, PairedInfo) {
    using namespace omnigraph::de;
    using Index = UnclusteredPairedInfoIndexT<Graph>;
    const auto &graph = CommonGraph();

    Index pi(graph);
    RandomPairedIndex<Index>(pi, 100).Generate(100);
    //Add more self-conjugates
    RawPoint p(0, 42);
    auto it = graph.ConstEdgeBegin(true);
    for (size_t i = 0; i < 5; ++i, ++it)
        pi.Add(*it, graph.conjugate(*it), p);

    Save(file_name, pi);

    Index ni(graph);
    Load(file_name, ni);

    EXPECT_EQ(pi.size(), ni.size());
    for (auto pit = omnigraph::de::pair_begin(pi), nit = omnigraph::de::pair_begin(ni);
         pit != omnigraph::de::pair_end(pi); ++pit, ++nit) {
        EXPECT_EQ(pit.first(), nit.first());
        EXPECT_EQ(pit.second(), nit.second());
        EXPECT_EQ(pit->size(), nit->size());

        for (auto ppit = pit->begin(), npit = nit->begin();
             ppit != pit->end(); ++ppit, ++npit) {
            EXPECT_EQ(ppit->weight, npit->weight);
            EXPECT_EQ(ppit->d, npit->d);
        }
    }
}

TEST(Io, KmerMapper) {
    const auto &graph = CommonGraph();

    KmerMapper<Graph> kmer_mapper(graph);
    RandomKmerMapper<Graph>(kmer_mapper).Generate(100);

    Save(file_name, kmer_mapper);

    KmerMapper<Graph> new_mapper(graph);
    Load(file_name, new_mapper);

    CompareContainers(kmer_mapper, new_mapper);
}

TEST(Io, GFADBG) {
    std::filesystem::path gfa_out_base("src/test/debruijn/graph_fragments/gfa_saves");

    CheckGFAInOut("src/test/debruijn/graph_fragments/tipobulge/tipobulge", "tipobulge", gfa_out_base);
    CheckGFAInOut("src/test/debruijn/graph_fragments/topology_ec/self_comp", "self_comp", gfa_out_base);
    //fixme support 0-in-2-out DBG vertices in GFAWriter
//    CheckGFAInOut("src/test/debruijn/graph_fragments/topology_ec/big_bad", "big_bad", gfa_out_base);
}
