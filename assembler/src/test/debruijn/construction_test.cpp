//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "io/reads/vector_reader.hpp"
#include "io/reads/read_stream_vector.hpp"
#include "io/reads/rc_reader_wrapper.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/filesystem/temporary.hpp"
#include "pipeline/graph_pack.hpp" // FIXME: get rid of it
#include "modules/graph_construction.hpp"
#include "modules/alignment/edge_index.hpp"

#include "test_utils.hpp"
#include "tmp_folder_fixture.hpp"

#include <vector>
#include <set>
#include <string>

#include <gtest/gtest.h>

using namespace debruijn_graph;
using namespace test_utils;

class GraphConstruction : public ::testing::Test, public TmpFolderFixture {};

//todo rename tests
TEST_F( GraphConstruction, SimpleThread ) {
    std::vector<std::string> reads = { "ACAAACCACCA" };
//    vector<string> edges = { "ACAAACCACCA" };
    AssertGraph (5, reads, reads);
}

TEST_F( GraphConstruction, SimpleThread2 ) {
    std::vector<std::string> reads = { "ACAAACCACCC", "AAACCACCCAC" };
    std::vector<std::string> edges = { "ACAAACCACCCAC" };
    AssertGraph (5, reads, edges);
}

TEST_F( GraphConstruction, SplitThread ) {
    std::vector<std::string> reads = { "ACAAACCACCA", "ACAAACAACCC" };
    std::vector<std::string> edges = { "ACAAAC", "CAAACCACCA", "CAAACAACCC" };
    AssertGraph (5, reads, edges);
}

TEST_F( GraphConstruction, SplitThread2 ) {
    std::vector<std::string> reads = { "ACAAACCACCA", "ACAAACAACCA" };
    std::vector<std::string> edges = { "AACCACCA", "ACAAAC", "CAAACCA", "CAAACAACCA" };
    AssertGraph (5, reads, edges);
}

TEST_F( GraphConstruction, Buldge ) {
    std::vector<std::string> reads = { "ACAAAACACCA", "ACAAACCACCA" };
//    vector<string> edges = { "ACAAAACACCA", "ACAAACCACCA" };
    AssertGraph (5, reads, reads);
}

TEST_F( GraphConstruction, CondenseSimple ) {
    std::vector<std::string> reads = { "CGAAACCAC", "CGAAAACAC", "AACCACACC", "AAACACACC" };
    std::vector<std::string> edges = { "CGAAAACACAC", "CACACC", "CGAAACCACAC" };
    AssertGraph (5, reads, edges);
}

void CheckIndex(const std::vector<std::string> &reads, const std::string &tmpdir, size_t k) {
    typedef io::VectorReadStream<io::SingleRead> RawStream;
    GraphPack gp(k, tmpdir, 0);
    auto workdir = fs::tmp::make_temp_dir(gp.workdir(), "tests");
    io::ReadStreamList<io::SingleRead> streams(io::RCWrap<io::SingleRead>(RawStream(MakeReads(reads))));
    auto &graph = gp.get_mutable<Graph>();
    auto &index = gp.get_mutable<EdgeIndex<Graph>>();
    ConstructGraphWithIndex(config::debruijn_config::construction(), workdir, streams, graph, index);
    auto &stream = streams.back();
    stream.reset();
    io::SingleRead read;
    while (!stream.eof()) {
        stream >> read;
        RtSeq kmer = read.sequence().start<RtSeq>(k + 1) >> 'A';
        for (size_t i = k; i < read.size(); i++) {
            kmer = kmer << read[i];
            EXPECT_TRUE(index.contains(kmer));
        }
    }
}

TEST_F( GraphConstruction, TestKmerStoringIndex ) {
    std::vector<std::string> reads = { "CGAAACCAC", "CGAAAACAC", "AACCACACC", "AAACACACC" };
    CheckIndex(reads, tmp_folder(), 5);
}

TEST_F( GraphConstruction, TestKmerFreeIndex ) {
    std::vector<std::string> reads = { "CGAAACCAC", "CGAAAACAC", "AACCACACC", "AAACACACC" };
    CheckIndex(reads, tmp_folder(), 5);
}

TEST_F( GraphConstruction, SimpleTestEarlyPairedInfo ) {
    std::vector<MyPairedRead> paired_reads = {{"CCCAC", "CCACG"}, {"ACCAC", "CCACA"}};
    std::vector<MyEdge> edges = {"CCCA", "ACCA", "CCAC", "CACG", "CACA"};
    CoverageInfo coverage_info = {{"CCCA", 1}, {"ACCA", 1}, {"CCAC", 4}, {"CACG", 1}, {"CACA", 1}};
    EdgePairInfo edge_pair_info = {{{"CCCA", "CACG"}, {2, 1.0}}, {{"ACCA", "CACA"}, {2, 1.0}}
        , {{"CCCA", "CCAC"}, {1, 1.0}}, {{"ACCA", "CCAC"}, {1, 1.0}}
        , {{"CCAC", "CACG"}, {1, 1.0}}, {{"CCAC", "CACA"}, {1, 1.0}}};

    AssertGraph(3, paired_reads, 5, 6, edges, coverage_info, edge_pair_info);
}
