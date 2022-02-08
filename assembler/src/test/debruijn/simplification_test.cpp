//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//#include "modules/simplification/parallel_simplification_algorithms.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "stages/simplification_pipeline/single_cell_simplification.hpp"
#include "stages/simplification_pipeline/rna_simplification.hpp"
#include "pipeline/graph_pack.hpp"

#include "graphio.hpp"
#include "tmp_folder_fixture.hpp"

#include <gtest/gtest.h>

using namespace debruijn_graph;
using namespace debruijn_graph::config;

debruijn_config::simplification::bulge_remover standard_br_config_generation() {
    debruijn_config::simplification::bulge_remover br_config;
    br_config.enabled = true;
    br_config.main_iteration_only = false;
    br_config.max_bulge_length_coefficient = 4;
    br_config.max_additive_length_coefficient = 0;
    br_config.max_coverage = 1000.;
    br_config.max_relative_coverage = 1.2;
    br_config.max_delta = 3;
    br_config.max_number_edges = std::numeric_limits<size_t>::max();
    br_config.dijkstra_vertex_limit = std::numeric_limits<size_t>::max();
    br_config.max_relative_delta = 0.1;
    //fixme test both
    br_config.parallel = false;//true;
    br_config.buff_size = 10000;
    br_config.buff_cov_diff = 2.;
    br_config.buff_cov_rel_diff = 0.2;
    return br_config;
}

//static size_t standard_read_length() {
//    return 100;
//}

debruijn_config::simplification::bulge_remover standard_br_config() {
    static debruijn_config::simplification::bulge_remover br_config = standard_br_config_generation();
    return br_config;
}

debruijn_config::simplification::erroneous_connections_remover standard_ec_config_generation() {
    debruijn_config::simplification::erroneous_connections_remover ec_config;
    ec_config.condition = "{ cb 30 , ec_lb 20 }";
    return ec_config;
}

debruijn_config::simplification::erroneous_connections_remover standard_ec_config() {
    static debruijn_config::simplification::erroneous_connections_remover ec_config = standard_ec_config_generation();
    return ec_config;
}

debruijn_config::simplification::topology_based_ec_remover topology_based_ec_config_generation() {
    debruijn_config::simplification::topology_based_ec_remover tec_config;
    tec_config.max_ec_length_coefficient = 20;
    tec_config.plausibility_length = 200;
    tec_config.uniqueness_length = 1500;
    return tec_config;
}

debruijn_config::simplification::max_flow_ec_remover max_flow_based_ec_config_generation() {
    debruijn_config::simplification::max_flow_ec_remover mfec_config;
    mfec_config.enabled = true;
    mfec_config.max_ec_length_coefficient = 20;
    mfec_config.plausibility_length = 200;
    mfec_config.uniqueness_length = 3000;
    return mfec_config;
}

debruijn_config::simplification::topology_based_ec_remover standard_tec_config() {
    static debruijn_config::simplification::topology_based_ec_remover tec_config = topology_based_ec_config_generation();
    return tec_config;
}

debruijn_config::simplification::max_flow_ec_remover standard_mfec_config() {
    static debruijn_config::simplification::max_flow_ec_remover tec_config = max_flow_based_ec_config_generation();
    return tec_config;
}

debruijn_config::simplification::tip_clipper standard_tc_config_generation() {
    debruijn_config::simplification::tip_clipper tc_config;
    tc_config.condition = "{ tc_lb 2.5 , cb 1000. , rctc 1.2 }";
    return tc_config;
}

debruijn_config::simplification::tip_clipper standard_tc_config() {
    static debruijn_config::simplification::tip_clipper tc_config = standard_tc_config_generation();
    return tc_config;
}

debruijn_config::simplification::relative_coverage_comp_remover standard_rcc_config() {
    debruijn_config::simplification::relative_coverage_comp_remover rcc;
    //rather unrealistic value =)
    rcc.enabled = true;
    rcc.coverage_gap = 2.;
    rcc.length_coeff = 2.;
    rcc.tip_allowing_length_coeff = 2.;
    rcc.max_ec_length_coefficient = 65;
    rcc.max_coverage_coeff = 10000.;
    rcc.vertex_count_limit = 10;
    return rcc;
}

debruijn::simplification::SimplifInfoContainer standard_simplif_relevant_info() {
    debruijn::simplification::SimplifInfoContainer info(config::pipeline_type::base);
    return info.set_read_length(100)
            .set_detected_coverage_bound(10.)
            .set_main_iteration(true)
            .set_chunk_cnt(1);
}

std::string graph_fragment_root() {
    return "./src/test/debruijn/graph_fragments/";
}

void PrintGraph(const Graph & g) {
    for (VertexId v: g) {
        for (EdgeId e: g.OutgoingEdges(v)) {
            std::cout << g.int_id(e) << ":" << g.int_id(g.EdgeStart(e)) << " " << g.int_id(g.EdgeEnd(e)) << std::endl;
        }
    }
    std::cout << std::endl;
}

int trivial_false(EdgeId, const std::vector<EdgeId>&){
    return 0;
};

void DefaultClipTips(Graph& graph) {
    debruijn::simplification::TipClipperInstance(graph, standard_tc_config(), standard_simplif_relevant_info())->Run();
}

void DefaultRemoveBulges(Graph& graph) {
    debruijn::simplification::BRInstance(graph, standard_br_config(), standard_simplif_relevant_info(), trivial_false)->Run();
}

class Simplification : public ::testing::Test, public TmpFolderFixture {
};

TEST_F( Simplification,  SimpleTipClipperTest ) {
    ConjugateDeBruijnGraph g(55);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/simpliest_tip/simpliest_tip", g));

    DefaultClipTips(g);

    EXPECT_EQ(4, g.size());
}

TEST_F( Simplification,  SimpleBulgeRemovalTest ) {
    Graph g(55);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/simpliest_bulge/simpliest_bulge", g));

    DefaultRemoveBulges(g);

    EXPECT_EQ(4, g.size());
}

TEST_F( Simplification,  TipobulgeTest ) {
    Graph g(55);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/tipobulge/tipobulge", g));

    DefaultClipTips(g);

    DefaultRemoveBulges(g);

    EXPECT_EQ(16, g.size());
}

TEST_F( Simplification,  SimpleECTest ) {
    Graph g(55);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g));

    debruijn_config::simplification::erroneous_connections_remover ec_config;
    ec_config.condition = "{ icb 7000 , ec_lb 20 }";

    debruijn::simplification::ECRemoverInstance(g, ec_config, standard_simplif_relevant_info())->Run();

    EXPECT_EQ(16, g.size());
}

TEST_F( Simplification,  SimpleIterECTest ) {
    Graph g(55);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g));

    debruijn_config::simplification::erroneous_connections_remover ec_config;
    ec_config.condition = "{ icb 7000 , ec_lb 20 }";

    auto ec_remover_ptr = debruijn::simplification::ECRemoverInstance(g, ec_config, standard_simplif_relevant_info());

    AlgorithmRunningHelper<Graph>::IterativeThresholdsRun(*ec_remover_ptr, 2);
    EXPECT_EQ(16, g.size());
}

TEST_F( Simplification,  IterECTest ) {
    Graph g(55);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g));

    debruijn_config::simplification::erroneous_connections_remover ec_config;
    ec_config.condition = "{ icb 7000 , ec_lb 20 }";

    auto ec_remover_ptr = debruijn::simplification::ECRemoverInstance(g, ec_config, standard_simplif_relevant_info());

    ec_remover_ptr->Run(false, 0.5);
    EXPECT_EQ(20, g.size());

    ec_remover_ptr->Run(false, 1.0);
    EXPECT_EQ(16, g.size());
}

TEST_F( Simplification,  IterUniquePath ) {
    Graph g(55);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g));

    auto tec_config = standard_tec_config();
    while(debruijn::simplification::TopologyRemoveErroneousEdges<Graph>(g, tec_config, 0)) {
    }

    EXPECT_EQ(16, g.size());
}

//TEST( Simplification,  MFIterUniquePath ) {
//    Graph g(55);
//    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g));
//
//    debruijn_config::simplification::max_flow_ec_remover mfec_config = standard_mfec_config();
//    mfec_config.uniqueness_length = 500;
//    MaxFlowRemoveErroneousEdges<Graph>(g, mfec_config);
//
//    EXPECT_EQ(g.size(), 16u);
//}

//todo very strange figure!!!
TEST_F( Simplification,  MFUniquePath ) {
    Graph g(55);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/unique_path", g));
    debruijn_config::simplification::max_flow_ec_remover mfec_config = standard_mfec_config();
    mfec_config.uniqueness_length = 400;
    debruijn::simplification::MaxFlowRemoveErroneousEdges<Graph>(g, mfec_config);

    EXPECT_EQ(12, g.size());
}

//TEST( Simplification,  TopologyTC ) {
//    Graph g(55);
//    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/unique_path", g));
//    debruijn_config::simplification::max_flow_ec_remover tec_config = standard_mfec_config();
//    tec_config.uniqueness_length = 400;
//
//    MaxFlowRemoveErroneousEdges<Graph>(g, tec_config, 0);
//
//    EXPECT_EQ(g.size(), 12u);
//}

//TEST( Simplification,  SelfComp ) {
//       Graph g(55);
//       ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/self_comp", g));
//       debruijn_config::simplification::max_flow_ec_remover mfec_config = standard_mfec_config();
//       mfec_config.uniqueness_length = 1500;
//       MaxFlowRemoveErroneousEdges<Graph>(g, mfec_config);
//
//       EXPECT_EQ(g.size(), 4u);
//}

TEST_F( Simplification,  ComplexBulgeRemoverOnSimpleBulge ) {
    Graph g(55);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/simpliest_bulge/simpliest_bulge", g));
    omnigraph::complex_br::ComplexBulgeRemover<Graph> remover(g, g.k() * 5, 5, nullptr, 1);
    remover.Run();
    INFO("Done");

//       WriteGraphPack(gp, string("./src/test/debruijn/graph_fragments/complex_bulge/complex_bulge_res.dot"));
    EXPECT_EQ(4, g.size());
}

TEST_F( Simplification,  ComplexBulge ) {
    graph_pack::GraphPack gp(55, tmp_folder(), 0);
    ASSERT_TRUE(graphio::ScanGraphPack("./src/test/debruijn/graph_fragments/complex_bulge/complex_bulge", gp));
    auto &graph = gp.get_mutable<Graph>();

    omnigraph::complex_br::ComplexBulgeRemover<Graph> remover(graph, graph.k() * 5, 5, nullptr, 1);
    remover.Run();

    EXPECT_EQ(8, graph.size());
}

TEST_F( Simplification,  BigComplexBulge ) {
    graph_pack::GraphPack gp(55, tmp_folder(), 0);
    ASSERT_TRUE(graphio::ScanGraphPack("./src/test/debruijn/graph_fragments/big_complex_bulge/big_complex_bulge", gp));
    auto &graph = gp.get_mutable<Graph>();

    omnigraph::complex_br::ComplexBulgeRemover<Graph> remover(graph, graph.k() * 5, 5, nullptr, 1);
    remover.Run();
    EXPECT_EQ(66, graph.size());
}

//Relative coverage removal tests

void TestRelativeCoverageRemover(const std::string &path, const std::string &tmp_folder, size_t graph_size) {
    graph_pack::GraphPack gp(55, tmp_folder, 0);
    ASSERT_TRUE(graphio::ScanGraphPack(path, gp));
    INFO("Relative coverage component removal:");
    auto &graph = gp.get_mutable<Graph>();
    auto &flanking_cov = gp.get_mutable<omnigraph::FlankingCoverage<Graph>>();

    auto algo = debruijn::simplification::RelativeCoverageComponentRemoverInstance(graph, flanking_cov,
                                                                                   standard_rcc_config(),
                                                                                   standard_simplif_relevant_info());
    algo->Run();
    EXPECT_EQ(graph_size, graph.size());
}

//todo review
TEST_F( Simplification,  RelativeCoverageRemover ) {
    TestRelativeCoverageRemover(graph_fragment_root() + "rel_cov_ec/constructed_graph", tmp_folder(), 12u);
}

TEST_F( Simplification,  RelativeCoverageRemover1 ) {
    TestRelativeCoverageRemover(graph_fragment_root() + "complex_bulge/complex_bulge", tmp_folder(), 12u);
}

TEST_F( Simplification,  RelativeCoverageRemover2 ) {
    TestRelativeCoverageRemover(graph_fragment_root() + "complex_bulge_2/graph", tmp_folder(), 4u);
}

TEST_F( Simplification,  RelativeCoverageRemover3 ) {
    TestRelativeCoverageRemover(graph_fragment_root() + "tipobulge/tipobulge", tmp_folder(), 4u);
}

TEST_F( Simplification,  RelativeCoverageRemover4 ) {
    TestRelativeCoverageRemover(graph_fragment_root() + "tipobulge_2/graph", tmp_folder(), 4u);
}

//End of relative coverage removal tests


TEST_F( Simplification,  CompressorTest ) {
    std::string path = "./src/test/debruijn/graph_fragments/compression/graph";
    size_t graph_size = 12;
    graph_pack::GraphPack gp(55, tmp_folder(), 0);
    ASSERT_TRUE(graphio::ScanGraphPack(path, gp));
    auto &graph = gp.get_mutable<Graph>();

    CompressAllVertices(graph, standard_simplif_relevant_info().chunk_cnt());
    EXPECT_EQ(graph.size(), graph_size);
}

#if 0
TEST_F( Simplification,  ParallelCompressor1 ) {
    std::string path = "./src/test/debruijn/graph_fragments/compression/graph";
    size_t graph_size = 12;
    graph_pack::GraphPack gp(55, tmp_folder(), 0);
    ASSERT_TRUE(graphio::ScanGraphPack(path, gp));
    auto &graph = gp.get_mutable<Graph>();

    debruijn::simplification::ParallelCompress(graph, standard_simplif_relevant_info().chunk_cnt(), false);
    EXPECT_EQ(graph_size, gp.g.size());
}

TEST_F( Simplification,  ParallelTipClipper1 ) {
    std::string path = "./src/test/debruijn/graph_fragments/tips/graph";
    size_t graph_size = 12;
    graph_pack::GraphPack gp(55, tmp_folder(), 0);
    ASSERT_TRUE(graphio::ScanGraphPack(path, gp));
    auto &graph = gp.get_mutable<Graph>();
    debruijn::simplification::ConditionParser<Graph> parser(graph, standard_tc_config().condition, standard_simplif_relevant_info());
    parser();
    debruijn::simplification::ParallelClipTips(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    EXPECT_EQ(graph_size, gp.g.size());
}

TEST_F( Simplification,  ParallelECRemover ) {
    std::string path = graph_fragment_root() + "complex_bulge/complex_bulge";
    graph_pack::GraphPack gp(55, tmp_folder(), 0);
    ASSERT_TRUE(graphio::ScanGraphPack(path, gp));
    std::string condition = "{ cb 1000 , ec_lb 20 }";
    auto &graph = gp.get_mutable<Graph>();
    debruijn::simplification::ConditionParser<Graph> parser(graph, condition, standard_simplif_relevant_info());
    parser();
    debruijn::simplification::ParallelEC(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    EXPECT_EQ(gp.g.size(), 16u);
    debruijn::simplification::ParallelEC(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    EXPECT_EQ(gp.g.size(), 12u);
}

TEST_F( Simplification,  ParallelECRemover1 ) {
    std::string path = graph_fragment_root() + "complex_bulge_2/graph";
    std::string condition = "{ cb 100 , ec_lb 20 }";
    graph_pack::GraphPack gp(55, tmp_folder(), 0);
    ASSERT_TRUE(graphio::ScanGraphPack(path, gp));
    auto &graph = gp.get_mutable<Graph>();
    debruijn::simplification::ConditionParser<Graph> parser(graph, condition, standard_simplif_relevant_info());
    parser();
    debruijn::simplification::ParallelEC(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    EXPECT_EQ(gp.g.size(), 20u);
    debruijn::simplification::ParallelEC(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    EXPECT_EQ(gp.g.size(), 4u);
    EXPECT_EQ(GraphComponent<Graph>::WholeGraph(gp.g).e_size(), 2u);
}

TEST_F( Simplification,  ParallelECRemover2 ) {
    std::string path = graph_fragment_root() + "rel_cov_ec/constructed_graph";
    std::string condition = "{ cb 100 , ec_lb 20 }";
    graph_pack::GraphPack gp(55, tmp_folder(), 0);
    ASSERT_TRUE(graphio::ScanGraphPack(path, gp));
    auto &graph = gp.get_mutable<Graph>();
    debruijn::simplification::ConditionParser<Graph> parser(graph, condition, standard_simplif_relevant_info());
    parser();
    debruijn::simplification::ParallelEC(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    EXPECT_EQ(gp.g.size(), 20u);
}
#endif

//TEST( Simplification,  ComplexTipRemover ) {
//    string path = "./src/test/debruijn/graph_fragments/ecs/graph";
//    size_t graph_size = 0;
//    GraphPack gp(55, tmp_folder(), 0);
//    graphio::ScanGraphPack(path, gp);
//    ParallelEC(graph, standard_ec_config().condition, standard_simplif_relevant_info());
//    EXPECT_EQ(gp.g.size(), graph_size);
//}
//
//TEST( Simplification,  ComplexTipRemover2 ) {
//    string path = "./src/test/debruijn/graph_fragments/ecs/graph";
//    size_t graph_size = 0;
//    GraphPack gp(55, tmp_folder(), 0);
//    graphio::ScanGraphPack(path, gp);
//    ParallelEC(graph, standard_ec_config().condition, standard_simplif_relevant_info());
//    EXPECT_EQ(gp.g.size(), graph_size);
//}
