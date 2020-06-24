//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "graphio.hpp"
#include "test_utils.hpp"
//#include "modules/simplification/parallel_simplification_algorithms.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "stages/simplification_pipeline/single_cell_simplification.hpp"
#include "stages/simplification_pipeline/rna_simplification.hpp"
#include <boost/test/unit_test.hpp>
//#include "repeat_resolving_routine.hpp"

namespace debruijn_graph {
using namespace config;

BOOST_FIXTURE_TEST_SUITE(graph_simplification_tests, fs::TmpFolderFixture)

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

BOOST_AUTO_TEST_CASE( SimpleTipClipperTest ) {
    ConjugateDeBruijnGraph g(55);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/simpliest_tip/simpliest_tip", g);

    DefaultClipTips(g);

    BOOST_CHECK_EQUAL(g.size(), 4u);
}

BOOST_AUTO_TEST_CASE( SimpleBulgeRemovalTest ) {
    Graph g(55);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/simpliest_bulge/simpliest_bulge", g);

    DefaultRemoveBulges(g);

    BOOST_CHECK_EQUAL(g.size(), 4u);
}

BOOST_AUTO_TEST_CASE( TipobulgeTest ) {
    Graph g(55);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/tipobulge/tipobulge", g);

    DefaultClipTips(g);

    DefaultRemoveBulges(g);

    BOOST_CHECK_EQUAL(g.size(), 16u);
}

BOOST_AUTO_TEST_CASE( SimpleECTest ) {
    Graph g(55);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g);

    debruijn_config::simplification::erroneous_connections_remover ec_config;
    ec_config.condition = "{ icb 7000 , ec_lb 20 }";

    debruijn::simplification::ECRemoverInstance(g, ec_config, standard_simplif_relevant_info())->Run();

    BOOST_CHECK_EQUAL(g.size(), 16u);
}

BOOST_AUTO_TEST_CASE( SimpleIterECTest ) {
    Graph g(55);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g);

    debruijn_config::simplification::erroneous_connections_remover ec_config;
    ec_config.condition = "{ icb 7000 , ec_lb 20 }";

    auto ec_remover_ptr = debruijn::simplification::ECRemoverInstance(g, ec_config, standard_simplif_relevant_info());

    AlgorithmRunningHelper<Graph>::IterativeThresholdsRun(*ec_remover_ptr, 2);
    BOOST_CHECK_EQUAL(g.size(), 16u);
}

BOOST_AUTO_TEST_CASE( IterECTest ) {
    Graph g(55);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g);

    debruijn_config::simplification::erroneous_connections_remover ec_config;
    ec_config.condition = "{ icb 7000 , ec_lb 20 }";

    auto ec_remover_ptr = debruijn::simplification::ECRemoverInstance(g, ec_config, standard_simplif_relevant_info());

    ec_remover_ptr->Run(false, 0.5);
    BOOST_CHECK_EQUAL(g.size(), 20u);

    ec_remover_ptr->Run(false, 1.0);
    BOOST_CHECK_EQUAL(g.size(), 16u);
}

BOOST_AUTO_TEST_CASE( IterUniquePath ) {
    Graph g(55);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g);

    auto tec_config = standard_tec_config();
    while(debruijn::simplification::TopologyRemoveErroneousEdges<Graph>(g, tec_config, 0)) {
    }

    BOOST_CHECK_EQUAL(g.size(), 16u);
}

//BOOST_AUTO_TEST_CASE( MFIterUniquePath ) {
//    Graph g(55);
//    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g);
//
//    debruijn_config::simplification::max_flow_ec_remover mfec_config = standard_mfec_config();
//    mfec_config.uniqueness_length = 500;
//    MaxFlowRemoveErroneousEdges<Graph>(g, mfec_config);
//
//    BOOST_CHECK_EQUAL(g.size(), 16u);
//}

//todo very strange figure!!!
BOOST_AUTO_TEST_CASE( MFUniquePath ) {
    Graph g(55);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/unique_path", g);
    debruijn_config::simplification::max_flow_ec_remover mfec_config = standard_mfec_config();
    mfec_config.uniqueness_length = 400;
    debruijn::simplification::MaxFlowRemoveErroneousEdges<Graph>(g, mfec_config);

    BOOST_CHECK_EQUAL(g.size(), 12u);
}

//BOOST_AUTO_TEST_CASE( TopologyTC ) {
//    Graph g(55);
//    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/unique_path", g);
//    debruijn_config::simplification::max_flow_ec_remover tec_config = standard_mfec_config();
//    tec_config.uniqueness_length = 400;
//
//    MaxFlowRemoveErroneousEdges<Graph>(g, tec_config, 0);
//
//    BOOST_CHECK_EQUAL(g.size(), 12u);
//}

//BOOST_AUTO_TEST_CASE( SelfComp ) {
//       Graph g(55);
//       graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/self_comp", g);
//       debruijn_config::simplification::max_flow_ec_remover mfec_config = standard_mfec_config();
//       mfec_config.uniqueness_length = 1500;
//       MaxFlowRemoveErroneousEdges<Graph>(g, mfec_config);
//
//       BOOST_CHECK_EQUAL(g.size(), 4u);
//}

BOOST_AUTO_TEST_CASE( ComplexBulgeRemoverOnSimpleBulge ) {
       Graph g(55);
       graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/simpliest_bulge/simpliest_bulge", g);
       omnigraph::complex_br::ComplexBulgeRemover<Graph> remover(g, g.k() * 5, 5, nullptr, 1);
       remover.Run();
       INFO("Done");

//       WriteGraphPack(gp, string("./src/test/debruijn/graph_fragments/complex_bulge/complex_bulge_res.dot"));
       BOOST_CHECK_EQUAL(g.size(), 4u);
}

BOOST_AUTO_TEST_CASE( ComplexBulge ) {
    GraphPack gp(55, "tmp", 0);
    graphio::ScanGraphPack("./src/test/debruijn/graph_fragments/complex_bulge/complex_bulge", gp);
    auto &graph = gp.get_mutable<Graph>();
    omnigraph::complex_br::ComplexBulgeRemover<Graph> remover(graph, graph.k() * 5, 5, nullptr, 1);
    remover.Run();
    BOOST_CHECK_EQUAL(graph.size(), 8u);
}

BOOST_AUTO_TEST_CASE( BigComplexBulge ) {
    GraphPack gp(55, "tmp", 0);
    graphio::ScanGraphPack("./src/test/debruijn/graph_fragments/big_complex_bulge/big_complex_bulge", gp);
    auto &graph = gp.get_mutable<Graph>();
    omnigraph::complex_br::ComplexBulgeRemover<Graph> remover(graph, graph.k() * 5, 5, nullptr, 1);
    remover.Run();
    BOOST_CHECK_EQUAL(graph.size(), 66u);
}

//Relative coverage removal tests

void TestRelativeCoverageRemover(const std::string &path, size_t graph_size) {
    GraphPack gp(55, "tmp", 0);
    graphio::ScanGraphPack(path, gp);
    INFO("Relative coverage component removal:");
    auto &graph = gp.get_mutable<Graph>();
    auto &index = gp.get_mutable<EdgeIndex<Graph>>();
    auto &flanking_cov = gp.get_mutable<omnigraph::FlankingCoverage<Graph>>();

    auto algo = debruijn::simplification::RelativeCoverageComponentRemoverInstance(graph, flanking_cov,
                                                                                   standard_rcc_config(),
                                                                                   standard_simplif_relevant_info());
    algo->Run();
    BOOST_CHECK_EQUAL(graph.size(), graph_size);
}

//todo review
BOOST_AUTO_TEST_CASE( RelativeCoverageRemover ) {
    TestRelativeCoverageRemover(graph_fragment_root() + "rel_cov_ec/constructed_graph", 12u);
}

BOOST_AUTO_TEST_CASE( RelativeCoverageRemover1 ) {
    TestRelativeCoverageRemover(graph_fragment_root() + "complex_bulge/complex_bulge", 12u);
}

BOOST_AUTO_TEST_CASE( RelativeCoverageRemover2 ) {
    TestRelativeCoverageRemover(graph_fragment_root() + "complex_bulge_2/graph", 4u);
}

BOOST_AUTO_TEST_CASE( RelativeCoverageRemover3 ) {
    TestRelativeCoverageRemover(graph_fragment_root() + "tipobulge/tipobulge", 4u);
}

BOOST_AUTO_TEST_CASE( RelativeCoverageRemover4 ) {
    TestRelativeCoverageRemover(graph_fragment_root() + "tipobulge_2/graph", 4u);
}

//End of relative coverage removal tests


BOOST_AUTO_TEST_CASE( CompressorTest ) {
    std::string path = "./src/test/debruijn/graph_fragments/compression/graph";
    size_t graph_size = 12;
    GraphPack gp(55, "tmp", 0);
    graphio::ScanGraphPack(path, gp);
    auto &graph = gp.get_mutable<Graph>();
    CompressAllVertices(graph, standard_simplif_relevant_info().chunk_cnt());
    BOOST_CHECK_EQUAL(graph.size(), graph_size);
}

#if 0
BOOST_AUTO_TEST_CASE( ParallelCompressor1 ) {
    std::string path = "./src/test/debruijn/graph_fragments/compression/graph";
    size_t graph_size = 12;
    GraphPack gp(55, "tmp", 0);
    graphio::ScanGraphPack(path, gp);
    auto &graph = gp.get_mutable<Graph>();
    debruijn::simplification::ParallelCompress(graph, standard_simplif_relevant_info().chunk_cnt(), false);
    BOOST_CHECK_EQUAL(graph.size(), graph_size);
}

BOOST_AUTO_TEST_CASE( ParallelTipClipper1 ) {
    std::string path = "./src/test/debruijn/graph_fragments/tips/graph";
    size_t graph_size = 12;
    GraphPack gp(55, "tmp", 0);
    graphio::ScanGraphPack(path, gp);
    auto &graph = gp.get_mutable<Graph>();
    debruijn::simplification::ConditionParser<Graph> parser(graph, standard_tc_config().condition, standard_simplif_relevant_info());
    parser();
    debruijn::simplification::ParallelClipTips(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    BOOST_CHECK_EQUAL(graph.size(), graph_size);
}

BOOST_AUTO_TEST_CASE( ParallelECRemover ) {
    std::string path = graph_fragment_root() + "complex_bulge/complex_bulge";
    GraphPack gp(55, "tmp", 0);
    graphio::ScanGraphPack(path, gp);
    std::string condition = "{ cb 1000 , ec_lb 20 }";
    auto &graph = gp.get_mutable<Graph>();
    debruijn::simplification::ConditionParser<Graph> parser(graph, condition, standard_simplif_relevant_info());
    parser();
    debruijn::simplification::ParallelEC(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    BOOST_CHECK_EQUAL(graph.size(), 16u);
    debruijn::simplification::ParallelEC(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    BOOST_CHECK_EQUAL(graph.size(), 12u);
}

BOOST_AUTO_TEST_CASE( ParallelECRemover1 ) {
    std::string path = graph_fragment_root() + "complex_bulge_2/graph";
    std::string condition = "{ cb 100 , ec_lb 20 }";
    GraphPack gp(55, "tmp", 0);
    graphio::ScanGraphPack(path, gp);
    auto &graph = gp.get_mutable<Graph>();
    debruijn::simplification::ConditionParser<Graph> parser(graph, condition, standard_simplif_relevant_info());
    parser();
    debruijn::simplification::ParallelEC(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    BOOST_CHECK_EQUAL(graph.size(), 20u);
    debruijn::simplification::ParallelEC(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    BOOST_CHECK_EQUAL(graph.size(), 4u);
    BOOST_CHECK_EQUAL(GraphComponent<Graph>::WholeGraph(graph).e_size(), 2u);
}

BOOST_AUTO_TEST_CASE( ParallelECRemover2 ) {
    std::string path = graph_fragment_root() + "rel_cov_ec/constructed_graph";
    std::string condition = "{ cb 100 , ec_lb 20 }";
    GraphPack gp(55, "tmp", 0);
    graphio::ScanGraphPack(path, gp);
    auto &graph = gp.get_mutable<Graph>();
    debruijn::simplification::ConditionParser<Graph> parser(graph, condition, standard_simplif_relevant_info());
    parser();
    debruijn::simplification::ParallelEC(graph, parser.max_length_bound(), parser.max_coverage_bound(), standard_simplif_relevant_info().chunk_cnt());
    BOOST_CHECK_EQUAL(graph.size(), 20u);
}
#endif

//BOOST_AUTO_TEST_CASE( ComplexTipRemover ) {
//    string path = "./src/test/debruijn/graph_fragments/ecs/graph";
//    size_t graph_size = 0;
//    conj_graph_pack gp(55, "tmp", 0);
//    graphio::ScanGraphPack(path, gp);
//    ParallelEC(gp.g, standard_ec_config().condition, standard_simplif_relevant_info());
//    BOOST_CHECK_EQUAL(gp.g.size(), graph_size);
//}
//
//BOOST_AUTO_TEST_CASE( ComplexTipRemover2 ) {
//    string path = "./src/test/debruijn/graph_fragments/ecs/graph";
//    size_t graph_size = 0;
//    conj_graph_pack gp(55, "tmp", 0);
//    graphio::ScanGraphPack(path, gp);
//    ParallelEC(gp.g, standard_ec_config().condition, standard_simplif_relevant_info());
//    BOOST_CHECK_EQUAL(gp.g.size(), graph_size);
//}

BOOST_AUTO_TEST_SUITE_END()}
