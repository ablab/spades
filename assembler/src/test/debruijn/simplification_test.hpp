//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "graph_simplification.hpp"
//#include "repeat_resolving_routine.hpp"

namespace debruijn_graph {

BOOST_FIXTURE_TEST_SUITE(graph_simplification_tests, TmpFolderFixture)

static debruijn_config::simplification::bulge_remover standard_br_config_generation() {
	debruijn_config::simplification::bulge_remover br_config;
	br_config.enabled = true;
	br_config.max_bulge_length_coefficient = 4;
	br_config.max_additive_length_coefficient = 0;
	br_config.max_coverage = 1000.;
	br_config.max_relative_coverage = 1.2;
	br_config.max_delta = 3;
	br_config.max_relative_delta = 0.1;
	return br_config;
}

static size_t standard_read_length() {
	return 100;
}

debruijn_config::simplification::bulge_remover standard_br_config() {
	static debruijn_config::simplification::bulge_remover br_config = standard_br_config_generation();
	return br_config;
}

static debruijn_config::simplification::erroneous_connections_remover standard_ec_config_generation() {
	debruijn_config::simplification::erroneous_connections_remover ec_config;
	ec_config.condition = "{ cb 30 , ec_lb 20 }";
	return ec_config;
}

debruijn_config::simplification::erroneous_connections_remover standard_ec_config() {
	static debruijn_config::simplification::erroneous_connections_remover ec_config = standard_ec_config_generation();
	return ec_config;
}

static debruijn_config::simplification::topology_based_ec_remover topology_based_ec_config_generation() {
	debruijn_config::simplification::topology_based_ec_remover tec_config;
	tec_config.max_ec_length_coefficient = 20;
	tec_config.plausibility_length = 200;
	tec_config.uniqueness_length = 1500;
	return tec_config;
}

static debruijn_config::simplification::max_flow_ec_remover max_flow_based_ec_config_generation() {
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

static debruijn_config::simplification::tip_clipper standard_tc_config_generation() {
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

void PrintGraph(const Graph & g) {
	FOREACH(VertexId v, g.vertices()) {
		FOREACH(EdgeId e, g.OutgoingEdges(v)) {
			cout << g.int_id(e) << ":" << g.int_id(g.EdgeStart(e)) << " " << g.int_id(g.EdgeEnd(e)) << endl;
		}
	}
	cout << endl;
}

void DefaultClipTips(Graph& graph) {
	ClipTips(graph, standard_tc_config(), 100);
}
/*
void DefaultRemoveBulges(Graph& graph) {
	auto factory = GetBulgeRemoverFactory(graph, standard_br_config());
	RunConcurrentAlgorithm(graph, factory, CoverageComparator<Graph>(graph));
}
*/

BOOST_AUTO_TEST_CASE( SimpleTipClipperTest ) {
    ConjugateDeBruijnGraph g(55);
	graphio::ScanBasicGraph<ConjugateDeBruijnGraph>("./src/test/debruijn/graph_fragments/simpliest_tip/simpliest_tip", g);

	DefaultClipTips(g);

	BOOST_CHECK_EQUAL(g.size(), 4u);
}

BOOST_AUTO_TEST_CASE( SimpleBulgeRemovalTest ) {
	Graph g(55);
	graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/simpliest_bulge/simpliest_bulge", g);

	RemoveBulges(g, standard_br_config());

	BOOST_CHECK_EQUAL(g.size(), 4u);
}

BOOST_AUTO_TEST_CASE( TipobulgeTest ) {
	Graph g(55);
	graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/tipobulge/tipobulge", g);

	DefaultClipTips(g);

	RemoveBulges(g, standard_br_config());

	BOOST_CHECK_EQUAL(g.size(), 16u);
}

BOOST_AUTO_TEST_CASE( SimpleECTest ) {
	Graph g(55);
	graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g);

	debruijn_config::simplification::erroneous_connections_remover ec_config;
	ec_config.condition = "{ icb 7000 , ec_lb 20 }";

	RemoveLowCoverageEdges(g, ec_config);

	BOOST_CHECK_EQUAL(g.size(), 16u);
}

BOOST_AUTO_TEST_CASE( IterECTest ) {
	Graph g(55);
	graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g);

	debruijn_config::simplification::erroneous_connections_remover ec_config;
	ec_config.condition = "{ icb 7000 , ec_lb 20 }";

	RemoveLowCoverageEdges<Graph>(g, ec_config, 0, 0, 0., 2, 0);

	BOOST_CHECK_EQUAL(g.size(), 20u);

	RemoveLowCoverageEdges<Graph>(g, ec_config, 0, 0, 0., 2, 1);

	BOOST_CHECK_EQUAL(g.size(), 16u);
}

BOOST_AUTO_TEST_CASE( IterUniquePath ) {
	Graph g(55);
	graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g);

	debruijn_config::simplification::topology_based_ec_remover tec_config = standard_tec_config();
	while(TopologyRemoveErroneousEdges<Graph>(g, tec_config, 0)) {
	}

	BOOST_CHECK_EQUAL(g.size(), 16u);
}

//BOOST_AUTO_TEST_CASE( MFIterUniquePath ) {
//	Graph g(55);
//	graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path", g);
//
//	debruijn_config::simplification::max_flow_ec_remover mfec_config = standard_mfec_config();
//	mfec_config.uniqueness_length = 500;
//	MaxFlowRemoveErroneousEdges<Graph>(g, mfec_config);
//
//	BOOST_CHECK_EQUAL(g.size(), 16u);
//}

//todo very strange figure!!!
BOOST_AUTO_TEST_CASE( MFUniquePath ) {
	Graph g(55);
	graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/unique_path", g);
	debruijn_config::simplification::max_flow_ec_remover mfec_config = standard_mfec_config();
	mfec_config.uniqueness_length = 400;
	MaxFlowRemoveErroneousEdges<Graph>(g, mfec_config);

	BOOST_CHECK_EQUAL(g.size(), 12u);
}

//BOOST_AUTO_TEST_CASE( TopologyTC ) {
//	Graph g(55);
//	graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/topology_ec/unique_path", g);
//	debruijn_config::simplification::max_flow_ec_remover tec_config = standard_mfec_config();
//	tec_config.uniqueness_length = 400;
//
//	MaxFlowRemoveErroneousEdges<Graph>(g, tec_config, 0);
//
//	BOOST_CHECK_EQUAL(g.size(), 12u);
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
//       OppositionLicvidator<Graph> licvidator(gp.g, gp.g.k() * 5, 5);
//       licvidator.Licvidate();
       omnigraph::complex_br::ComplexBulgeRemover<Graph> remover(g, g.k() * 5, 5);
       remover.Run();
       INFO("Done");

//       WriteGraphPack(gp, string("./src/test/debruijn/graph_fragments/complex_bulge/complex_bulge_res.dot"));
       BOOST_CHECK_EQUAL(g.size(), 4u);
}

BOOST_AUTO_TEST_CASE( ComplexBulge ) {
    conj_graph_pack gp(55, tmp_folder, 0);
       graphio::ScanGraphPack("./src/test/debruijn/graph_fragments/complex_bulge/complex_bulge", gp);
//       OppositionLicvidator<Graph> licvidator(gp.g, gp.g.k() * 5, 5);
//       licvidator.Licvidate();

       omnigraph::complex_br::ComplexBulgeRemover<Graph> remover(gp.g, gp.g.k() * 5, 5);
       remover.Run();

//       WriteGraphPack(gp, string("./src/test/debruijn/graph_fragments/complex_bulge/complex_bulge_res.dot"));
       BOOST_CHECK_EQUAL(gp.g.size(), 8u);
}

BOOST_AUTO_TEST_CASE( BigComplexBulge ) {
    conj_graph_pack gp(55, tmp_folder, 0);
       graphio::ScanGraphPack("./src/test/debruijn/graph_fragments/big_complex_bulge/big_complex_bulge", gp);
//       OppositionLicvidator<Graph> licvidator(gp.g, gp.g.k() * 5, 5);
//       licvidator.Licvidate();
       omnigraph::complex_br::ComplexBulgeRemover<Graph> remover(gp.g, gp.g.k() * 5, 5);
       remover.Run();
//       WriteGraphPack(gp, string("./src/test/debruijn/graph_fragments/big_complex_bulge/big_complex_bulge_res.dot"));
       BOOST_CHECK_EQUAL(gp.g.size(), 66u);
}

template<class Graph, class InnerIndex>
void FillKmerCoverageWithAvg(const Graph& g, InnerIndex& idx) {
    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        Sequence nucls = g.EdgeNucls(e);
        double cov = g.coverage(e);
        auto kpomer = idx.ConstructKWH(runtime_k::RtSeq(g.k() + 1, nucls));
        kpomer >>= 0;
        for (size_t i = 0; i < g.length(e); ++i) {
            kpomer <<= nucls[i + g.k()];
            idx.get_raw_value_reference(kpomer).count = unsigned(math::floor(cov));
        }
    }
}

BOOST_AUTO_TEST_CASE( RelativeCoverageRemover ) {
    typedef graph_pack<ConjugateDeBruijnGraph, runtime_k::RtSeq, 
            KmerStoringEdgeIndex<ConjugateDeBruijnGraph, runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq>, SimpleStoring>> gp_t;
    gp_t gp(55, tmp_folder, 0);
    graphio::ScanGraphPack("./src/test/debruijn/graph_fragments/rel_cov_ec/constructed_graph", gp);
    INFO("Relative coverage component removal:");
    FillKmerCoverageWithAvg(gp.g, gp.index.inner_index());
    gp.flanking_cov.Fill(gp.index.inner_index());
    RemoveRelativelyLowCoverageComponents(gp.g, gp.flanking_cov, standard_rcc_config(), 10., 100);
    BOOST_CHECK_EQUAL(gp.g.size(), 12u/*28u*/);
}

BOOST_AUTO_TEST_SUITE_END()}
