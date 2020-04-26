//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/test/unit_test.hpp>

#include "standard.hpp"
#include "pipeline/graph_pack.hpp"
#include "pipeline/graphio.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"

namespace debruijn_graph {

BOOST_FIXTURE_TEST_SUITE(detail_coverage_tests, TmpFolderFixture)

BOOST_AUTO_TEST_CASE( AgreementTest ) {
    string filename = "./data/debruijn/ECOLI_IS220_QUAKE_400K/K55/latest/saves/construction";

    GraphPack gp(55, "tmp", 0,
                       Sequence(), 0,
                       false, false,
                       /*flanking_range*/50);
    graphio::ScanBasicGraph(filename, gp.g);
    bool loaded = graphio::LoadEdgeIndex(filename, gp.index.inner_index());
    VERIFY(loaded);
    gp.index.Update();
    gp.index.Attach();
    FlankingCoverage<Graph, Index::InnerIndexT> old_flanking(gp.g, gp.index.inner_index(), 50);
    gp.flanking_cov.Fill(gp.index.inner_index());
    for (auto it = gp.g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        double tolerance = 1.;//one percent
        BOOST_CHECK_CLOSE(gp.flanking_cov.GetInCov(e), old_flanking.GetInCov(e), tolerance);
        BOOST_CHECK_CLOSE(gp.flanking_cov.GetOutCov(e), old_flanking.GetOutCov(e), tolerance);
//        cout << gp.flanking_cov.GetInCov(e) << endl;
//        BOOST_CHECK(math::eq(gp.flanking_cov.GetOutCov(e), old_flanking.GetOutCov(e)));
//        BOOST_CHECK(math::eq(gp.flanking_cov.GetInCov(e), old_flanking.GetInCov(e)));
    }
}

BOOST_AUTO_TEST_CASE( AgreementTest2 ) {
    string filename = "./data/debruijn/ECOLI_IS220_QUAKE_400K/K55/latest/saves/simplification";

    GraphPack gp(55, "tmp", 0,
                       Sequence(), 0,
                       false, false,
                       /*flanking_range*/50);
    graphio::ScanBasicGraph(filename, gp.g);
    bool loaded = graphio::LoadEdgeIndex("./data/debruijn/ECOLI_IS220_QUAKE_400K/K55/latest/saves/construction", gp.index.inner_index());
    VERIFY(loaded);
    gp.index.Update();
    gp.index.Attach();
    FlankingCoverage<Graph, Index::InnerIndexT> old_flanking(gp.g, gp.index.inner_index(), 50);
    gp.flanking_cov.Fill(gp.index.inner_index());
    for (auto it = gp.g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        double tolerance = 1.;//one percent
        cout << "diff " << (double) gp.flanking_cov.GetInCov(e) - (double) old_flanking.GetInCov(e) << " ; edge "
                << gp.g.str(e) << " ; abs new " << (double) gp.flanking_cov.GetInCov(e) << endl;
//        BOOST_CHECK_CLOSE(gp.flanking_cov.GetInCov(e), old_flanking.GetInCov(e), tolerance);
//        BOOST_CHECK_CLOSE(gp.flanking_cov.GetOutCov(e), old_flanking.GetOutCov(e), tolerance);
//        cout << gp.flanking_cov.GetInCov(e) << endl;
//        BOOST_CHECK(math::eq(gp.flanking_cov.GetOutCov(e), old_flanking.GetOutCov(e)));
//        BOOST_CHECK(math::eq(gp.flanking_cov.GetInCov(e), old_flanking.GetInCov(e)));
    }
}

BOOST_AUTO_TEST_SUITE_END()

}
