//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * path_extend.hpp
 *
 *  Created on: Nov 29, 2013
 *      Author: andrey
 */


#pragma once

#include <boost/test/unit_test.hpp>

#include "test_utils.hpp"
#include "modules/path_extend/path_visualizer.hpp"
#include "modules/path_extend/pe_resolver.hpp"

namespace path_extend {

BOOST_FIXTURE_TEST_SUITE(overlap_removal, fs::TmpFolderFixture)

typedef std::initializer_list<size_t> PathInit;
typedef std::initializer_list<PathInit> PathsInit;

inline BidirectionalPath AsPath(const Graph &g,
                                const GraphElementFinder<Graph> &finder,
                                PathInit ids) {
//    INFO("Converting ids " << utils::join(ids));
    std::vector<EdgeId> edges;
    std::transform(ids.begin(), ids.end(), std::back_inserter(edges),
                   [&] (size_t id) {
                       EdgeId e = finder.ReturnEdgeId(id);
                       VERIFY_MSG(e != EdgeId(), "Couldn't convert id " << id);
//                       INFO("Edge found " << g.str(e));
                       return e;});
    return BidirectionalPath(g, edges);
}

inline void FormPaths(const Graph &g,
               const GraphElementFinder<Graph> &finder,
               GraphCoverageMap &cov_map,
               PathContainer &paths,
               PathsInit path_ids) {
    for (const auto &ids: path_ids) {
        AddPath(paths, AsPath(g, finder, ids), cov_map);
    }
}

inline std::set<size_t> FindPath(const PathContainer &paths, const BidirectionalPath &path) {
    std::set<size_t> answer;
    for (size_t i = 0; i < paths.size(); ++i)
        if (*paths.Get(i) == path || *paths.GetConjugate(i) == path)
            answer.insert(i);
    return answer;
}

inline void CheckPaths(const Graph &g,
               const GraphElementFinder<Graph> &finder,
               const PathContainer &paths,
               PathsInit path_ids) {
//    INFO("Resulting paths");
//    for (size_t i = 0; i < paths.size(); ++i)
//        paths.Get(i)->PrintINFO();

    std::set<size_t> poss;
    for (const auto &ids: path_ids) {
        auto pos = FindPath(paths, AsPath(g, finder, ids));
        BOOST_TEST(!pos.empty());
        utils::insert_all(poss, pos);
    }
    BOOST_TEST(poss.size() == paths.size());
    BOOST_TEST(path_ids.size() == paths.size());
}

inline PathContainer RemoveOverlaps(const Graph &g,
                                    const GraphElementFinder<Graph> &finder,
                                    PathsInit path_ids,
                                    size_t min_edge_len, size_t max_diff,
                                    bool end_start_only = false,
                                    bool retain_one = true) {
    GraphCoverageMap cov_map(g);
    PathContainer container;

    FormPaths(g, finder, cov_map, container, path_ids);

    OverlapRemover overlap_remover(g, container, cov_map,
                                         min_edge_len, max_diff);

    overlap_remover.MarkOverlaps(end_start_only, retain_one);

    PathSplitter splitter(overlap_remover.overlaps(), container, cov_map);
    splitter.Split();

    return container;
}

BOOST_AUTO_TEST_CASE( TestBasicNoDiff ) {
    Graph g(13);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
    GraphCoverageMap cov_map(g);
    PathContainer container;

    auto path1 = AsPath(g, finder, {26, 157, 70, 23, 130, 68});
//    AddPath(container, path1, cov_map);
    auto path2 = AsPath(g, finder, {26, 157, 70, 23, 130, 68});
//    AddPath(container, path2, cov_map);
    auto path3 = AsPath(g, finder, {157, 70, 23, 130, 68});
//    AddPath(container, path2, cov_map);
    auto path4 = AsPath(g, finder, {26, 157, 70, 23, 130});
    size_t min_edge_len = 0;
    size_t max_diff = 0;
    OverlapFindingHelper finding_helper(g, cov_map, min_edge_len, max_diff);
    BOOST_TEST(finding_helper.IsEqual(path1, path2));
    BOOST_TEST(!finding_helper.IsEqual(path1, path3));
    BOOST_TEST(!finding_helper.IsEqual(path1, path4));
    BOOST_TEST(finding_helper.IsSubpath(path3, path1));
    BOOST_TEST(finding_helper.IsSubpath(path4, path1));
    BOOST_TEST(finding_helper.IsSubpath(path2, path1));
    BOOST_TEST(!finding_helper.IsSubpath(path3, path4));
}

BOOST_AUTO_TEST_CASE( TestBasicDiff ) {
    Graph g(13);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
    GraphCoverageMap cov_map(g);
    PathContainer container;

    auto path1 = AsPath(g, finder, {26, 157, 70, 23, 130, 68});
//    AddPath(container, path1, cov_map);
    auto path2 = AsPath(g, finder, {157, 70, 23, 130});
//    AddPath(container, path2, cov_map);
    auto path3 = AsPath(g, finder, {157, 70, 23, 130, 68});
//    AddPath(container, path2, cov_map);
    auto path4 = AsPath(g, finder, {26, 157, 70, 23, 130});
    auto path5 = AsPath(g, finder, {70, 23, 130, 68});
    size_t min_edge_len = 0;
    size_t max_diff = 1;
    OverlapFindingHelper finding_helper(g, cov_map, min_edge_len, max_diff);
    BOOST_TEST(finding_helper.IsEqual(path1, path2));
    BOOST_TEST(finding_helper.IsEqual(path2, path1));

    BOOST_TEST(finding_helper.IsEqual(path1, path3));
    BOOST_TEST(finding_helper.IsEqual(path3, path1));

    BOOST_TEST(finding_helper.IsEqual(path1, path4));
    BOOST_TEST(finding_helper.IsEqual(path4, path1));

    BOOST_TEST(finding_helper.IsSubpath(path2, path1));
    BOOST_TEST(finding_helper.IsSubpath(path1, path2));
    BOOST_TEST(finding_helper.IsSubpath(path3, path1));
    BOOST_TEST(finding_helper.IsSubpath(path1, path3));
    BOOST_TEST(finding_helper.IsSubpath(path4, path1));
    BOOST_TEST(finding_helper.IsSubpath(path1, path4));

    BOOST_TEST(finding_helper.IsSubpath(path3, path4));
    BOOST_TEST(finding_helper.IsSubpath(path4, path3));

    BOOST_TEST(finding_helper.IsSubpath(path5, path1));
    BOOST_TEST(finding_helper.IsSubpath(path5, path2));
    BOOST_TEST(!finding_helper.IsSubpath(path3, path5));
    BOOST_TEST(!finding_helper.IsEqual(path2, path5));
    BOOST_TEST(!finding_helper.IsEqual(path5, path2));
}

BOOST_AUTO_TEST_CASE( TestBasicRepeatConjugateEdges ) {
    Graph g(13);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
    GraphCoverageMap cov_map(g);
    PathContainer container;

    PathsInit path_ids = {{158, 30},
                           {127, 30},
                           {30, 126},
                           {30, 159}};

    PathsInit result_ids = {{30}, {30}, {30}, {30},
                            {158}, {158}, {127}, {127}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff,
                       /*end_start_only*/ true, /*retain one*/ false),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicRepeat ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{17572, 1565},
                     {18066, 1565},
                     {1565,  19391},
                     {1565,  20042}};

    PathsInit result_ids = {{20042}, {19391}, {18066}, {17572},
                            {1565}, {1565}, {1565}, {1565}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff,
                              /*end_start_only*/ true, /*retain one*/ false),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicRepeatUneven ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{17572, 1565},
                          {18066, 1565},
                          {1565,  19391}};

    PathsInit result_ids = {{19391}, {18066}, {17572},
                            {1565}, {1565}, {1565}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff,
                       /*end_start_only*/ true, /*retain one*/ false),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicRepeatAllDisconnected ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{17572, 1565},
                          {18066, 1565},
                          {1565,  19391}};

//    PathsInit result_ids = {{19391}, {18066}, {17572},
//                            {1565}, {1565}, {1565}};
    PathsInit result_ids = {{19391}, {18066}, {17572},
                            {1565}, {1565}, {1565}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff,
                       /*end_start_only*/ false, /*retain one*/ true),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicRepeatNoEffect ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {
            {17572, 1565},
            {18066, 1565}};

    PathsInit result_ids = {{17572, 1565}, {18066, 1565}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff,
                       /*end_start_only*/ true, /*retain one*/ false),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicRepeatEnd ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {
               {17572, 1565},
               {18066, 1565}};

    PathsInit result_ids = {{18066}, {17572}, {1565}, {1565}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff,
                              /*end_start_only*/ false, /*retain one*/ false),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicRepeatEndRetain ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {
            {17572, 1565},
            {18066, 1565}};

    PathsInit result_ids = {{18066, 1565}, {17572}, {1565}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff,
                       /*end_start_only*/ false, /*retain one*/ true),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicRepeatStart ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {
               {1565, 19391},
               {1565, 20042}};

    PathsInit result_ids = {{20042}, {19391}, {1565}, {1565}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff,
                              /*end_start_only*/ false, /*retain one*/ false),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicRepeatStartRetain ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {
            {1565, 19391},
            {1565, 20042}};

    PathsInit result_ids = {{1565, 20042}, {1565}, {19391}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff,
                       /*end_start_only*/ false, /*retain one*/ true),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicDeduplication ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {
            {17572, 1565},
            {17572, 1565},
            {1565, 19391},
            {1565}, {1565}};

    GraphCoverageMap cov_map(g);
    PathContainer container;

    FormPaths(g, finder, cov_map, container, path_ids);

    PathsInit result_ids = {{17572, 1565}, {1565, 19391}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;
    //equal_only = false
    Deduplicate(g, container, cov_map, min_edge_len, max_diff);


    CheckPaths(g, finder,
               container,
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicDeduplicationEqual ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {
            {17572, 1565},
            {17572, 1565},
            {1565, 19391},
            {1565}, {1565}};

    GraphCoverageMap cov_map(g);
    PathContainer container;

    FormPaths(g, finder, cov_map, container, path_ids);

    PathsInit result_ids = {{17572, 1565}, {1565, 19391}, {1565}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;
    Deduplicate(g, container, cov_map, min_edge_len, max_diff, /*equal_only*/true);


    CheckPaths(g, finder,
               container,
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestDeduplicateDiff ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {
            {20129, 17993, 20131, 17993, 19506},
            {20129, 17993, 19506}};

    GraphCoverageMap cov_map(g);
    PathContainer container;

    FormPaths(g, finder, cov_map, container, path_ids);

    //PathsInit result_ids = {{20129, 17993, 20131, 17993, 19506}};
    PathsInit result_ids = {{20129, 17993, 19506}};

    size_t min_edge_len = 0;
    size_t max_diff = 100;
    Deduplicate(g, container, cov_map, min_edge_len, max_diff);

    CheckPaths(g, finder,
               container,
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestDeduplicateDiffGap ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {
            {20129, 17993, 20131, 17993, 19506},
            {20129, 17993, 19506}};

    GraphCoverageMap cov_map(g);
    PathContainer container;

    FormPaths(g, finder, cov_map, container, path_ids);
    container.Get(1)->SetGapAt(2, Gap(100));

    //PathsInit result_ids = {{20129, 17993, 20131, 17993, 19506}};
    PathsInit result_ids = {{20129, 17993, 19506}};

    size_t min_edge_len = 0;
    size_t max_diff = 100;
    Deduplicate(g, container, cov_map, min_edge_len, max_diff);

    CheckPaths(g, finder,
               container,
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestDeduplicateMinEdgeLen1 ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{2142, 4373, 20044},
                          {4373, 20044},
                          {2142, 4373, 20044, 18816}};

    GraphCoverageMap cov_map(g);
    PathContainer container;

    FormPaths(g, finder, cov_map, container, path_ids);

    PathsInit result_ids = {{2142, 4373, 20044},
                            {4373, 20044},
                            {2142, 4373, 20044, 18816}};

    size_t min_edge_len = 130;
    size_t max_diff = 0;
    Deduplicate(g, container, cov_map, min_edge_len, max_diff);

    CheckPaths(g, finder,
               container,
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestDeduplicateMinEdgeLen2 ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{2142, 4373, 20044},
                          {2142, 4373, 20044, 18816}};

    GraphCoverageMap cov_map(g);
    PathContainer container;

    FormPaths(g, finder, cov_map, container, path_ids);

    PathsInit result_ids = {{2142, 4373, 20044, 18816}};

    size_t min_edge_len = 126;
    size_t max_diff = 0;
    Deduplicate(g, container, cov_map, min_edge_len, max_diff);

    CheckPaths(g, finder,
               container,
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestEndEndMinEdgeLen ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{2142, 4373, 20044, 20385},
                          {2142, 4373, 20044, 18816}};

    PathsInit result_ids = {{2142, 4373, 20044}, {20385},
                            {2142, 4373, 20044, 18816}};

    size_t min_edge_len = 126;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestTricky0 ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{2142, 4373, 20044, 20385},
                          {4373, 20044, 18816}};

    PathsInit result_ids = {{2142, 4373, 20044, 20385},
                            {4373, 20044, 18816}};

    size_t min_edge_len = 127;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestTricky1 ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{2142, 4373, 20044, 20385},
                          {4373, 20044, 18816}};

    PathsInit result_ids = {{2142, 4373, 20044, 20385},
                            {4373, 20044}, {18816}};

    size_t min_edge_len = 126;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestTricky2 ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{2142, 4373, 20044, 20385},
                          {4373, 20044, 18816}};

    PathsInit result_ids = {{2142, 4373, 20044}, {20385},
                            {4373, 20044, 18816}};

    size_t min_edge_len = 126;
    size_t max_diff = 2;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestBasicMinEdgeLen ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{20449, 2142, 4373, 20044},
                          {2142, 4373, 20044, 20385}};

    PathsInit result_ids = {{20449}, {2142, 4373, 20044},
                            {2142, 4373, 20044}, {20385}};

    size_t min_edge_len = 100;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff),
               result_ids);
}

BOOST_AUTO_TEST_CASE( TestOverlapWithMiddle ) {
    Graph g(55);
    GraphElementFinder<Graph> finder(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g);

    PathsInit path_ids = {{18953, 20449, 20129},
                          {20385, 20449},
                          {20449, 2142}};

    PathsInit result_ids = {{18953, 20449, 20129},
                            {20385}, {20449},
                            {20449}, {2142}};

    size_t min_edge_len = 0;
    size_t max_diff = 0;

    CheckPaths(g, finder,
               RemoveOverlaps(g, finder, path_ids,
                              min_edge_len, max_diff),
               result_ids);
}

//TODO add more tricky tests on whole the process

BOOST_AUTO_TEST_SUITE_END()

}
