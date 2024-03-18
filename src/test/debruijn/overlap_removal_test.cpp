//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/handlers/id_track_handler.hpp"
#include "modules/path_extend/pe_resolver.hpp"
#include "modules/path_extend/pe_utils.hpp"
#include "modules/path_extend/path_extender.hpp"
#include "modules/path_extend/overlap_remover.hpp"

#include "graphio.hpp"

#include <gtest/gtest.h>

using namespace path_extend;
using namespace debruijn_graph;

typedef std::initializer_list<size_t> PathInit;
typedef std::initializer_list<PathInit> PathsInit;

std::unique_ptr<BidirectionalPath> AsPath(const Graph &g,
                                          const omnigraph::GraphElementFinder<Graph> &finder,
                                          PathInit ids) {
//    INFO("Converting ids " << utils::join(ids));
    std::vector<EdgeId> edges;
    std::transform(ids.begin(), ids.end(), std::back_inserter(edges),
                   [&] (size_t id) {
                       EdgeId e = finder.ReturnEdgeId(id);
                       VERIFY_MSG(e != EdgeId(), "Couldn't convert id " << id);
//                       INFO("Edge found " << g.str(e));
                       return e;});
    return BidirectionalPath::create(g, edges);
}

inline void FormPaths(const Graph &g,
               const omnigraph::GraphElementFinder<Graph> &finder,
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
        if (paths.Get(i) == path || paths.GetConjugate(i) == path)
            answer.insert(i);
    return answer;
}

inline void CheckPaths(const Graph &g,
               const omnigraph::GraphElementFinder<Graph> &finder,
               const PathContainer &paths,
               PathsInit path_ids) {
//    INFO("Resulting paths");
//    for (size_t i = 0; i < paths.size(); ++i)
//        paths.Get(i)->PrintINFO();

    std::set<size_t> poss;
    for (const auto &ids: path_ids) {
        auto path = AsPath(g, finder, ids);
        auto pos = FindPath(paths, *path);
        EXPECT_FALSE(pos.empty());
        utils::insert_all(poss, pos);
    }
    EXPECT_EQ(poss.size(), paths.size());
    EXPECT_EQ(path_ids.size(), paths.size());
}

inline PathContainer RemoveOverlaps(const Graph &g,
                                    const omnigraph::GraphElementFinder<Graph> &finder,
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

TEST( OverlapRemoval, BasicNoDiff ) {
    Graph g(13);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
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
    EXPECT_TRUE(finding_helper.IsEqual(*path1, *path2));
    EXPECT_FALSE(finding_helper.IsEqual(*path1, *path3));
    EXPECT_FALSE(finding_helper.IsEqual(*path1, *path4));
    EXPECT_TRUE(finding_helper.IsSubpath(*path3, *path1));
    EXPECT_TRUE(finding_helper.IsSubpath(*path4, *path1));
    EXPECT_TRUE(finding_helper.IsSubpath(*path2, *path1));
    EXPECT_FALSE(finding_helper.IsSubpath(*path3, *path4));
}

TEST( OverlapRemoval, BasicDiff ) {
    Graph g(13);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
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
    EXPECT_TRUE(finding_helper.IsEqual(*path1, *path2));
    EXPECT_TRUE(finding_helper.IsEqual(*path2, *path1));

    EXPECT_TRUE(finding_helper.IsEqual(*path1, *path3));
    EXPECT_TRUE(finding_helper.IsEqual(*path3, *path1));

    EXPECT_TRUE(finding_helper.IsEqual(*path1, *path4));
    EXPECT_TRUE(finding_helper.IsEqual(*path4, *path1));

    EXPECT_TRUE(finding_helper.IsSubpath(*path2, *path1));
    EXPECT_TRUE(finding_helper.IsSubpath(*path1, *path2));
    EXPECT_TRUE(finding_helper.IsSubpath(*path3, *path1));
    EXPECT_TRUE(finding_helper.IsSubpath(*path1, *path3));
    EXPECT_TRUE(finding_helper.IsSubpath(*path4, *path1));
    EXPECT_TRUE(finding_helper.IsSubpath(*path1, *path4));

    EXPECT_TRUE(finding_helper.IsSubpath(*path3, *path4));
    EXPECT_TRUE(finding_helper.IsSubpath(*path4, *path3));

    EXPECT_TRUE(finding_helper.IsSubpath(*path5, *path1));
    EXPECT_TRUE(finding_helper.IsSubpath(*path5, *path2));
    EXPECT_FALSE(finding_helper.IsSubpath(*path3, *path5));
    EXPECT_FALSE(finding_helper.IsEqual(*path2, *path5));
    EXPECT_FALSE(finding_helper.IsEqual(*path5, *path2));
}

TEST( OverlapRemoval, BasicRepeatConjugateEdges ) {
    Graph g(13);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
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

TEST( OverlapRemoval, BasicRepeat ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, BasicRepeatUneven ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, BasicRepeatAllDisconnected ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, BasicRepeatNoEffect ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, BasicRepeatEnd ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, BasicRepeatEndRetain ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, BasicRepeatStart ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, BasicRepeatStartRetain ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, BasicDeduplication ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, BasicDeduplicationEqual ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, DeduplicateDiff ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, DeduplicateDiffGap ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

    PathsInit path_ids = {
            {20129, 17993, 20131, 17993, 19506},
            {20129, 17993, 19506}};

    GraphCoverageMap cov_map(g);
    PathContainer container;

    FormPaths(g, finder, cov_map, container, path_ids);
    container.Get(1).SetGapAt(2, Gap(100));

    //PathsInit result_ids = {{20129, 17993, 20131, 17993, 19506}};
    PathsInit result_ids = {{20129, 17993, 19506}};

    size_t min_edge_len = 0;
    size_t max_diff = 100;
    Deduplicate(g, container, cov_map, min_edge_len, max_diff);

    CheckPaths(g, finder,
               container,
               result_ids);
}

TEST( OverlapRemoval, DeduplicateMinEdgeLen1 ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, DeduplicateMinEdgeLen2 ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, EndEndMinEdgeLen ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, Tricky0 ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, Tricky1 ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, Tricky2 ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, BasicMinEdgeLen ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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

TEST( OverlapRemoval, OverlapWithMiddle ) {
    Graph g(55);
    omnigraph::GraphElementFinder<Graph> finder(g);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/ecoli_400k/distance_estimation", g));

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
