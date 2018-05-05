#include <gtest/gtest.h>
#include "graph.hpp"
#include "hmmpath.hpp"

TEST(DepthAtLeast, Depth) {
  for (size_t n = 10; n <= 10000000; n *= 10) {
    std::string s(n, 'A');
    auto graph = Graph({s});
    impl::DepthAtLeast<Graph::GraphCursor> depthAtLeast;
    auto cursor = graph.begins()[0];
    EXPECT_TRUE(depthAtLeast.depth_at_least(cursor, std::min<size_t>(n, depthAtLeast.STACK_LIMIT / 2)));
    if (n < depthAtLeast.STACK_LIMIT / 2) {
      EXPECT_FALSE(depthAtLeast.depth_at_least(cursor, n + 1));
    }
  }
}

TEST(TestStackLimit, Depth) {
  for (size_t n = 10; n <= 10000000; n *= 10) {
    std::string s(n, 'A');
    auto graph = Graph({s});
    impl::Depth<Graph::GraphCursor> depth;
    auto cursor = graph.begins()[0];
    EXPECT_EQ(depth.depth(cursor), n);
  }
}

