#include <gtest/gtest.h>
#include "graph.hpp"
#include "hmmpath.hpp"

TEST(TestStackLimit, Depth) {
  for (size_t n = 10; n <= 10000000; n *= 10) {
    std::string s(n, 'A');
    auto graph = Graph({s});
    impl::Depth<Graph::GraphCursor> depth;
    auto cursor = graph.begins()[0];
    EXPECT_EQ(depth.depth(cursor), n);
  }
}

