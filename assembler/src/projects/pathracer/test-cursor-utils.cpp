#include <gtest/gtest.h>
#include "graph.hpp"
#include "hmmpath.hpp"

#include "cursor_utils.hpp"

TEST(ChechIsolatedLoopDetection, cursor_utils_hpp) {
    size_t n = 1000, k = 21;
    std::string s(n, 'A');
    auto graph = DBGraph(k, {s});
    auto vcursors = vertex_cursors(graph.all(), nullptr);
    EXPECT_EQ(vcursors.size(), 1);
}
