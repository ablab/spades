#include <gtest/gtest.h>
#include "graph.hpp"
#include "hmmpath.hpp"

TEST(InfinteDepthAtLeast, Depth) {
    size_t n = 1000, k = 21;
    std::string s0(n, 'A');
    auto s1 = std::string("ASDASDADADSDADADADADADADASDASDADADADASDA") + std::string(k, 'A');
    auto suf2 = std::string("dsfsfsdfsfsfsdfs");
    auto suf3 = std::string("21313123213131231231312312312*");
    auto s2 = std::string(k, 'A') + suf2;
    auto s3 = std::string(k, 'A') + suf3;
    auto graph = DBGraph(k, {s0, s1, s2, s3});
    depth_filter::DepthInt<DBGraph::GraphCursor> depth;
    for (size_t i = 0; i < s0.size(); ++i) {
        auto cursor = graph.get_pointer(0, i);
        EXPECT_EQ(depth.depth(cursor, nullptr), depth.INF);
    }
    for (size_t i = 0; i < s1.size(); ++i) {
        auto cursor = graph.get_pointer(1, i);
        EXPECT_EQ(depth.depth(cursor, nullptr), depth.INF);
    }

    for (size_t i = 0; i < suf2.size(); ++i) {
        auto cursor = graph.get_pointer(2, s2.size() - 1 - i);
        EXPECT_EQ(depth.depth(cursor, nullptr), i + 1);
    }

    for (size_t i = 0; i < suf3.size(); ++i) {
        auto cursor = graph.get_pointer(3, s3.size() - 1 - i);
        EXPECT_EQ(depth.depth(cursor, nullptr), i);
    }
}
