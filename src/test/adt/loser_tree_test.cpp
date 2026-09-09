
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "adt/loser_tree.hpp"
#include "utils/stl_utils.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <algorithm>

template <typename LoserTree>
auto get(LoserTree &lt, size_t n = size_t(-1)) {
    std::vector<typename LoserTree::value_type> result;
    lt.multi_merge_unique(std::back_inserter(result), n);
    return result;
}

template <class T0, class... Ts>
auto make_vector(T0 &&first, Ts &&... args) {
    using first_type = std::decay_t<T0>;
    return std::vector<first_type>{std::forward<T0>(first), std::forward<Ts>(args)...};
}

TEST(LoserTree, empty_test) {
    std::vector<int> v1 = {};
    std::vector<int> v2 = {};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});
    EXPECT_TRUE(lt.empty());

    EXPECT_EQ(get(lt, 0), std::vector<int>());
    EXPECT_EQ(get(lt, 0), std::vector<int>());
    EXPECT_EQ(get(lt, 1), std::vector<int>());
    EXPECT_EQ(get(lt, 100), std::vector<int>());
    EXPECT_TRUE(lt.empty());
}

TEST(LoserTree, one_empty) {
    std::vector<int> v1 = {1, 2, 2, 2, 3, 5};
    std::vector<int> v2 = {};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});

    EXPECT_EQ(get(lt, 0), std::vector<int>());
    EXPECT_EQ(get(lt, 0), std::vector<int>());
    EXPECT_EQ(get(lt, 1), std::vector<int>({1}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({2, 3}));
    EXPECT_EQ(get(lt, 100), std::vector<int>({5}));
    EXPECT_TRUE(lt.empty());
}

TEST(LoserTree, common1) {
    std::vector<int> v1 = {1, 1, 1, 1, 5, 5, 5};
    std::vector<int> v2 = {1, 2, 5};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});

    EXPECT_EQ(get(lt, 0), std::vector<int>());
    EXPECT_EQ(get(lt, 0), std::vector<int>());
    EXPECT_EQ(get(lt, 1), std::vector<int>({1}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({2, 5}));
    EXPECT_EQ(get(lt, 100), std::vector<int>({}));
    EXPECT_TRUE(lt.empty());
}

TEST(LoserTree, repeated) {
    std::vector<int> v1 = {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5};
    std::vector<int> v2 = {1, 2, 5};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});

    EXPECT_EQ(get(lt, 0), std::vector<int>());
    EXPECT_EQ(get(lt, 0), std::vector<int>());
    EXPECT_EQ(get(lt, 1), std::vector<int>({1}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({2, 3}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({5}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({}));
    EXPECT_TRUE(lt.empty());
}

TEST(LoserTree, get_all) {
    std::vector<int> v1 = {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5};
    std::vector<int> v2 = {1, 2, 5};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});

    EXPECT_EQ(get(lt), std::vector<int>({1, 2, 3, 5}));
    EXPECT_TRUE(lt.empty());

    EXPECT_EQ(get(lt), std::vector<int>({}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({}));
    EXPECT_EQ(get(lt), std::vector<int>({}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({}));
    EXPECT_TRUE(lt.empty());
}

TEST(LoserTree, get_all2) {
    std::vector<int> v1 = {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5};
    std::vector<int> v2 = {2, 5};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});

    EXPECT_EQ(get(lt, 4), std::vector<int>({1, 2, 3, 5}));
    EXPECT_TRUE(lt.empty());

    EXPECT_EQ(get(lt), std::vector<int>({}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({}));
    EXPECT_EQ(get(lt), std::vector<int>({}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({}));
    EXPECT_TRUE(lt.empty());
}

TEST(LoserTree, threeway) {
    std::vector<int> v1 = {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5, 30};
    std::vector<int> v2 = {2, 5, 29, 29, 29};
    std::vector<int> v3 = {-1, 4, 28, 28, 29, 30};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()),
                                    adt::make_range(v2.cbegin(), v2.cend()),
                                    adt::make_range(v3.cbegin(), v3.cend())});

    EXPECT_EQ(get(lt, 4), std::vector<int>({-1, 1, 2, 3}));
    EXPECT_EQ(get(lt, 4), std::vector<int>({4, 5, 28, 29}));
    EXPECT_EQ(get(lt, 4), std::vector<int>({30}));
    EXPECT_TRUE(lt.empty());
    EXPECT_EQ(get(lt), std::vector<int>({}));
    EXPECT_EQ(get(lt, 2), std::vector<int>({}));
}

TEST(LoserTree, threeway2) {
    std::vector<int> v1 = {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5, 30};
    std::vector<int> v2 = {2, 5, 29, 29, 29};
    std::vector<int> v3 = {-1, 4, 28, 28, 29, 30};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()),
                                    adt::make_range(v2.cbegin(), v2.cend()),
                                    adt::make_range(v3.cbegin(), v3.cend())});

    EXPECT_EQ(get(lt, 1), std::vector<int>({-1}));
    EXPECT_EQ(get(lt, 1), std::vector<int>({1}));
    EXPECT_EQ(get(lt, 1), std::vector<int>({2}));
    EXPECT_EQ(get(lt, 1), std::vector<int>({3}));
    EXPECT_EQ(get(lt, 1), std::vector<int>({4}));
    EXPECT_EQ(get(lt, 1), std::vector<int>({5}));
    EXPECT_EQ(get(lt, 1), std::vector<int>({28}));
    EXPECT_EQ(get(lt, 1), std::vector<int>({29}));
    EXPECT_EQ(get(lt, 1), std::vector<int>({30}));
    EXPECT_EQ(get(lt, 1), std::vector<int>({}));
    EXPECT_TRUE(lt.empty());
}
