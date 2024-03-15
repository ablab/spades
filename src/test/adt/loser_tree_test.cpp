#define BOOST_TEST_MODULE loser_tree_test
#include <boost/test/included/unit_test.hpp>

#include "adt/loser_tree.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include "utils/stl_utils.hpp"

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

BOOST_AUTO_TEST_CASE(empty_test) {
    std::vector<int> v1 = {};
    std::vector<int> v2 = {};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});
    BOOST_CHECK(lt.empty());

    BOOST_CHECK_EQUAL(get(lt, 0), std::vector<int>());
    BOOST_CHECK_EQUAL(get(lt, 0), std::vector<int>());
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>());
    BOOST_CHECK_EQUAL(get(lt, 100), std::vector<int>());
    BOOST_CHECK(lt.empty());
}

BOOST_AUTO_TEST_CASE(one_empty) {
    std::vector<int> v1 = {1, 2, 2, 2, 3, 5};
    std::vector<int> v2 = {};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});

    BOOST_CHECK_EQUAL(get(lt, 0), std::vector<int>());
    BOOST_CHECK_EQUAL(get(lt, 0), std::vector<int>());
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({1}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({2, 3}));
    BOOST_CHECK_EQUAL(get(lt, 100), std::vector<int>({5}));
    BOOST_CHECK(lt.empty());
}

BOOST_AUTO_TEST_CASE(common1) {
    std::vector<int> v1 = {1, 1, 1, 1, 5, 5, 5};
    std::vector<int> v2 = {1, 2, 5};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});

    BOOST_CHECK_EQUAL(get(lt, 0), std::vector<int>());
    BOOST_CHECK_EQUAL(get(lt, 0), std::vector<int>());
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({1}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({2, 5}));
    BOOST_CHECK_EQUAL(get(lt, 100), std::vector<int>({}));
    BOOST_CHECK(lt.empty());
}

BOOST_AUTO_TEST_CASE(repeated) {
    std::vector<int> v1 = {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5};
    std::vector<int> v2 = {1, 2, 5};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});

    BOOST_CHECK_EQUAL(get(lt, 0), std::vector<int>());
    BOOST_CHECK_EQUAL(get(lt, 0), std::vector<int>());
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({1}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({2, 3}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({5}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({}));
    BOOST_CHECK(lt.empty());
}

BOOST_AUTO_TEST_CASE(get_all) {
    std::vector<int> v1 = {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5};
    std::vector<int> v2 = {1, 2, 5};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});

    BOOST_CHECK_EQUAL(get(lt), std::vector<int>({1, 2, 3, 5}));
    BOOST_CHECK(lt.empty());

    BOOST_CHECK_EQUAL(get(lt), std::vector<int>({}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({}));
    BOOST_CHECK_EQUAL(get(lt), std::vector<int>({}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({}));
    BOOST_CHECK(lt.empty());
}

BOOST_AUTO_TEST_CASE(get_all2) {
    std::vector<int> v1 = {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5};
    std::vector<int> v2 = {2, 5};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()), adt::make_range(v2.cbegin(), v2.cend())});

    BOOST_CHECK_EQUAL(get(lt, 4), std::vector<int>({1, 2, 3, 5}));
    BOOST_CHECK(lt.empty());

    BOOST_CHECK_EQUAL(get(lt), std::vector<int>({}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({}));
    BOOST_CHECK_EQUAL(get(lt), std::vector<int>({}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({}));
    BOOST_CHECK(lt.empty());
}

BOOST_AUTO_TEST_CASE(threeway) {
    std::vector<int> v1 = {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5, 30};
    std::vector<int> v2 = {2, 5, 29, 29, 29};
    std::vector<int> v3 = {-1, 4, 28, 28, 29, 30};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()),
                                    adt::make_range(v2.cbegin(), v2.cend()),
                                    adt::make_range(v3.cbegin(), v3.cend())});

    BOOST_CHECK_EQUAL(get(lt, 4), std::vector<int>({-1, 1, 2, 3}));
    BOOST_CHECK_EQUAL(get(lt, 4), std::vector<int>({4, 5, 28, 29}));
    BOOST_CHECK_EQUAL(get(lt, 4), std::vector<int>({30}));
    BOOST_CHECK(lt.empty());
    BOOST_CHECK_EQUAL(get(lt), std::vector<int>({}));
    BOOST_CHECK_EQUAL(get(lt, 2), std::vector<int>({}));
}

BOOST_AUTO_TEST_CASE(threeway2) {
    std::vector<int> v1 = {1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5, 30};
    std::vector<int> v2 = {2, 5, 29, 29, 29};
    std::vector<int> v3 = {-1, 4, 28, 28, 29, 30};
    auto lt = adt::make_loser_tree({adt::make_range(v1.cbegin(), v1.cend()),
                                    adt::make_range(v2.cbegin(), v2.cend()),
                                    adt::make_range(v3.cbegin(), v3.cend())});

    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({-1}));
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({1}));
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({2}));
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({3}));
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({4}));
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({5}));
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({28}));
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({29}));
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({30}));
    BOOST_CHECK_EQUAL(get(lt, 1), std::vector<int>({}));
    BOOST_CHECK(lt.empty());
}
