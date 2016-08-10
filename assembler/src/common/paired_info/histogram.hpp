//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <btree/btree_set.h>
#include "common/adt/flat_set.hpp"
#include "common/adt/small_pod_vector.hpp"
#include "index_point.hpp"

namespace omnigraph {

namespace de {

template<class Point>
class Histogram {
    typedef Histogram<Point> self_type;
    typedef typename std::less<Point> key_compare;
    typedef typename std::allocator<Point> allocator_type;
    typedef typename adt::flat_set<Point, key_compare, adt::SmallPODVector> Tree;

public:
    typedef typename Tree::key_type key_type;
    typedef typename Tree::value_type value_type;
    typedef typename Tree::pointer pointer;
    typedef typename Tree::const_pointer const_pointer;
    typedef typename Tree::reference reference;
    typedef typename Tree::const_reference const_reference;
    typedef typename Tree::size_type size_type;
    typedef typename Tree::difference_type difference_type;
    typedef typename Tree::iterator iterator;
    typedef typename Tree::const_iterator const_iterator;
    typedef typename Tree::reverse_iterator reverse_iterator;
    typedef typename Tree::const_reverse_iterator const_reverse_iterator;

    enum {
        kValueSize = sizeof(Point)
    };

public:
    // Default constructor.
    Histogram() = default;

    // Copy constructor.
    Histogram(const self_type &x)
            : tree_(x.tree_) {}

    template <class InputIterator>
    Histogram(InputIterator b, InputIterator e) {
        insert(b, e);
    }

    Histogram(std::initializer_list<Point> l) {
        insert(l.begin(), l.end());
    }

    // Iterator routines.
    iterator begin() { return tree_.begin(); }
    const_iterator begin() const { return tree_.begin(); }
    iterator end() { return tree_.end(); }
    const_iterator end() const { return tree_.end(); }
    reverse_iterator rbegin() { return tree_.rbegin(); }
    const_reverse_iterator rbegin() const { return tree_.rbegin(); }
    reverse_iterator rend() { return tree_.rend(); }
    const_reverse_iterator rend() const { return tree_.rend(); }

    // Lookup routines.
    iterator lower_bound(const key_type &key) { return tree_.lower_bound(key); }
    const_iterator lower_bound(const key_type &key) const { return tree_.lower_bound(key); }
    iterator upper_bound(const key_type &key) { return tree_.upper_bound(key); }
    const_iterator upper_bound(const key_type &key) const { return tree_.upper_bound(key); }
    std::pair<iterator,iterator> equal_range(const key_type &key) { return tree_.equal_range(key); }
    std::pair<const_iterator,const_iterator> equal_range(const key_type &key) const { return tree_.equal_range(key); }

    // Utility routines.
    void clear() { tree_.clear(); }
    void swap(self_type &x) { tree_.swap(x.tree_); }

    // Size routines.
    size_type size() const { return tree_.size(); }
    size_type max_size() const { return tree_.max_size(); }
    bool empty() const { return tree_.empty(); }
    size_type bytes_used() const { return tree_.bytes_used(); }

    // Lookup routines.
    iterator find(const key_type &key) { return tree_.find(key); }
    const_iterator find(const key_type &key) const { return tree_.find(key); }
    size_type count(const key_type &key) const { return tree_.count(key); }

    // Insertion routines.
    std::pair<iterator,bool> insert(const value_type &x) { return tree_.insert(x); }
    iterator insert(iterator position, const value_type &x) { return tree_.insert(position, x); }
    template <typename InputIterator>
    void insert(InputIterator b, InputIterator e) { tree_.insert(b, e); }

    // Deletion routines.
    size_type erase(const key_type &key) { return tree_.erase(key); }
    // Erase the specified iterator from the btree. The iterator must be valid
    // (i.e. not equal to end()).  Return an iterator pointing to the node after
    // the one that was erased (or end() if none exists).
    iterator erase(const iterator &iter) { return tree_.erase(iter); }
    void erase(const iterator &first, const iterator &last) { tree_.erase(first, last); }

    bool operator==(const self_type& x) const {
        if (size() != x.size())
            return false;

        for (const_iterator i = begin(), xi = x.begin(); i != end(); ++i, ++xi)
            if (*i != *xi)
                return false;

        return true;
    }

    bool operator!=(const self_type& other) const {
        return !operator==(other);
    }

protected:
    Tree tree_;

private:
    // This is template voodoo which creates function overload depending on
    // whether Point has const operator+= or not.
    template<class>
    struct true_helper : std::true_type {};
    template<class T = Point>
    static auto test_can_merge(int) -> true_helper<decltype(std::declval<const T>().operator+=(std::declval<const T>()))>;
    template<class>
    static auto test_can_merge(long) -> std::false_type;
    template<class T = Point>
    struct can_merge : decltype(test_can_merge<T>(0)) {};

public:
    // This function overload is enabled only when Point has const operator+= (e.g. RawPoint)
    // and therefore we can update it inplace.
    template<class U = Point>
    typename std::enable_if<can_merge<U>::value, size_t>::type
    merge_point(const U &new_point) {
        // First, try to insert a point
        const auto &result = insert(new_point);
        if (result.second)
            return 1;
        // We already having something there. Try to merge stuff in.
        *result.first += new_point;
        return 0;
    }

    // Otherwise this overload is used, which removes the point from set,
    // updates it and re-inserts back.
    template<class U = Point>
    typename std::enable_if<!can_merge<U>::value, size_t>::type
    merge_point(const U &new_point) {
        auto result = insert(new_point);
        if (result.second)
            return 1;
        Point updated = *result.first + new_point;
        auto after_removed = erase(result.first);
        insert(after_removed, updated);
        return 0;
    }

    template<class OtherHist>
    size_t merge(const OtherHist &other) {
        // If histogram is empty, we could simply insert everything
        if (size() == 0) {
            insert(other.begin(), other.end());
            return size();
        }

        size_t old_size = size();
        for (const auto &new_point : other)
            merge_point(new_point);
        return size() - old_size;
    }
};

template<typename T>
inline std::ostream &operator<<(std::ostream &os, const Histogram<T> &b) {
    os << "{";
    for (const auto& e : b)
        os << e << "; ";
    os << "}";
    return os;
}

typedef Histogram<RawGapPoint> RawGapHistogram;
typedef Histogram<GapPoint> GapHistogram;

typedef Histogram<RawPoint> RawHistogram;
typedef Histogram<Point> HistogramWithWeight;

}

}
