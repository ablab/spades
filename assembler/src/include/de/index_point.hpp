//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <btree/btree_set.h>

namespace omnigraph {

namespace de {

// Define several storage-only POD types which can be
// implicitly converted to / from double.

class DEDistance {
public:
    DEDistance() = default;
    DEDistance(int d)
            : d_((float)d) {}
    DEDistance(double d)
            : d_((float)d) {}
    DEDistance(size_t d)
            : d_((float)d) {}
    operator float() const { return d_; }
    DEDistance operator+= (double d) {
        d_ += (float)d;
        return *this;
    }
    DEDistance operator*= (double d) {
        d_ *= (float)d;
        return *this;
    }
private:
    float d_;
};

class DEWeight {
public:
    DEWeight() = default;
    DEWeight(double d)
            : d_((float)d) {}
    operator float() const { return d_; }
    DEWeight operator+= (double d) {
        d_ += (float)d;
        return *this;
    }
    DEWeight operator*= (double d) {
        d_ *= (float)d;
        return *this;
    }
private:
    float d_;
};

struct RawPoint {
    DEDistance d;
    mutable DEWeight weight;

    RawPoint()
            : d(0.0), weight(0.0) {}

    RawPoint(DEDistance distance, DEWeight weight)
            : d(distance), weight(weight) {}

    RawPoint(DEDistance distance, DEWeight weight, DEDistance)
            : d(distance), weight(weight) {}


    const RawPoint operator+=(const RawPoint &rhs) const {
        weight += rhs.weight;
        return *this;
    }

    std::string str() const {
        stringstream ss;
        ss << "Point: " << " distance = " << this->d
           << ", weight = " << this->weight;
        return ss.str();
    }

    RawPoint& operator=(const RawPoint& rhs) {
        using namespace math;
        update_value_if_needed<DEDistance>(d, rhs.d);
        update_value_if_needed<DEWeight>(weight, rhs.weight);
        return *this;
    }

    bool operator<(const RawPoint& rhs) const {
        return math::ls(this->d, rhs.d);
    }

    bool operator==(const RawPoint& rhs) const {
        return math::eq(this->d, rhs.d);
    }

    bool operator!=(const RawPoint& rhs) const {
        return !(operator==(rhs));
    }

    RawPoint operator-() const {
        return RawPoint(-d, weight);
    }

    RawPoint operator+(const RawPoint &rhs) const {
        return RawPoint(d, rhs.weight + weight);
    }

    DEWeight variation() const {
        return 0;
    }

    RawPoint Conjugate(size_t l1, size_t l2) const
    {
        return RawPoint(d + DEDistance(l2) - DEDistance(l1), weight);
    }
};

struct Point : public RawPoint {
    DEDistance var;
    Point()
            : var(0.0) {}

    Point(DEDistance distance, DEWeight weight, DEDistance variance)
            : RawPoint(distance, weight), var(variance) {}

    Point(const Point &rhs)
            : RawPoint(rhs), var(rhs.var) {}

    Point(const RawPoint &rhs)
            : RawPoint(rhs), var(0.0) {}

    Point& operator=(const Point& rhs) {
        using namespace math;
        update_value_if_needed<DEDistance>(d, rhs.d);
        update_value_if_needed<DEWeight>(weight, rhs.weight);
        update_value_if_needed<DEDistance>(var, rhs.var);
        return *this;
    }

    bool operator<(const Point& rhs) const {
        return math::ls(this->d, rhs.d);
    }

    bool operator==(const Point& rhs) const {
        return math::eq(this->d, rhs.d);
    }

    bool operator!=(const Point& rhs) const {
        return !(operator==(rhs));
    }

    Point operator-() const {
        return Point(-d, weight, var);
    }

    Point operator+(const Point &rhs) const {
        auto weight_rhs = rhs.weight;
        // counting new bounds in the case, when we are merging pair infos with var != 0
        auto left_bound = std::min(d - var, rhs.d - rhs.var);
        auto right_bound = std::max(d + var, rhs.d + rhs.var);
        auto new_dist = (left_bound + right_bound) * 0.5f;
        auto new_weight = weight + weight_rhs;
        auto new_variance = (right_bound - left_bound) * 0.5f;

        return Point(new_dist, new_weight, new_variance);
    }

    DEDistance variation() const {
        return var;
    }

    Point Conjugate(size_t l1, size_t l2) const
    {
        return Point(d + DEDistance(l2) - DEDistance(l1), weight, var);
    }
};

inline int rounded_d(const RawPoint& p) {
    return math::round_to_zero(p.d);
}

inline std::ostream& operator<<(std::ostream& os, const Point &point) {
    return os << point.str();
}

inline std::ostream& operator<<(std::ostream& os, const RawPoint &point) {
    return os << point.str();
}

template<class Point>
class Histogram {
    typedef Histogram<Point> self_type;
    typedef typename std::less<Point> key_compare;
    typedef typename std::allocator<Point> allocator_type;
    typedef typename btree::btree<btree::btree_set_params<Point, key_compare, allocator_type, 1024> > Tree;

  public:
    typedef typename Tree::params_type params_type;
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

 public:
    // Default constructor.
    Histogram()
            : tree_(key_compare(), allocator_type()) {
    }

    // Copy constructor.
    Histogram(const self_type &x)
            : tree_(x.tree_) {
    }

    template <class InputIterator>
    Histogram(InputIterator b, InputIterator e)
            : tree_(key_compare(), allocator_type()) {
        insert(b, e);
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
    void dump(std::ostream &os) const { tree_.dump(os); }
    void verify() const { tree_.verify(); }

    // Size routines.
    size_type size() const { return tree_.size(); }
    size_type max_size() const { return tree_.max_size(); }
    bool empty() const { return tree_.empty(); }
    size_type height() const { return tree_.height(); }
    size_type internal_nodes() const { return tree_.internal_nodes(); }
    size_type leaf_nodes() const { return tree_.leaf_nodes(); }
    size_type nodes() const { return tree_.nodes(); }
    size_type bytes_used() const { return tree_.bytes_used(); }
    static double average_bytes_per_value() { return Tree::average_bytes_per_value(); }
    double fullness() const { return tree_.fullness(); }
    double overhead() const { return tree_.overhead(); }

    // Lookup routines.
    iterator find(const key_type &key) { return tree_.find_unique(key); }
    const_iterator find(const key_type &key) const { return tree_.find_unique(key); }
    size_type count(const key_type &key) const { return tree_.count_unique(key); }

    // Insertion routines.
    std::pair<iterator,bool> insert(const value_type &x) { return tree_.insert_unique(x); }
    iterator insert(iterator position, const value_type &x) { return tree_.insert_unique(position, x); }
    template <typename InputIterator>
    void insert(InputIterator b, InputIterator e) { tree_.insert_unique(b, e); }

    // Deletion routines.
    int erase(const key_type &key) { return tree_.erase_unique(key); }
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
    template<class OtherHist, class U = Point>
    typename std::enable_if<can_merge<U>::value, size_t>::type
    merge(const OtherHist &other) {
        size_t added = 0;
        for (const auto& new_point : other) {
            // First, try to insert a point
            const auto& result = insert(new_point);
            if (!result.second) {
                // We already having something there. Try to merge stuff in.
                *result.first += new_point;
            } else
                added += 1;
        }

        return added;
    }

    // Otherwise this overload is used, which removes the point from set,
    // updates it and re-inserts back.
    template<class OtherHist, class U = Point>
    typename std::enable_if<!can_merge<U>::value, size_t>::type
    merge(const OtherHist &other) {
        size_t added = 0;
        for (const auto& new_point : other) {
            // First, try to insert a point
            auto result = insert(new_point);
            if (!result.second) {
                Point updated = *result.first + new_point;
                auto after_removed = erase(result.first);
                insert(after_removed, updated);
            } else
                added += 1;
        }

        return added;
    }
};

template <typename T>
inline std::ostream& operator<<(std::ostream &os, const Histogram<T> &b) {
  b.dump(os);
  return os;
}

typedef Histogram<RawPoint> RawHistogram;
typedef Histogram<Point> HistogramWithWeight;

inline bool ClustersIntersect(Point p1, Point p2) {
    return math::le(p1.d, p2.d + p1.var + p2.var) &&
           math::le(p2.d, p1.d + p1.var + p2.var);
}

// tuple of a pair of edges @first, @second, and a @point
template<typename EdgeId>
struct PairInfo {
    EdgeId first;
    EdgeId second;
    Point point;

    PairInfo()
            : first(), second(), point() {}


    PairInfo(const PairInfo& pair_info)
            : first(pair_info.first), second(pair_info.second), point(pair_info.point) {}

    PairInfo(EdgeId first, EdgeId second, DEDistance d, DEWeight weight, DEDistance var)
            : first(first), second(second), point(d, weight, var) {}

    PairInfo(EdgeId first, EdgeId second, Point point)
            : first(first), second(second), point(point) {}

    // Two paired infos are considered equal
    // if they coincide in all parameters except for weight and variance.
    bool operator==(const PairInfo& rhs) const {
        const PairInfo &lhs = *this;
        return lhs.first == rhs.first && lhs.second == rhs.second && lhs.point == rhs.point;
    }

    bool operator!=(const PairInfo& rhs) const {
        return !(*this == rhs);
    }

    bool operator<(const PairInfo<EdgeId>& rhs) const {
        const PairInfo<EdgeId>& lhs = *this;
        return lhs.first == rhs.first ?
               (lhs.second == rhs.second ? lhs.point < rhs.point : lhs.second < rhs.second)
               : lhs.first < rhs.first;
    }

    double d() const      { return point.d;      }
    double weight() const { return point.weight; }
    double var() const    { return point.var;    }
};

template<typename EdgeId>
ostream& operator<<(ostream& os, const PairInfo<EdgeId>& info) {
    return os << "PairInfo: first = " << info.first << ", second = " << info.second
           << "Point : " << info.point;
}

template<typename EdgeId>
const PairInfo<EdgeId> MinPairInfo(EdgeId id) {
    return PairInfo<EdgeId>(id, EdgeId(typename EdgeId::pointer_type(1)),
                            -10000000000, 0., 0.);
}

template<typename EdgeId>
const PairInfo<EdgeId> MaxPairInfo(EdgeId id) {
    return PairInfo<EdgeId>(id, EdgeId(typename EdgeId::pointer_type(-1)),
                            10000000000, 0., 0.);
}

template<typename EdgeId>
const PairInfo<EdgeId> MinPairInfo(EdgeId e1, EdgeId e2) {
    PairInfo<EdgeId> info = MinPairInfo(e1);
    info.second = e2;
    return info;
}

template<typename EdgeId>
const PairInfo<EdgeId> MaxPairInfo(EdgeId e1, EdgeId e2) {
    PairInfo<EdgeId> info = MaxPairInfo(e1);
    info.second = e2;
    return info;
}

/**
 * Method returns approximate distance between occurrences of edges in genome rounded to the nearest
 * integer. In case of a tie closest to 0 value is chosen thus one can assume that distance
 * is rounded the same way as opposite one.
 * todo check that written here is true
 */
template<typename EdgeId>
inline int rounded_d(PairInfo<EdgeId> const& pi) {
    return math::round_to_zero(pi.d());
}

template<typename EdgeId>
inline PairInfo<EdgeId> BackwardInfo(const PairInfo<EdgeId>& pi) {
    return PairInfo<EdgeId>(pi.second, pi.first, -pi.point);
}

template<typename EdgeId>
inline bool IsSymmetric(PairInfo<EdgeId> const& pi) {
    return pi.first == pi.second && math::eq(pi.d(), 0.);
}

}

}
