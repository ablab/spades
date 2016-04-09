//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>

namespace omnigraph {

namespace de {

/**
 * @brief Type for measuring distance between edges.
 *        ---edge1--->      --edge2->
 *        |----distance-----|
 *        Can be negative if edges are overlapped.
 *        For paired reads, distance roughly equals to insert size.
 */
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

/**
 * @brief Scalar type for measuring the weight of a point.
 *        Should not be negative.
 */
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

/**
 * @brief Scalar type for measuring the variance of point distance in clustered index.
 *        |---- max distance ----|
 *             |min distance|
 *        |var-|
 *        Should not be negative.
 */
typedef float DEVariance;

/**
 * @brief A proof-of-concept bicyclish wrapper for small integer types with saturated arithmetics.
 *        All operations are in the interval of [-OFFSET .. MAX_VAL - OFFSET]
 *        where MAX_VAL is 2^sizeof(T) - 1.
 */
template<typename T, T OFFSET = 0>
class DESat {
public:
    DESat(): _val(OFFSET) {}
    DESat(int val): _val(shrink(val)) {}
    DESat(float val): _val(shrink((int)val)) {}

    operator float() const {
        return (float)_val - OFFSET;
    }

    DESat& operator+=(DESat rhs) {
        _val = sadd(_val, rhs._val);
        if (_val > OFFSET) //subtract the second offset came from the rhs
            _val -= OFFSET;
        else
            _val = 0;
    }

    static DESat Max() {
        return DESat(std::numeric_limits<T>::max());
    }

private:
    DESat(T val): _val(val) {}

    T _val;
    static const T MAX_VAL = std::numeric_limits<T>::max();

    template<typename F>
    static T shrink(F d) {
        //Saturate to the allowed interval
        //TODO: optimize
        F t = d + OFFSET;
        t = std::max(F(0), t);
        t = std::min(F(MAX_VAL), t);
        return (T)t;
    }

    static T sadd(T a, T b) {
        //TODO: check
        return (a > MAX_VAL - b) ? MAX_VAL : a + b;
    }
};

/**
 * @brief Type for measuring a small gap between edges.
 *        Gap equals to distance minus the length of the first edge, and can be slightly negative in case of overlaps.
 *        ---edge1--->      --edge2->
 *                   |-gap--|
 */
//typedef DESat<uint16_t, 512> DEGap;
typedef float DEGap;

/**
 * @brief Type for weighting points in a paired info index.
 */
//typedef DESat<uint16_t> DECropWeight;
typedef float DECropWeight;

/**
 * @brief Raw point of unclustered index. Parameterized by distance and weight types.
 */
template<typename D, typename W>
struct __attribute((aligned(sizeof(D) + sizeof(W)))) RawPointT {
    typedef RawPointT<D, W> Self;
    D d;
    mutable W weight;

    RawPointT()
        : RawPointT(0.0, 0.0) {}

    RawPointT(D distance, W weight)
        : d(distance), weight(weight) {}

    Self& operator+=(const Self &rhs) {
        weight += rhs.weight;
        return *this;
    }

    std::string str() const {
        std::ostringstream ss;
        ss << "Point: " << " distance = " << this->d
           << ", weight = " << this->weight;
        return ss.str();
    }

    bool operator<(const Self& rhs) const {
        return math::ls(this->d, rhs.d);
    }

    bool operator==(const Self& rhs) const {
        return math::eq(this->d, rhs.d);
    }

    bool operator!=(const Self& rhs) const {
        return !(operator==(rhs));
    }

    Self operator-() const {
        return Self(-d, weight);
    }

    Self operator+(const Self &rhs) const {
        return Self(d, rhs.weight + this->weight);
    }

    DEVariance variance() { return 0; } //TODO: remove
};

typedef RawPointT<DEDistance, DEWeight> RawPoint;
typedef RawPointT<DEGap, DECropWeight> RawGapPoint;

inline int rounded_d(const RawPoint& p) {
    return math::round_to_zero(p.d);
}

inline std::ostream& operator<<(std::ostream& os, const RawPoint &point) {
    return os << point.str();
}

/**
 * @brief Clustered index point. Parameterized by distance and weight types, also has variance.
 */
template<typename D, typename W>
struct PointT : public RawPointT<D, W> {
    typedef PointT<D, W> Self;
    DEVariance var;
    PointT()
        : PointT(0.0, 0.0, 0.0) {}

    PointT(D distance, W weight, DEVariance variance)
        : RawPointT<D, W>(distance, weight), var(variance) {}

    PointT(const RawPointT<D, W> &rhs)
        : RawPointT<D, W>(rhs), var(0.0) {}

    bool operator<(const Self& rhs) const {
        return math::ls(this->d, rhs.d);
    }

    bool operator==(const Self& rhs) const {
        return math::eq(this->d, rhs.d);
    }

    bool operator!=(const Self& rhs) const {
        return !(operator==(rhs));
    }

    Self operator+(const Self &rhs) const {
        // counting new bounds in the case, when we are merging pair infos with var != 0
        auto left_bound = std::min(this->d - var, rhs.d - rhs.var);
        auto right_bound = std::max(this->d + var, rhs.d + rhs.var);
        auto new_dist = DEDistance((left_bound + right_bound) * 0.5f);
        auto new_weight = this->weight + rhs.weight; //TODO: crop
        auto new_variance = (right_bound - left_bound) * 0.5f;

        return Self(new_dist, new_weight, new_variance);
    }

    bool lt(const Self &rhs) const {
        return math::ls(this->weight, rhs.weight);
    }

    DEVariance variance() { return this->var; } //TODO: remove
};

typedef PointT<DEDistance, DEWeight> Point;
typedef PointT<DEGap, DECropWeight> GapPoint;

inline std::ostream& operator<<(std::ostream& os, const Point &point) {
    return os << point.str();
}

/**
 * @brief Policy-like type which provides associated point types for unclustered index
 *        and static methods for converting between them.
 */
struct RawPointTraits {
    typedef RawGapPoint Gapped;
    typedef RawPoint Expanded;

    static Gapped Shrink(Expanded p, DEDistance edge) {
        DEGap gap = DEGap(p.d - edge);
        return RawGapPoint(gap, p.weight);
    }

    static Expanded Expand(Gapped p, DEDistance edge) {
        RawPoint res(p.d, p.weight);
        res.d += edge;
        return res;
    }
};

/**
 * @brief Policy-like type which provides associated point types for clustered index
 *        and static methods for converting between them.
 */
struct PointTraits {
    typedef GapPoint Gapped;
    typedef Point Expanded;

    static Gapped Shrink(const Expanded &p, DEDistance edge) {
        DEGap gap = DEGap(p.d - edge);
        return GapPoint(gap, p.weight, p.var);
    }

    static Expanded Expand(Gapped p, DEDistance edge) {
        Point res(p.d, p.weight, p.var);
        res.d += edge;
        return res;
    }
};

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

}

}

namespace std {
template<>
class numeric_limits<omnigraph::de::DEDistance> : public numeric_limits<float> {};
template<>
class numeric_limits<omnigraph::de::DEWeight> : public numeric_limits<float> {};
}
