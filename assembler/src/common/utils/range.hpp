#pragma once

#include "utils/verify.hpp"

namespace omnigraph {

struct Range {
private:
    bool inside(size_t left, size_t right, size_t point) const {
        return left <= point && point <= right;
    }

public:
    //inclusive
    size_t start_pos;
    //exclusive
    size_t end_pos;

    size_t size() const {
        VERIFY(end_pos >= start_pos);
        return end_pos - start_pos;
    }

    void shift(int shift) {
        VERIFY(shift > 0 || size_t(-shift) <= start_pos);
        start_pos += shift;
        end_pos += shift;
    }

    Range(): start_pos(0), end_pos(0) {
        VERIFY(end_pos >= start_pos);
    }

    Range(size_t start_pos, size_t end_pos)
            : start_pos(start_pos),
              end_pos(end_pos) {
        VERIFY(end_pos >= start_pos);
    }

    bool operator<(const Range &other) const {
      if (start_pos != other.start_pos)
        return start_pos < other.start_pos;
      return end_pos < other.end_pos;
    }

    bool contains(const Range& that) const {
        return start_pos <= that.start_pos && end_pos >= that.end_pos;
    }

    Range Merge(const Range &other) const {
        return Range(this->start_pos, other.end_pos);
    }

    Range Invert(size_t base_length) const {
        VERIFY(base_length >= end_pos);
        return Range(base_length - end_pos, base_length - start_pos);
    }

    Range& operator=(const Range& other) {
        start_pos = other.start_pos;
        end_pos = other.end_pos;
        return *this;
    }

    bool empty() const {
        return start_pos == end_pos;
    }

    bool Intersect(const Range &other) const {
        return inside(start_pos, end_pos, other.start_pos) || inside(start_pos, end_pos, other.end_pos) ||
                inside(other.start_pos, other.end_pos, start_pos);
    }

    bool IntersectLeftOf(const Range &other) const {
        return inside(start_pos, end_pos, other.start_pos) && inside(other.start_pos, other.end_pos, end_pos);
    }

    bool operator==(const Range &that) const {
        return start_pos == that.start_pos && end_pos == that.end_pos;
    }

    bool operator!=(const Range &that) const {
        return !(*this == that);
    }
};

inline std::ostream& operator<<(std::ostream& os, const Range& range) {
    os << "[" << (range.start_pos + 1) << " - " << range.end_pos << "]";
    return os;
}

}
