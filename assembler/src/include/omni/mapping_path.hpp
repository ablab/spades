#ifndef __OMNI_MAPPING_PATH_HPP__
#define __OMNI_MAPPING_PATH_HPP__

namespace omnigraph {

/**
 * This class is a representation of how certain sequence is mapped to genome. Needs further adjustment.
 */
template<typename ElementId>
class Path {
    std::vector<ElementId> sequence_;
    size_t start_pos_;
    size_t end_pos_;
 public:
    typedef typename vector<ElementId>::const_iterator iterator;

    Path(const vector<ElementId>& sequence, size_t start_pos, size_t end_pos)
            : sequence_(sequence), start_pos_(start_pos),  end_pos_( end_pos) {
    }

    Path()
            : sequence_(),
              start_pos_(-1ul),
              end_pos_(-1ul) {
    }

    size_t start_pos() const { return start_pos_; }
    size_t end_pos() const { return end_pos_; }

    size_t size() const { return sequence_.size(); }

    const std::vector<ElementId>& sequence() const { return sequence_; }
    ElementId operator[](size_t index) const { return sequence_[index]; }

    iterator begin() const { return sequence_.begin(); }
    iterator end() const { return sequence_.end(); }
};

struct Range {
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

    bool contains(const Range& that) {
        return start_pos <= that.start_pos && end_pos >= that.end_pos;
    }
};

inline std::ostream& operator<<(std::ostream& os, const Range& range) {
    os << "[" << range.start_pos << ", " << range.end_pos << "]";
    return os;
}

struct MappingRange {
    Range initial_range;
    Range mapped_range;

    MappingRange(Range initial_range, Range mapped_range)
            : initial_range(initial_range), mapped_range(mapped_range) {}
};

inline std::ostream& operator<<(std::ostream& os, const MappingRange& map_range) {
    os << map_range.initial_range << " --> " << map_range.mapped_range;
    return os;
}

template<typename ElementId>
class MappingPath {
 public:
    MappingPath() {}

    MappingPath(const std::vector<ElementId>& edges,
                const std::vector<MappingRange> range_mappings)
            : edges_(edges),
              range_mappings_(range_mappings) {}

    size_t size() const { return edges_.size(); }

    std::pair<const ElementId, const MappingRange> operator[](size_t idx) const {
        return std::make_pair(edges_[idx], range_mappings_[idx]);
    }

    std::pair<const ElementId, const MappingRange> front() const {
        return std::make_pair(edges_.front(), range_mappings_.front());
    }

    std::pair<const ElementId, const MappingRange> back() const {
        return std::make_pair(edges_.back(), range_mappings_.back());
    }

    size_t start_pos() const {
        return range_mappings_.front().mapped_range.start_pos;
    }

    size_t end_pos() const {
        return range_mappings_.back().mapped_range.end_pos;
    }

    Path<ElementId> simple_path() const {
        if (edges_.size() != 0)
            return Path<ElementId>(edges_,
                                   range_mappings_[0].mapped_range.start_pos,
                                   range_mappings_[range_mappings_.size() - 1].mapped_range.end_pos);
        else
            return Path<ElementId>();
    }

    void join(const MappingPath<ElementId>& that) {
        for (size_t i = 0; i < that.size(); ++i) {
            edges_.push_back(that.edges_[i]);
            range_mappings_.push_back(that.range_mappings_[i]);
        }
    }

 private:
    std::vector<ElementId> edges_;
    std::vector<MappingRange> range_mappings_;
};

}

#endif
