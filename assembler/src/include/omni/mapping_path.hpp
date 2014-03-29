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
    	return start_pos == that.start_pos || end_pos == that.end_pos;
    }

    bool operator!=(const Range &that) const {
    	return !(*this == that);
    }
};

inline std::ostream& operator<<(std::ostream& os, const Range& range) {
    os << "[" << (range.start_pos + 1) << " - " << range.end_pos << "]";
    return os;
}

struct MappingRange {
    Range initial_range;
    Range mapped_range;

    MappingRange() {
    }

    MappingRange(Range initial_range, Range mapped_range)
            : initial_range(initial_range), mapped_range(mapped_range) {}

    MappingRange(size_t i_start, size_t i_end, size_t m_start, size_t m_end)
            : initial_range(i_start, i_end), mapped_range(m_start, m_end) {}

    MappingRange Merge(const MappingRange &other) const {
    	return MappingRange(initial_range.Merge(other.initial_range), mapped_range.Merge(other.mapped_range));
    }

    MappingRange Shift(int shift) const {
    	VERIFY(initial_range.end_pos >= initial_range.start_pos);
    	if(empty())
    		return MappingRange();
    	MappingRange result(*this);
    	if(int(result.mapped_range.end_pos) <= -shift)
    		return MappingRange();
    	result.mapped_range.end_pos += shift;
    	if(int(result.mapped_range.start_pos) <= -shift) {
    		result.initial_range.start_pos -= result.mapped_range.start_pos + shift;
    		if(result.initial_range.start_pos >= result.initial_range.end_pos)
    			result.initial_range.start_pos = result.initial_range.end_pos - 1;
    		result.mapped_range.start_pos = 0;
    	} else {
    		result.mapped_range.start_pos += shift;
    	}
    	return result;
    }

    MappingRange Fit(size_t length) const {
    	VERIFY(initial_range.end_pos >= initial_range.start_pos);
    	if(empty())
    		return MappingRange();
    	MappingRange result(*this);
    	if(result.mapped_range.start_pos >= length)
    		return MappingRange();
    	if(result.mapped_range.end_pos >= length) {
    		if(result.initial_range.end_pos + length < result.mapped_range.end_pos)
        		return MappingRange();
    		result.initial_range.end_pos -= result.mapped_range.end_pos - length;
    		result.mapped_range.end_pos = length;
    	}
    	return result;
    }

    bool empty() const {
    	return initial_range.empty() || mapped_range.empty();
    }

    bool operator<(const MappingRange &other) const {
    	if(this->initial_range != other.initial_range)
    		return this->initial_range < other.initial_range;
    	return this->mapped_range < other.mapped_range;
    }

    bool Intersect(const MappingRange &other) {
    	return initial_range.Intersect(other.initial_range) && mapped_range.Intersect(other.mapped_range);
    }

    bool IntersectLeftOf(const MappingRange &other) const {
    	return initial_range.IntersectLeftOf(other.initial_range) && mapped_range.IntersectLeftOf(other.mapped_range);
    }

    bool StrictlyContinuesWith(const MappingRange &other, size_t max_gap, size_t gap_diff = 0) const {
        return this->initial_range.end_pos <= other.initial_range.start_pos 
                && this->mapped_range.end_pos <= other.mapped_range.start_pos 
                && other.initial_range.start_pos - this->initial_range.end_pos 
                    <= other.mapped_range.start_pos - this->mapped_range.end_pos + gap_diff
                && other.mapped_range.start_pos - this->mapped_range.end_pos 
                    <= other.initial_range.start_pos - this->initial_range.end_pos + gap_diff
                && other.initial_range.start_pos - this->initial_range.end_pos <= max_gap;
    }

    bool operator==(const MappingRange &that) const {
    	return initial_range == that.initial_range || mapped_range == that.mapped_range;
    }

    bool operator!=(const MappingRange &that) const {
    	return !(*this == that);
    }

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

    Path<ElementId> path() const {
        if (edges_.size() != 0)
            return Path<ElementId>(edges_,
                                   range_mappings_[0].mapped_range.start_pos,
                                   range_mappings_[range_mappings_.size() - 1].mapped_range.end_pos);
        else
            return Path<ElementId>();
    }

    const std::vector<ElementId>& simple_path() const {
        return edges_;
    }

    void join(const MappingPath<ElementId>& that) {
        for (size_t i = 0; i < that.size(); ++i) {
            edges_.push_back(that.edges_[i]);
            range_mappings_.push_back(that.range_mappings_[i]);
        }
    }

    void push_back(ElementId id, MappingRange range) {
        edges_.push_back(id);
        range_mappings_.push_back(range);
    }

 private:
    std::vector<ElementId> edges_;
    std::vector<MappingRange> range_mappings_;
};

}

#endif
