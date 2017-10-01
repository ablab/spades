//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "sequence/sequence.hpp"
#include "sequence/range.hpp"
//FIXME bad dependency
#include "io/reads/single_read.hpp"
#include <boost/iterator/iterator_facade.hpp>

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

    Path() : sequence_(),
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

struct MappingRange {
// on genome/contig/whatever
    Range initial_range;
//on edge
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

    MappingRange ShiftInitial(int shift) const {
        MappingRange result(*this);
        result.initial_range.shift(shift);
        return result;
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
    MappingRange operator = (const MappingRange & other) {
        initial_range = other.initial_range;
        mapped_range = other.mapped_range;
        return *this;
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

    MappingPath(const ElementId &edge,
                const MappingRange &range_mapping)
            : edges_({ edge }),
              range_mappings_({ range_mapping }) {}
    
    MappingPath(const std::vector<ElementId>& edges,
                const std::vector<MappingRange> range_mappings)
            : edges_(edges),
              range_mappings_(range_mappings) {}

    MappingPath(const std::vector<pair<ElementId, MappingRange>>& edge_mappings) {
        edges_.reserve(edge_mappings.size());
        range_mappings_.reserve(edge_mappings.size());
        for (const auto &em : edge_mappings) {
            edges_.push_back(em.first);
            range_mappings_.push_back(em.second);
        }
    }

    class mapping_path_iter : public boost::iterator_facade<mapping_path_iter,
            std::pair<const ElementId, const MappingRange>,
            boost::random_access_traversal_tag,
            std::pair<const ElementId, const MappingRange>> {
        friend class boost::iterator_core_access;

        const MappingPath &mapping_path_;
        size_t pos_;

        std::pair<const ElementId, const MappingRange> dereference() const {
            return mapping_path_[pos_];
        };

        bool equal(const mapping_path_iter &that) const {
            return &mapping_path_ == &that.mapping_path_ && pos_ == that.pos_;
        }

        ptrdiff_t distance_to(const mapping_path_iter &that) const {
            return that.pos_ - pos_;
        }

        void advance(ptrdiff_t n) {
            pos_ += n;
        }

        void increment() {
            advance(1);
        }

        void decrement() {
            advance(-1);
        }

    public:
        mapping_path_iter(const MappingPath &mapping_path, size_t pos) :
                mapping_path_(mapping_path),
                pos_(pos) {}
    };

    mapping_path_iter begin() const {
        return mapping_path_iter(*this, 0);
    }

    mapping_path_iter end() const {
        return mapping_path_iter(*this, size());
    };

    size_t size() const { return edges_.size(); }

    size_t empty() const { return edges_.empty(); }

    ElementId edge_at(size_t idx) const {
       return edges_[idx];
    };

    MappingRange mapping_at(size_t idx) const {
        return range_mappings_[idx];
    };

    std::pair<const ElementId, const MappingRange> operator[](size_t idx) const {
        VERIFY(idx < size());
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

    void join(const MappingPath<ElementId>& that, int pos_shift = 0) {
        for (size_t i = 0; i < that.size(); ++i) {
            edges_.push_back(that.edges_[i]);
            range_mappings_.push_back(that.range_mappings_[i].ShiftInitial(pos_shift));
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

template <typename ElementId>
inline std::ostream& operator<<(std::ostream& os, const MappingPath<ElementId>& mp) {
    os << "MappingPath ( ";
    for(size_t i = 0; i < mp.size(); i++) {
        os << mp[i] << " ";
    }
    os << " )";
    return os;
}

inline boost::optional<Sequence> Subseq(const io::SingleRead& read, size_t start, size_t end) {
    //DEBUG("Requesting subseq of read length " << read.size() << " from " << start << " to " << end);
    VERIFY(end > start);
    //if (end == start) {
    //    DEBUG("Returning empty sequence");
    //    return boost::make_optional(Sequence());
    //}
    auto subread = read.Substr(start, end);
    if (subread.IsValid()) {
        //DEBUG("Gap seq valid. Length " << subread.size());
        return boost::make_optional(subread.sequence());
    } else {
        //DEBUG("Gap seq invalid. Length " << subread.size());
        //DEBUG("sequence: " << subread.GetSequenceString());
        return boost::none;
    }
}

inline boost::optional<Sequence> Subseq(const io::SingleReadSeq& read, size_t start, size_t end) {
    return boost::make_optional(read.sequence().Subseq(start, end));
}

inline boost::optional<Sequence> Subseq(const Sequence& s, size_t start, size_t end) {
    return boost::make_optional(s.Subseq(start, end));
}

//FIXME move somewhere else
template<class Graph>
class GapDescription {
    typedef typename Graph::EdgeId EdgeId;

    //Edges on the left and on the right of the gap
    EdgeId left_;
    EdgeId right_;

    //Estimated nucleotide gap/overlap between the edges !after trimming! (see further).
    // Negative values indicate the overlap between edges.
    // Should be non-negative for proper final joinings.
    int estimated_dist_;

    //Number of nucleotides to trim from the (end of the left)/(beginning of the right) edge
    size_t left_trim_;
    size_t right_trim_;

    //Optional "filling" sequence, giving "additional" nucleotides which
    // should be added while closing the gap.
    // Length guaranteed to be equal to estimated_gap (if present).
    boost::optional<Sequence> filling_seq_;

    GapDescription(EdgeId left, EdgeId right,
                   int estimated_dist,
                   size_t left_trim, size_t right_trim,
                   boost::optional<Sequence> filling_seq) :
            left_(left),
            right_(right),
            estimated_dist_(estimated_dist),
            left_trim_(left_trim),
            right_trim_(right_trim),
            filling_seq_(filling_seq) {
    }

    auto AsTuple() const ->
    decltype(std::make_tuple(left_, right_, left_trim_, right_trim_, estimated_dist_, filling_seq_)) {
        return std::make_tuple(left_, right_, left_trim_, right_trim_, estimated_dist_, filling_seq_);
    }

public:
    static const int INVALID_GAP = std::numeric_limits<int>::min();
    //static const GapDescription INVALID_GAP_DESC; //= GapDescription();

    template<class ReadT>
    static GapDescription CreateFixOverlap(const Graph &g, const ReadT& read, 
                                           size_t seq_start, size_t seq_end,
                                           EdgeId left, size_t left_offset,
                                           EdgeId right, size_t right_offset) {
        VERIFY(left_offset > 0 && right_offset >= 0 && 
                left_offset <= g.length(left) && right_offset < g.length(right));

        DEBUG("Creating gap description");

        //trying to shift on the left edge
        if (seq_start >= seq_end) {
            //+1 is a trick to avoid empty gap sequences
            size_t overlap = seq_start - seq_end + 1;
            DEBUG("Overlap of size " << overlap << " detected. Fixing.");
            size_t left_shift = std::min(overlap, left_offset - 1);
            VERIFY(seq_start >= left_shift);
            seq_start -= left_shift;
            left_offset -= left_shift;
        }

        //trying to shift on the right edge
        if (seq_start >= seq_end) {
            //+1 is a trick to avoid empty gap sequences
            size_t overlap = seq_start - seq_end + 1;
            DEBUG("Overlap of size " << overlap << " remained. Fixing.");
            size_t right_shift = std::min(overlap, g.length(right) - right_offset - 1);
            VERIFY(seq_end + right_shift <= read.size());
            seq_end += right_shift;
            right_offset += right_shift;
        }

        if (seq_start < seq_end) {
            auto gap_seq = Subseq(read, seq_start, seq_end);
            if (gap_seq) {
                DEBUG("Gap info successfully created");
                VERIFY(left_offset > 0 && right_offset >= 0 && 
                        left_offset <= g.length(left) && right_offset < g.length(right));
                return GapDescription(left, right,
                                      *gap_seq,
                                      g.length(left) - left_offset,
                                      right_offset);
            } else {
                DEBUG("Something wrong with read subsequence");
            }
        } else {
            size_t overlap = seq_start - seq_end + 1;
            //FIXME remove warn
            WARN("Failed to fix overlap of size " << overlap);
        }
        return GapDescription();
    }

    GapDescription(EdgeId left, EdgeId right,
                   int estimated_dist,
                   size_t left_trim = 0, size_t right_trim = 0) :
            GapDescription(left, right,
                           estimated_dist,
                           left_trim, right_trim,
                           boost::none) {
    }

    GapDescription() : GapDescription(EdgeId(0), EdgeId(0), INVALID_GAP) {
    }

    GapDescription(EdgeId left, EdgeId right,
                   const Sequence &filling_seq,
                   size_t left_trim = 0, size_t right_trim = 0) :
            left_(left),
            right_(right),
            estimated_dist_(int(filling_seq.size())),
            left_trim_(left_trim),
            right_trim_(right_trim),
            filling_seq_(filling_seq) {
    }

    EdgeId left() const {
        return left_;
    }

    EdgeId right() const {
        return right_;
    }

    size_t left_trim() const {
        return left_trim_;
    }

    size_t right_trim() const {
        return right_trim_;
    }

    bool no_trim() const {
        return left_trim_ == 0 && right_trim() == 0;
    }

    int estimated_dist() const {
        return estimated_dist_;
    }

    bool has_filling() const {
        return filling_seq_;
    }

    Sequence filling_seq() const {
        return *filling_seq_;
    }

    void set_left(EdgeId e) {
        left_ = e;
    }

    void set_right(EdgeId e) {
        right_ = e;
    }

    void set_estimated_dist(int dist) {
        VERIFY_MSG(!filling_seq_, "Filling sequence specified");
        estimated_dist_ = dist;
    }

    void set_filling_seq(Sequence fill_seq) {
        estimated_dist_ = fill_seq.size();
        filling_seq_ = boost::make_optional(fill_seq);
    }

    GapDescription<Graph> conjugate(const Graph &g) const {
        GapDescription<Graph> res(g.conjugate(right_),
                                  g.conjugate(left_),
                                  estimated_dist_,
                                  right_trim_,
                                  left_trim_,
                                  filling_seq_ ? boost::make_optional(!*filling_seq_) : boost::none);
        return res;
    }

    string str(const Graph &g) const {
        stringstream s;
        s << "left: " << g.str(left_)
          << "; right: " << g.str(right_)
          << "; estimated distance : " << estimated_dist_
          << "; left trim: " << left_trim_
          << "; right trim: " << right_trim_
          << "; sequence " << (filling_seq_ ? filling_seq_->str() : "no_sequence") << endl;
        return s.str();
    }

    bool operator<(const GapDescription &rhs) const {
        return AsTuple() < rhs.AsTuple();
    }

    bool operator!=(const GapDescription rhs) const {
        return AsTuple() != rhs.AsTuple();
    }

    bool operator==(const GapDescription rhs) const {
        return !(*this != rhs);
    }
private:
    DECL_LOGGER("GapDescription");

};

}
