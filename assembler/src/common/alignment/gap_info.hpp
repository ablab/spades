#pragma once

#include "utils/logger/logger.hpp"
#include "utils/verify.hpp"
#include "sequence/sequence.hpp"

#include <boost/optional.hpp>

namespace io {
class SingleRead;
class SingleReadSeq;
}

namespace omnigraph {

Sequence Subseq(const io::SingleRead& read, size_t start, size_t end);

Sequence Subseq(const io::SingleReadSeq& read, size_t start, size_t end);

Sequence Subseq(const Sequence& s, size_t start, size_t end);

//TODO rename to GapInfo
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
                   boost::optional <Sequence> filling_seq) :
            left_(left),
            right_(right),
            estimated_dist_(estimated_dist),
            left_trim_(left_trim),
            right_trim_(right_trim),
            filling_seq_(filling_seq) {
    }

    auto AsTuple() const {
        return std::make_tuple(left_, right_, left_trim_, right_trim_, estimated_dist_, filling_seq_);
    }

public:
    static const int INVALID_GAP = std::numeric_limits<int>::min();
    //static const GapDescription INVALID_GAP_DESC; //= GapDescription();

    GapDescription(EdgeId left, EdgeId right,
                   int estimated_dist,
                   size_t left_trim = 0, size_t right_trim = 0) :
            GapDescription(left, right,
                           estimated_dist,
                           left_trim, right_trim,
                           boost::none) {
    }

    GapDescription() : GapDescription(EdgeId(), EdgeId(), INVALID_GAP) {
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
        return static_cast<bool>(filling_seq_);
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

    std::string str(const Graph &g) const {
        std::stringstream s;
        s << "left: " << g.str(left_)
          << "; right: " << g.str(right_)
          << "; estimated distance : " << estimated_dist_
          << "; left trim: " << left_trim_
          << "; right trim: " << right_trim_
          << "; sequence " << (filling_seq_ ? filling_seq_->str() : "no_sequence") << std::endl;
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

//Attempts to "shrink" alignments to make them non-overlapping in the query (delivering non-empty "filling" sequence),
// while remaining within graph edges
template<class Graph, class ReadT>
static GapDescription<Graph> CreateGapInfoTryFixOverlap(const Graph &g, const ReadT &read,
                                                        size_t seq_start, size_t seq_end,
                                                        typename Graph::EdgeId left, size_t left_offset,
                                                        typename Graph::EdgeId right, size_t right_offset) {
    VERIFY(left_offset > 0 && right_offset >= 0 &&
           left_offset <= g.length(left) && right_offset < g.length(right));

    TRACE("Creating gap description");

    //trying to shift on the left edge
    if (seq_start >= seq_end) {
        //+1 is a trick to avoid empty gap sequences
        size_t overlap = seq_start - seq_end + 1;
        TRACE("Overlap of size " << overlap << " detected. Fixing.");
        size_t left_shift = std::min(overlap, left_offset - 1);
        VERIFY(seq_start >= left_shift);
        seq_start -= left_shift;
        left_offset -= left_shift;
    }

    //trying to shift on the right edge
    if (seq_start >= seq_end) {
        //+1 is a trick to avoid empty gap sequences
        const size_t overlap = seq_start - seq_end + 1;
        TRACE("Overlap of size " << overlap << " remained. Fixing.");
        size_t right_shift = std::min(overlap, g.length(right) - right_offset - 1);
        VERIFY(seq_end + right_shift <= read.size());
        seq_end += right_shift;
        right_offset += right_shift;
    }

    if (seq_start < seq_end) {
        auto gap_seq = Subseq(read, seq_start, seq_end);
        if (!gap_seq.empty()) {
            TRACE("Gap info successfully created");
            VERIFY(left_offset > 0 && right_offset >= 0 &&
                   left_offset <= g.length(left) && right_offset < g.length(right));
            return GapDescription<Graph>(left, right,
                                  gap_seq,
                                  g.length(left) - left_offset,
                                  right_offset);
        } else {
            TRACE("Something wrong with read subsequence");
        }
    } else {
        TRACE("Failed to fix overlap of size " << seq_start - seq_end + 1);
    }
    return GapDescription<Graph>();
}

}
