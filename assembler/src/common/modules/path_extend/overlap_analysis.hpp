#pragma once

#include "utils/logger/logger.hpp"
#include "sequence/range.hpp"
#include "ssw/ssw_cpp.h"

namespace debruijn_graph {

struct OverlapInfo {
    Range r1;
    Range r2;
    size_t match_cnt;

    OverlapInfo(const Range& r1_, const Range& r2_, size_t match_cnt_)
            : r1(r1_),
              r2(r2_),
              match_cnt(match_cnt_) {
        VERIFY(match_cnt <= std::min(r1.size(), r2.size()));
    }

    OverlapInfo()
            : match_cnt(0) {
    }

    double identity() const {
        if (match_cnt == 0)
            return 0.;
        return (double)match_cnt / (double)size();
    }

    size_t size() const {
        return std::max(r1.size(), r2.size());
    }

    bool operator==(const OverlapInfo &that) const {
        return r1 == that.r1 && r2 == that.r2 && match_cnt == that.match_cnt;
    }

    bool operator!=(const OverlapInfo &that) const {
        return !(*this == that);
    }
};

inline std::ostream& operator<<(std::ostream& os, const OverlapInfo& info) {
    return os << "R1: [" << info.r1.start_pos << ", " << info.r1.end_pos
            << "]; R2: [" << info.r2.start_pos << ", " << info.r2.end_pos << "]"
            << "; match_cnt: " << info.match_cnt;
}

class SWOverlapAnalyzer {
    static const uint32_t CIGAR_FLAG_MASK = (1 << 4) - 1;
    static const uint32_t CIGAR_MATCH_FLAG = 7;
    typedef typename Graph::EdgeId EdgeId;
    size_t flank_length_;

    const StripedSmithWaterman::Aligner aligner_;
    const StripedSmithWaterman::Filter filter_;

    size_t CountMatches(std::vector<uint32_t> cigar) const {
        size_t match_cnt = 0;
        for (uint32_t entry : cigar) {
            if ((entry & CIGAR_FLAG_MASK) == CIGAR_MATCH_FLAG) {
                match_cnt += (entry >> 4);
            }
        }
        return match_cnt;
    }

    OverlapInfo InnerAnalyze(const Sequence& s1, const Sequence& s2) const {
        if (s1.size() == 0 || s2.size() == 0) {
            return OverlapInfo();
        }
        StripedSmithWaterman::Alignment alignment;
        if (aligner_.Align(s1.str().c_str(), s2.str().c_str(), int(s2.size()), filter_, &alignment)) {
            if (alignment.sw_score > 0) {
                return OverlapInfo(Range(alignment.query_begin, alignment.query_end + 1),
                            Range(alignment.ref_begin, alignment.ref_end + 1),
                            CountMatches(alignment.cigar));
            }
        }
        return OverlapInfo();
    }

public:
    SWOverlapAnalyzer(size_t flank_length)
            : flank_length_(flank_length),
              aligner_(/*match_score*/1,
              /*mismatch_penalty*/3,
                       /*gap_opening_penalty*/4,
                       /*gap_extending_penalty*/3) {
        DEBUG("Considered max overlap " << flank_length);
    }


    OverlapInfo AnalyzeOverlap(const Sequence& s1, const Sequence& s2) const {
        DEBUG("Analysis started");
        size_t start1 = flank_length_ > s1.size() ? 0 : s1.size() - flank_length_;
        size_t end2 = flank_length_ > s2.size() ? s2.size() : flank_length_;

        DEBUG("s1 " << s1.Subseq(start1, s1.size()));
        DEBUG("s2 " << s2.Subseq(0, end2));
        OverlapInfo result = InnerAnalyze(s1.Subseq(start1, s1.size()), s2.Subseq(0, end2));
        if (result == OverlapInfo()) {
            DEBUG("Empty overlap")
            return result;
        }

        result.r1.shift(int(start1));
        DEBUG("Result " << result)
        return result;
    }

    template<class Graph>
    OverlapInfo AnalyzeOverlap(const Graph& g, EdgeId e1, EdgeId e2) const {
        DEBUG("Analyzing edges " << g.str(e1) << " and " << g.str(e2));
        return AnalyzeOverlap(g.EdgeNucls(e1), g.EdgeNucls(e2));
    }

private:
    DECL_LOGGER("SWOverlapAnalyzer");
};

}
