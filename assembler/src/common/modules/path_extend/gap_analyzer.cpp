//***************************************************************************
//* Copyright (c) 2014-2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "gap_analyzer.hpp"
#include "overlap_analysis.hpp"

namespace path_extend {

using namespace debruijn_graph;

static size_t HammingDistance(const Sequence& s1, const Sequence& s2) {
    VERIFY(s1.size() == s2.size());
    size_t dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i]) {
            dist++;
        }
    }
    return dist;
}

double HammingGapAnalyzer::ScoreGap(const Sequence& s1, const Sequence& s2) const {
    VERIFY(s1.size() == s2.size());
    return 1.0 - (double) HammingDistance(s1, s2) / (double) s1.size();
}

GapDescription HammingGapAnalyzer::FixGap(const GapDescription &gap) const {
    VERIFY_MSG(gap.no_trim(), "Trims not supported yet");

    size_t max_overlap = basic_overlap_length_;
    if (gap.estimated_dist() < 0) {
        max_overlap -= gap.estimated_dist();
    }

    max_overlap = std::min(max_overlap, g_.k() + std::min(g_.length(gap.left()), g_.length(gap.right())));

    DEBUG("Corrected max overlap " << max_overlap);

    double best_score = min_gap_score_;
    int fixed_gap = GapDescription::INVALID_GAP;
    
    size_t min_overlap = 1;
    if (gap.estimated_dist() < 0) {
        min_overlap = std::max(min_overlap, size_t(math::round(MIN_OVERLAP_COEFF * double(-gap.estimated_dist()))));
    }
    //todo better usage of estimated overlap
    DEBUG("Min overlap " << min_overlap);

    for (size_t l = max_overlap; l >= min_overlap; --l) {
        //TRACE("Sink: " << g_.EdgeNucls(sink).Subseq(g_.length(sink) + g_.k() - l).str());
        //TRACE("Source: " << g_.EdgeNucls(source).Subseq(0, l));
        double score = 0;
        score = ScoreGap(g_.EdgeNucls(gap.left()).Subseq(g_.length(gap.left()) + g_.k() - l),
                         g_.EdgeNucls(gap.right()).Subseq(0, l));
        if (math::gr(score, best_score)) {
            TRACE("Curr overlap " << l);
            TRACE("Score: " << score);
            best_score = score;
            fixed_gap = -int(l);
        }

        if (l == short_overlap_threshold_ && fixed_gap != GapDescription::INVALID_GAP) {
            //look at "short" overlaps only if long overlaps couldn't be found
            DEBUG("Not looking at short overlaps");
            break;
        }
    }
    
    if (fixed_gap != INVALID_GAP) {
        DEBUG("Found candidate gap length with score " << best_score);
        DEBUG("Estimated gap: " << gap.estimated_dist() <<
              ", fixed gap: " << fixed_gap << " (overlap " << (-fixed_gap) << ")");
        
        auto answer = gap;
        answer.set_estimated_dist(fixed_gap);
        return answer;
    } else {
        return GapDescription();
    }
}

GapDescription LAGapAnalyzer::FixGap(const GapDescription &gap) const {
    VERIFY_MSG(gap.no_trim(), "Trims not supported yet");
    //estimated_gap is in k-mers

    size_t estimated_overlap = gap.estimated_dist() < 0 ? size_t(abs(gap.estimated_dist())) : 0;
    DEBUG("SW analyzer");
    DEBUG(size_t(math::round(double(estimated_overlap) * ESTIMATED_GAP_MULTIPLIER))
          + GAP_ADDITIONAL_COEFFICIENT);

    SWOverlapAnalyzer overlap_analyzer(size_t(math::round(double(estimated_overlap) * ESTIMATED_GAP_MULTIPLIER))
                                       + GAP_ADDITIONAL_COEFFICIENT);
    
    auto overlap_info = overlap_analyzer.AnalyzeOverlap(g_, gap.left(), gap.right());
    DEBUG(overlap_info);
    
    if (overlap_info.size() < min_la_length_) {
        DEBUG("Low alignment size");
        DEBUG(min_la_length_);
        DEBUG(min_la_length_);
        return GapDescription();
    }

    size_t max_flank_length = std::max(overlap_info.r2.start_pos,
                                       g_.length(gap.left()) + g_.k() - overlap_info.r1.end_pos);
    DEBUG("Max flank length - " << max_flank_length);
    
    if (int(math::round(double(max_flank_length) * flank_multiplication_coefficient_))
        + flank_addition_coefficient_ > int(overlap_info.size())) {
        DEBUG("Too long flanks for such alignment");
        return GapDescription();
    }

    if (math::ls(overlap_info.identity(), IDENTITY_RATIO)) {
        DEBUG("Low identity score");
        return GapDescription();
    }

    if (overlap_info.r1.end_pos <= g_.k() || overlap_info.r2.start_pos >= g_.length(gap.right())) {
        DEBUG("Less than k+1 nucleotides were left of one of the edges");
        return GapDescription();
    }

    //TODO Is it ok to have a non-symmetric overlap gap description
    return GapDescription(gap.left(), gap.right(),
                          -int(overlap_info.r2.size()),
                          g_.length(gap.left()) + g_.k() - overlap_info.r1.end_pos,
                          overlap_info.r2.start_pos);
}

GapDescription CompositeGapAnalyzer::FixGap(const GapDescription &gap) const {
    VERIFY_MSG(gap.right_trim() == 0 && gap.left_trim() == 0, "Not supported yet");
    DEBUG("Trying to fix estimated gap " << gap.estimated_dist() <<
          " between " << g_.str(gap.left()) << " and " << g_.str(gap.right()));
    
    if (gap.estimated_dist() > int(may_overlap_threshold_)) {
        DEBUG("Edges are supposed to be too far to check overlaps");
        return gap;
    }

    for (auto joiner : joiners_) {
        GapDescription fixed_gap = joiner->FixGap(gap);
        if (fixed_gap != GapDescription()) {
            return fixed_gap;
        }
    }

    //couldn't find decent overlap
    if (gap.estimated_dist() < must_overlap_threshold_) {
        DEBUG("Estimated gap looks unreliable");
        return GapDescription();
    } else {
        DEBUG("Overlap was not found");
        auto answer = gap;
        answer.set_estimated_dist(std::max(gap.estimated_dist(), int(artificial_gap_)));
        return answer;
    }
}

}
