//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * path_extender.hpp
 *
 *  Created on: Mar 5, 2012
 *      Author: andrey
 */

#pragma once

#include "extension_chooser.hpp"
#include "common/modules/alignment/gap_info.hpp"
#include "assembly_graph/paths/bidirectional_path_container.hpp"
#include "path_filter.hpp"
#include "overlap_analysis.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include <cmath>

namespace path_extend {

inline BidirectionalPath OptimizedConjugate(const BidirectionalPath &path) {
    return path.GetConjPath() ? *path.GetConjPath() : path.Conjugate();
}

//TODO think about symmetry and what if it breaks?
class OverlapFindingHelper {
    const Graph &g_;
    const GraphCoverageMap &coverage_map_;
    const size_t min_edge_len_;
    const size_t max_diff_;
    const bool try_extend_;

    //TODO think of the cases when (gap + length) < 0
    //Changes second argument on success
    void TryExtendToEnd(const BidirectionalPath &path, size_t &pos) const {
        if (pos < path.Size() &&
                path.GapAt(pos).gap + path.LengthAt(pos) <= max_diff_)
            pos = path.Size();
    }

    //Changes second argument on success
    void TryExtendToStart(const BidirectionalPath &path, size_t &pos) const {
        if (pos > 0 && path.Length() - path.LengthAt(pos) <= max_diff_)
            pos = 0;
    }

    pair<Range, Range> ComparePaths(const BidirectionalPath &path1,
                                    const BidirectionalPath &path2,
                                    size_t start2) const {
        TRACE("Comparing paths " << path1.GetId() << " and " << path2.GetId());
        //TODO change to edit distance?
        int shift1 = 0;
        //path1 is always matched from the start
        const size_t start1 = 0;
        size_t end1 = start1;
        size_t end2 = start2;

        for (size_t i = start1; i < path1.Size(); ++i) {
            if (abs(shift1) > int(max_diff_))
                break;

            bool match = false;
            size_t j = end2;
            int shift2 = 0;
            for (; j < path2.Size(); ++j) {
                if (end1 == 0) {
                    //Force first match to start with pos2
                    if (j > start2) {
                        break;
                    }
                }

                if (abs(shift2) > int(max_diff_))
                    break;
                if (path1.At(i) == path2.At(j) &&
                        (end1 == 0 ||
                            abs(shift1 + path1.GapAt(i).gap - shift2 - path2.GapAt(j).gap) <= int(max_diff_))) {
                    match = true;
                    break;
                } else {
                    shift2 += path2.ShiftLength(j);
                }
            }
            if (match) {
                end1 = i+1;
                end2 = j+1;
                shift1 = 0;
            } else {
                shift1 += path1.ShiftLength(i);
            }
        }

        //Extending the ends of the paths if possible
        if (try_extend_ && end1 > 0) {
            TryExtendToEnd(path1, end1);
            TryExtendToEnd(path2, end2);
            //no need to extend path1 left
            VERIFY(start1 == 0);
            TryExtendToStart(path2, start2);
        }

        return make_pair(Range(start1, end1), Range(start2, end2));
    }

public:
    OverlapFindingHelper(const Graph &g,
                         const GraphCoverageMap &coverage_map,
                         size_t min_edge_len,
                         size_t max_diff) :
            g_(g),
            coverage_map_(coverage_map),
            min_edge_len_(min_edge_len),
            max_diff_(max_diff),
            //had to enable try_extend, otherwise equality lost symmetry
            try_extend_(max_diff_ > 0) {
    }

    bool IsSubpath(const BidirectionalPath &path,
                   const BidirectionalPath &other) const {
        for (size_t j = 0; j < other.Size(); ++j) {
            auto range_pair = ComparePaths(path, other, j);
            if (range_pair.first.end_pos == path.Size()) {
                return true;
            }
        }
        return false;
    }

    //NB! Equality is not transitive if max_diff is > 0
    bool IsEqual(const BidirectionalPath &path,
                 const BidirectionalPath &other) const {
        auto ends_pair = CommonPrefix(path, other);
        return ends_pair.first == path.Size()
               && ends_pair.second == other.Size();
    }


    pair<size_t, size_t> CommonPrefix(const BidirectionalPath &path1,
                                      const BidirectionalPath &path2) const {
        auto answer = make_pair(0, 0);
        size_t cum = 0;
        size_t max_overlap = 0;
        for (size_t j = 0; j < path2.Size(); ++j) {
            auto range_pair = ComparePaths(path1, path2, j);
            if (range_pair.second.start_pos == 0 && range_pair.first.size() > max_overlap) {
                answer = make_pair(range_pair.first.end_pos, range_pair.second.end_pos);
                max_overlap = range_pair.first.size();
            }

            if (!try_extend_)
                break;

            cum += path2.ShiftLength(j);
            if (cum > max_diff_)
                break;
        }
        return answer;
    };

    //overlap is forced to start from the beginning of path1
    pair<Range, Range> FindOverlap(const BidirectionalPath &path1,
                                   const BidirectionalPath &path2,
                                   bool end_start_only) const {
        size_t max_overlap = 0;
        pair<Range, Range> matching_ranges;
        for (size_t j = 0; j < path2.Size(); ++j) {
            auto range_pair = ComparePaths(path1, path2, j);
            VERIFY(range_pair.first.start_pos == 0);
            //checking if overlap is valid
            if (end_start_only && range_pair.second.end_pos != path2.Size())
                continue;

            size_t overlap_size = range_pair.first.size();
            if (overlap_size > max_overlap ||
                //prefer overlaps with end of path2
                (overlap_size == max_overlap &&
                 range_pair.second.end_pos == path2.Size())) {
                max_overlap = overlap_size;
                matching_ranges = range_pair;
            }
        }
        return matching_ranges;
    }

    vector<const BidirectionalPath*> FindCandidatePaths(const BidirectionalPath &path) const {
        set<const BidirectionalPath*> candidates;
        size_t cum_len = 0;
        for (size_t i = 0; i < path.Size(); ++i) {
            if (cum_len > max_diff_)
                break;
            EdgeId e = path.At(i);
            if (g_.length(e) >= min_edge_len_) {
                utils::insert_all(candidates, coverage_map_.GetCoveringPaths(e));
                cum_len += path.ShiftLength(i);
            }
        }
        return vector<const BidirectionalPath*>(candidates.begin(), candidates.end());
    }

private:
    DECL_LOGGER("OverlapFindingHelper");
};

inline void SubscribeCoverageMap(BidirectionalPath * path, GraphCoverageMap &coverage_map) {
    path->Subscribe(&coverage_map);
    for (size_t i = 0; i < path->Size(); ++i) {
        coverage_map.BackEdgeAdded(path->At(i), path, path->GapAt(i));
    }
}

inline BidirectionalPath* AddPath(PathContainer &paths,
                                  const BidirectionalPath &path,
                                  GraphCoverageMap &coverage_map) {
    BidirectionalPath* p = new BidirectionalPath(path);
    BidirectionalPath* conj_p = new BidirectionalPath(OptimizedConjugate(path));
    SubscribeCoverageMap(p, coverage_map);
    SubscribeCoverageMap(conj_p, coverage_map);
    paths.AddPair(p, conj_p);
    return p;
}

class ShortLoopEstimator {
public:
    //Path must end with forward cycle edge, contain at least 2 edges and must not contain backward cycle edges
    //Returns 0 (for no loops), 1 (for a single loop) or 2 (for many loops)
    virtual size_t EstimateSimpleCycleCount(const BidirectionalPath& path, EdgeId backward_edge, EdgeId exit_edge) const = 0;

    virtual ~ShortLoopEstimator() {};
};

class ShortLoopResolver {
public:
    static const size_t BASIC_N_CNT = 100;

    ShortLoopResolver(const Graph& g, shared_ptr<ShortLoopEstimator> loop_estimator)
            : g_(g), loop_estimator_(loop_estimator) { }

    void ResolveShortLoop(BidirectionalPath& path) const {
        EdgeId back_cycle_edge;
        EdgeId loop_outgoing;
        EdgeId loop_incoming;
        if (path.Size() >=1 && GetLoopAndExit(g_, path.Back(), back_cycle_edge, loop_outgoing, loop_incoming)) {
            DEBUG("Resolving short loop...");
            MakeBestChoice(path, back_cycle_edge, loop_outgoing, loop_incoming);
            DEBUG("Resolving short loop done");
        }
    }

private:
    DECL_LOGGER("PathExtender")
    const Graph& g_;
    shared_ptr<ShortLoopEstimator> loop_estimator_;

    void UndoCycles(BidirectionalPath& p, EdgeId back_cycle_edge, EdgeId loop_incoming) const {
        if (p.Size() <= 2) {
            return;
        }
        EdgeId forward_cycle_edge = p.Back();
        size_t loop_start_index = p.Size();
        while (loop_start_index > 2) {
            if (p.At(loop_start_index - 1) == forward_cycle_edge && p.At(loop_start_index - 2) == back_cycle_edge) {
                loop_start_index -= 2;
            } else {
                break;
            }
        }

        if (p.At(loop_start_index - 1) == forward_cycle_edge) {
            p.PopBack(p.Size() - loop_start_index);
        } else if (loop_start_index != p.Size()) {
            //Means we jumped into the loop
            DEBUG("Jumped inside the loop, back loop edge " << g_.int_id(back_cycle_edge) << ", forward loop edge " << g_.int_id(forward_cycle_edge) << ", loop starts @ " << loop_start_index);
            p.PrintDEBUG();
            Gap gap = p.GapAt(loop_start_index);
            p.PopBack(p.Size() - loop_start_index);

            if (p.Back() != loop_incoming) {
                p.PushBack(loop_incoming, Gap(max(0, gap.gap - (int) g_.length(loop_incoming)  - (int) g_.length(forward_cycle_edge)), {gap.trash.previous, 0}));
            }
            p.PushBack(forward_cycle_edge);
            DEBUG("Restored the path");
            p.PrintDEBUG();
        }
    }

    //edges -- first edge is loop's back edge, second is loop exit edge
    void MakeBestChoice(BidirectionalPath& path, EdgeId back_cycle_edge, EdgeId loop_outgoing, EdgeId loop_incoming) const {
        EdgeId forward_cycle_edge = path.Back();
        UndoCycles(path, back_cycle_edge, loop_incoming);

        //Expects 0 (for no loops), 1 (for a single loop) or 2 (for many loops, will insert back_cycle_edge and Ns)
        size_t loop_count = loop_estimator_->EstimateSimpleCycleCount(path, back_cycle_edge, loop_outgoing);
        if (loop_count > 0) {
            path.PushBack(back_cycle_edge);
            if (loop_count == 1) {
                DEBUG("Single loop");
                path.PushBack(forward_cycle_edge);
                path.PushBack(loop_outgoing);
            }
            else {
                DEBUG("Multiple cycles");
                //If the forward edge is shorter than K, avoid overlapping bases between backward edge and outgoing edge
                //Make sure that the N-stretch will be exactly 100 bp
                uint32_t overlapping_bases = (uint32_t) max(int(g_.k()) - int(g_.length(forward_cycle_edge)), 0);
                path.PushBack(loop_outgoing, Gap(int(g_.k() + BASIC_N_CNT - overlapping_bases), {0, overlapping_bases}));
            }
        }
        else {
            path.PushBack(loop_outgoing);
        }
    }
};

class CoverageLoopEstimator : public ShortLoopEstimator {
public:
    CoverageLoopEstimator(const Graph& g, const FlankingCoverage<Graph>& flanking_cov)
            : g_(g), flanking_cov_(flanking_cov) {

    }

    //Path must end with forward cycle edge, contain at least 2 edges and must not contain backward cycle edges
    //Returns 0 (for no loops), 1 (for a single loop) or 2 (for many loops)
    size_t EstimateSimpleCycleCount(const BidirectionalPath& path, EdgeId backward_edge, EdgeId exit_edge) const override {
        VERIFY(path.Size() > 1);
        EdgeId forward_edge = path.Back();
        EdgeId incoming_edge = path[path.Size() - 2];
        double in_cov = flanking_cov_.GetOutCov(incoming_edge);
        double out_cov = flanking_cov_.GetInCov(exit_edge);
        double avg_coverage = (in_cov + out_cov) / 2.0;

        double fwd_count = math::round(g_.coverage(forward_edge) / avg_coverage);
        double back_count = math::round(g_.coverage(backward_edge) / avg_coverage);
        size_t result = (size_t) math::round(std::max(0.0, std::min(fwd_count - 1.0, back_count)));

        DEBUG("loop with start " << g_.int_id(incoming_edge)
                <<" e1 " << g_.int_id(forward_edge)
                << " e2 " << g_.int_id(backward_edge)
                << " out " <<g_.int_id(exit_edge)
                << " cov in = " << in_cov
                << " cov out " << out_cov
                << " cov " << avg_coverage
                << " cov e1 = " << g_.coverage(forward_edge)
                << " cov e2 = " << g_.coverage(backward_edge)
                << " fwd_count = " << fwd_count
                << " back_count = " << back_count
                << " result = " << result);

        return result;
    }

private:
    const Graph& g_;
    const FlankingCoverage<Graph>& flanking_cov_;
};

class PairedInfoLoopEstimator: public ShortLoopEstimator {
    const Graph& g_;
    shared_ptr<WeightCounter> wc_;
    double weight_threshold_;

public:
    PairedInfoLoopEstimator(const Graph& g, shared_ptr<WeightCounter> wc, double weight_threshold = 0.0)
            : g_(g),
              wc_(wc),
              weight_threshold_(weight_threshold) { }

    //Path must end with forward cycle edge, contain at least 2 edges and must not contain backward cycle edges
    //Returns 0 (for no loops), 1 (for a single loop) or 2 (for many loops)
    size_t EstimateSimpleCycleCount(const BidirectionalPath& path, EdgeId backward_edge, EdgeId /*exit_edge*/) const override {
        VERIFY(path.Size() > 1);
        VERIFY(wc_ != nullptr);
        EdgeId forward_cycle_edge = path.Back();

        size_t result = 0;
        double lopp_edge_weight = wc_->CountWeight(path, backward_edge);
        if (math::gr(lopp_edge_weight, weight_threshold_)) {
            //Paired information on loop back edges exits => at leat one iteration
            //Looking for paired information supporting more than 1 cycle
            if (NoSelfPairedInfo(backward_edge, forward_cycle_edge)) {
                //More likely to be a single cycle
                DEBUG("Single loop");
                result = 1;
            }
            else {
                DEBUG("Multiple cycles");
                //More likely to be a 2 or more cycles
                result = 2;
            }
        }
        return result;
    }

private:

    bool NoSelfPairedInfo(EdgeId back_cycle_edge, EdgeId forward_cycle_edge) const {
        size_t is = wc_->PairedLibrary().GetISMax();
        int forward_len = (int) g_.length(forward_cycle_edge);
        bool exists_pi = true;

        BidirectionalPath cycle(g_, back_cycle_edge);
        while (cycle.Length() < is + g_.length(back_cycle_edge)) {
            auto w = wc_->CountWeight(cycle, back_cycle_edge, std::set<size_t>(), forward_len);
            if (math::gr(w, weight_threshold_)) {
                //Paired information found within loop
                DEBUG("Found PI with back weight " << w << ", weight threshold " << weight_threshold_);
                exists_pi = false;
                break;
            }
            cycle.PushBack(back_cycle_edge, Gap(forward_len));
        }

        return exists_pi;
    }
};

class CombinedLoopEstimator: public ShortLoopEstimator {
public:
    CombinedLoopEstimator(const Graph& g,
                          const FlankingCoverage<Graph>& flanking_cov,
                          shared_ptr<WeightCounter> wc,
                          double weight_threshold = 0.0)
        : pi_estimator_(g, wc, weight_threshold),
          cov_estimator_(g, flanking_cov) {}

    //Path must end with forward cycle edge, contain at least 2 edges and must not contain backward cycle edges
    //Returns 0 (for no loops), 1 (for a single loop) or 2 (for many loops)
    size_t EstimateSimpleCycleCount(const BidirectionalPath& path, EdgeId backward_edge, EdgeId exit_edge) const override {
        size_t result = pi_estimator_.EstimateSimpleCycleCount(path, backward_edge, exit_edge);
        if (result == 1) {
            //Verify using coverage
            if (cov_estimator_.EstimateSimpleCycleCount(path, backward_edge, exit_edge) > 1)
                result = 2;
        }
        return result;
    }

private:
    PairedInfoLoopEstimator pi_estimator_;
    CoverageLoopEstimator cov_estimator_;
};



//TODO move to gap_closing.hpp
typedef omnigraph::GapDescription<Graph> GapDescription;
class GapAnalyzer {

public:
    static const int INVALID_GAP = GapDescription::INVALID_GAP;
    GapAnalyzer(const Graph& g)
            : g_(g) { }

    virtual GapDescription FixGap(const GapDescription &gap) const = 0;

    virtual ~GapAnalyzer() { }
protected:
    const Graph& g_;
};

class HammingGapAnalyzer: public GapAnalyzer {
    const double min_gap_score_;
    const size_t short_overlap_threshold_;
    const size_t basic_overlap_length_;

    static constexpr double MIN_OVERLAP_COEFF = 0.05;

    size_t HammingDistance(const Sequence& s1, const Sequence& s2) const {
        VERIFY(s1.size() == s2.size());
        size_t dist = 0;
        for (size_t i = 0; i < s1.size(); ++i) {
            if (s1[i] != s2[i]) {
                dist++;
            }
        }
        return dist;
    }

    double ScoreGap(const Sequence& s1, const Sequence& s2) const {
        VERIFY(s1.size() == s2.size());
        return 1.0 - (double) HammingDistance(s1, s2) / (double) s1.size();
    }

public:

    //todo review parameters in usages
    HammingGapAnalyzer(const Graph& g,
            double min_gap_score,
            size_t short_overlap_threshold,
            size_t basic_overlap_length):
                GapAnalyzer(g),
                min_gap_score_(min_gap_score),
                short_overlap_threshold_(short_overlap_threshold),
                basic_overlap_length_(basic_overlap_length)
    {
        DEBUG("HammingGapAnalyzer params: \n min_gap_score " << min_gap_score_ <<
              "\n short_overlap_threshold " << short_overlap_threshold_ <<
              "\n basic_overlap_length " << basic_overlap_length_);
    }

    GapDescription FixGap(const GapDescription &gap) const override {
        VERIFY_MSG(gap.no_trim(), "Trims not supported yet");

        size_t max_overlap = basic_overlap_length_;
        if (gap.estimated_dist() < 0) {
            max_overlap -= gap.estimated_dist();
        }

        max_overlap = min(max_overlap,
                                      g_.k() + min(g_.length(gap.left()), g_.length(gap.right())));

        DEBUG("Corrected max overlap " << max_overlap);

        double best_score = min_gap_score_;
        int fixed_gap = GapDescription::INVALID_GAP;

        size_t min_overlap = 1;
        if (gap.estimated_dist() < 0) {
            min_overlap = max(min_overlap, size_t(math::round(MIN_OVERLAP_COEFF * double(-gap.estimated_dist()))));
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

private:
    DECL_LOGGER("HammingGapAnalyzer");
};

//LA stands for Local Alignment
//TODO if current setting will work -- get rid of flank_*_coefficient params
class LAGapAnalyzer: public GapAnalyzer {
public:
    LAGapAnalyzer(const Graph& g, size_t min_la_length,
            double flank_multiplication_coefficient,
            int flank_addition_coefficient) :
            GapAnalyzer(g),
            min_la_length_(min_la_length),
            flank_multiplication_coefficient_(flank_multiplication_coefficient),
            flank_addition_coefficient_(flank_addition_coefficient) {
        DEBUG("flank_multiplication_coefficient - " << flank_multiplication_coefficient_);
        DEBUG("flank_addition_coefficient  - " << flank_addition_coefficient_ );
    }

    GapDescription FixGap(const GapDescription &gap) const override {
        VERIFY_MSG(gap.no_trim(), "Trims not supported yet");
        //estimated_gap is in k-mers

        size_t estimated_overlap = gap.estimated_dist() < 0 ? size_t(abs(gap.estimated_dist())) : 0;
        SWOverlapAnalyzer overlap_analyzer(size_t(math::round(double(estimated_overlap) * ESTIMATED_GAP_MULTIPLIER))
                                           + GAP_ADDITIONAL_COEFFICIENT);

        auto overlap_info = overlap_analyzer.AnalyzeOverlap(g_, gap.left(), gap.right());
        DEBUG(overlap_info);

        if (overlap_info.size() < min_la_length_) {
            DEBUG("Low alignment size");
            return GapDescription();
        }

        size_t max_flank_length = max(overlap_info.r2.start_pos,
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

private:
    DECL_LOGGER("LAGapAnalyzer");
    const size_t min_la_length_;
    const double flank_multiplication_coefficient_;
    const int flank_addition_coefficient_;

    static constexpr double IDENTITY_RATIO = 0.9;
    static constexpr double ESTIMATED_GAP_MULTIPLIER = 2.0;
    static constexpr size_t GAP_ADDITIONAL_COEFFICIENT = 30;
};


class CompositeGapAnalyzer: public GapAnalyzer {
public:

    CompositeGapAnalyzer(const Graph& g,
                       const vector<shared_ptr<GapAnalyzer>>& joiners,
                       size_t may_overlap_threshold,
                       int must_overlap_threshold,
                       size_t artificial_gap) :
            GapAnalyzer(g),
            joiners_(joiners),
            may_overlap_threshold_(may_overlap_threshold),
            must_overlap_threshold_(must_overlap_threshold),
            artificial_gap_(artificial_gap)
    {  }

    GapDescription FixGap(const GapDescription &gap) const override {
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
            answer.set_estimated_dist(max(gap.estimated_dist(), int(artificial_gap_)));
            return answer;
        }
    }

private:
    vector<shared_ptr<GapAnalyzer>> joiners_;
    const size_t may_overlap_threshold_;
    const int must_overlap_threshold_;
    const size_t artificial_gap_;

    DECL_LOGGER("CompositeGapAnalyzer");
};

//Detects a cycle as a minsuffix > IS present earlier in the path. Overlap is allowed.
class InsertSizeLoopDetector {
protected:
    GraphCoverageMap visited_cycles_coverage_map_;
    PathContainer path_storage_;
    size_t min_cycle_len_;

public:
    InsertSizeLoopDetector(const Graph& g, size_t is):
        visited_cycles_coverage_map_(g),
        path_storage_(),
        min_cycle_len_(is) {
    }

    ~InsertSizeLoopDetector() {
        path_storage_.DeleteAllPaths();
    }

    bool CheckCycledNonIS(const BidirectionalPath& path) const {
        if (path.Size() <= 2) {
            return false;
        }
        BidirectionalPath last = path.SubPath(path.Size() - 2);
        int pos = path.FindFirst(last);
        VERIFY(pos >= 0);
        return size_t(pos) != path.Size() - 2;
    }

    bool CheckCycled(const BidirectionalPath& path) const {
        return FindCycleStart(path) != -1;
    }
//first suffix longer than min_cycle_len
    int FindPosIS(const BidirectionalPath& path) const {
        int i = (int) path.Size() - 1;
        while (i >= 0 && path.LengthAt(i) < min_cycle_len_) {
            --i;
        }
        return i;
    }
    int FindCycleStart(const BidirectionalPath& path) const {
        TRACE("Looking for IS cycle " << min_cycle_len_);
        int i = FindPosIS(path);
        TRACE("last is pos " << i);
        if (i < 0) return -1;
//Tail
        BidirectionalPath last = path.SubPath(i);
        //last.Print();

        int pos = path.FindFirst(last);
// not cycle
        if (pos == i) pos = -1;
        TRACE("looking for 1sr IS cycle " << pos);
        return pos;
    }

//After cycle detected, removes min suffix > IS.
//returns the beginning of the cycle.
    int RemoveCycle(BidirectionalPath& path) const {
        int pos = FindCycleStart(path);
        DEBUG("Found IS cycle " << pos);
        if (pos == -1) {
            return -1;
        }

        int last_edge_pos = FindPosIS(path);
        VERIFY(last_edge_pos > -1);
        DEBUG("last edge pos " << last_edge_pos);
        VERIFY(last_edge_pos > pos);
        for (int i = (int) path.Size() - 1; i >= last_edge_pos; --i) {
            path.PopBack();
        }
        VERIFY((int) path.Size() == last_edge_pos);
        VERIFY(pos < (int) path.Size());
        DEBUG("result pos " <<pos);
        return pos;
    }

    //seems that it is outofdate
    bool InExistingLoop(const BidirectionalPath& path) {
        DEBUG("Checking existing loops");
        auto visited_cycles = visited_cycles_coverage_map_.GetEdgePaths(path.Back());
        for (auto cycle : *visited_cycles) {
            DEBUG("checking  cycle ");
            int pos = path.FindLast(*cycle);
            if (pos == -1)
                continue;

            int start_cycle_pos = pos + (int) cycle->Size();
            bool only_cycles_in_tail = true;
            int last_cycle_pos = start_cycle_pos;
            DEBUG("start_cycle pos "<< last_cycle_pos);
            for (int i = start_cycle_pos; i < (int) path.Size() - (int) cycle->Size(); i += (int) cycle->Size()) {
                if (!path.CompareFrom(i, *cycle)) {
                    only_cycles_in_tail = false;
                    break;
                } else {
                    last_cycle_pos = i + (int) cycle->Size();
                    DEBUG("last cycle pos changed " << last_cycle_pos);
                }
            }
            DEBUG("last_cycle_pos " << last_cycle_pos);
            only_cycles_in_tail = only_cycles_in_tail && cycle->CompareFrom(0, path.SubPath(last_cycle_pos));
            if (only_cycles_in_tail) {
// seems that most of this is useless, checking
                VERIFY (last_cycle_pos == start_cycle_pos);
                DEBUG("find cycle " << last_cycle_pos);
                DEBUG("path");
                path.PrintDEBUG();
                DEBUG("last subpath");
                path.SubPath(last_cycle_pos).PrintDEBUG();
                DEBUG("cycle");
                cycle->PrintDEBUG();
                DEBUG("last_cycle_pos " << last_cycle_pos << " path size " << path.Size());
                VERIFY(last_cycle_pos <= (int)path.Size());
                DEBUG("last cycle pos + cycle " << last_cycle_pos + (int)cycle->Size());
                VERIFY(last_cycle_pos + (int)cycle->Size() >= (int)path.Size());

                return true;
            }
        }
        return false;
    }

    void AddCycledEdges(const BidirectionalPath& path, size_t pos) {
        if (pos >= path.Size()) {
            DEBUG("Wrong position in IS cycle");
            return;
        }
        BidirectionalPath * p = new BidirectionalPath(path.SubPath(pos));
        BidirectionalPath * cp = new BidirectionalPath(p->Conjugate());
        visited_cycles_coverage_map_.Subscribe(p);
        visited_cycles_coverage_map_.Subscribe(cp);
        DEBUG("add cycle");
        p->PrintDEBUG();
    }
};

class PathExtender {
public:
    explicit PathExtender(const Graph &g):
        g_(g) { }

    virtual ~PathExtender() { }

    virtual bool MakeGrowStep(BidirectionalPath& path, PathContainer* paths_storage = nullptr) = 0;

protected:
    const Graph &g_;
    DECL_LOGGER("PathExtender")
};

class CompositeExtender {
public:

    CompositeExtender(const Graph &g, GraphCoverageMap& cov_map,
                      UsedUniqueStorage &unique,
                      const vector<shared_ptr<PathExtender>> &pes)
            : g_(g),
              cover_map_(cov_map),
              used_storage_(unique),
              extenders_(pes) {}

    void GrowAll(PathContainer& paths, PathContainer& result) {
        result.clear();
        GrowAllPaths(paths, result);
        result.FilterEmptyPaths();
    }

    void GrowPath(BidirectionalPath& path, PathContainer* paths_storage) {
        while (MakeGrowStep(path, paths_storage)) { }
    }


private:
    const Graph &g_;
    GraphCoverageMap &cover_map_;
    UsedUniqueStorage &used_storage_;
    vector<shared_ptr<PathExtender>> extenders_;

    bool MakeGrowStep(BidirectionalPath& path, PathContainer* paths_storage) {
        DEBUG("make grow step composite extender");

        size_t current = 0;
        while (current < extenders_.size()) {
            DEBUG("step " << current << " of total " << extenders_.size());
            if (extenders_[current]->MakeGrowStep(path, paths_storage)) {
                return true;
            }
           ++current;
        }
        return false;
    }
    
    void GrowAllPaths(PathContainer& paths, PathContainer& result) {
        for (size_t i = 0; i < paths.size(); ++i) {
            VERBOSE_POWER_T2(i, 100, "Processed " << i << " paths from " << paths.size() << " (" << i * 100 / paths.size() << "%)");
            if (paths.size() > 10 && i % (paths.size() / 10 + 1) == 0) {
                INFO("Processed " << i << " paths from " << paths.size() << " (" << i * 100 / paths.size() << "%)");
            }
            //In 2015 modes do not use a seed already used in paths.
            //FIXME what is the logic here?
            if (used_storage_.UniqueCheckEnabled()) {
                bool was_used = false;
                for (size_t ind =0; ind < paths.Get(i)->Size(); ind++) {
                    EdgeId eid = paths.Get(i)->At(ind);
                    if (used_storage_.IsUsedAndUnique(eid)) {
                        DEBUG("Used edge " << g_.int_id(eid));
                        was_used = true;
                        break;
                    } else {
                        used_storage_.insert(eid);
                    }
                }
                if (was_used) {
                    DEBUG("skipping already used seed");
                    continue;
                }
            }

            if (!cover_map_.IsCovered(*paths.Get(i))) {
                AddPath(result, *paths.Get(i), cover_map_);
                BidirectionalPath * path = new BidirectionalPath(*paths.Get(i));
                BidirectionalPath * conjugatePath = new BidirectionalPath(*paths.GetConjugate(i));
                SubscribeCoverageMap(path, cover_map_);
                SubscribeCoverageMap(conjugatePath, cover_map_);
                result.AddPair(path, conjugatePath);
                size_t count_trying = 0;
                size_t current_path_len = 0;
                do {
                    current_path_len = path->Length();
                    count_trying++;
                    GrowPath(*path, &result);
                    GrowPath(*conjugatePath, &result);
                } while (count_trying < 10 && (path->Length() != current_path_len));
                DEBUG("result path " << path->GetId());
                path->PrintDEBUG();
            }
        }
    }

};

//All Path-Extenders inherit this one
class LoopDetectingPathExtender : public PathExtender {
    const bool use_short_loop_cov_resolver_;
    ShortLoopResolver cov_loop_resolver_;

    InsertSizeLoopDetector is_detector_;
    UsedUniqueStorage &used_storage_;

protected:
    const bool investigate_short_loops_;
    const GraphCoverageMap &cov_map_;

    bool TryUseEdge(BidirectionalPath &path, EdgeId e, const Gap &gap) {
        bool success = used_storage_.TryUseEdge(path, e, gap);
        if (success) {
            DEBUG("Adding edge. PathId: " << path.GetId() << " path length: " << path.Length() - 1 << ", fixed gap : "
                                          << gap.gap << ", trash length: " << gap.trash.previous << "-" << gap.trash.current);
        }
        return success;
    }

    bool DetectCycle(BidirectionalPath& path) {
        DEBUG("detect cycle");
        if (is_detector_.CheckCycled(path)) {
            DEBUG("Checking IS cycle");
            int loop_pos = is_detector_.RemoveCycle(path);
            DEBUG("Removed IS cycle");
            if (loop_pos != -1) {
                is_detector_.AddCycledEdges(path, loop_pos);
                return true;
            }
        }
        return false;
    }

    bool DetectCycleScaffolding(BidirectionalPath& path, EdgeId e) {
        BidirectionalPath temp_path(path);
        temp_path.PushBack(e);
        return is_detector_.CheckCycledNonIS(temp_path);
    }

    virtual bool MakeSimpleGrowStep(BidirectionalPath& path, PathContainer* paths_storage = nullptr) = 0;

    virtual bool ResolveShortLoopByCov(BidirectionalPath& path) {
        LoopDetector loop_detector(&path, cov_map_);
        size_t init_len = path.Length();
        bool result = false;
        while (path.Size() >= 1 && loop_detector.EdgeInShortLoop(path.Back())) {
            cov_loop_resolver_.ResolveShortLoop(path);
            if (init_len == path.Length()) {
                return result;
            } else {
                result = true;
            }
            init_len = path.Length();
        }
        return true;
    }

    virtual bool ResolveShortLoopByPI(BidirectionalPath& path) = 0;

    virtual bool CanInvestigateShortLoop() const {
        return false;
    }

public:
    LoopDetectingPathExtender(const conj_graph_pack &gp,
                              const GraphCoverageMap &cov_map,
                              UsedUniqueStorage &unique,
                              bool investigate_short_loops,
                              bool use_short_loop_cov_resolver,
                              size_t is)
            : PathExtender(gp.g),
              use_short_loop_cov_resolver_(use_short_loop_cov_resolver),
              cov_loop_resolver_(gp.g, make_shared<CoverageLoopEstimator>(gp.g, gp.flanking_cov)),
              is_detector_(gp.g, is),
              used_storage_(unique),
              investigate_short_loops_(investigate_short_loops),
              cov_map_(cov_map) {

    }


    bool MakeGrowStep(BidirectionalPath& path, PathContainer* paths_storage) override {
        if (is_detector_.InExistingLoop(path)) {
            DEBUG("in existing loop");
            return false;
        }
        DEBUG("un ch enabled " << used_storage_.UniqueCheckEnabled());
        bool result;
        LoopDetector loop_detector(&path, cov_map_);
        if (DetectCycle(path)) {
            result = false;
        } else if (path.Size() >= 1 && InvestigateShortLoop() && loop_detector.EdgeInShortLoop(path.Back()) && use_short_loop_cov_resolver_) {
            DEBUG("edge in short loop");
            result = ResolveShortLoop(path);
        } else if (InvestigateShortLoop() && loop_detector.PrevEdgeInShortLoop() && use_short_loop_cov_resolver_) {
            DEBUG("Prev edge in short loop");
            path.PopBack();
            result = ResolveShortLoop(path);
        } else {
            DEBUG("Making step");
            result = MakeSimpleGrowStep(path, paths_storage);
            DEBUG("Made step");
            if (DetectCycle(path)) {
                result = false;
            } else if (path.Size() >= 1 && InvestigateShortLoop() && loop_detector.EdgeInShortLoop(path.Back())) {
                DEBUG("Edge in short loop");
                result = ResolveShortLoop(path);
            } else if (InvestigateShortLoop() && loop_detector.PrevEdgeInShortLoop()) {
                DEBUG("Prev edge in short loop");
                path.PopBack();
                result = ResolveShortLoop(path);
            }
        }
        return result;
    }

private:
    bool ResolveShortLoop(BidirectionalPath& p) {
        if (use_short_loop_cov_resolver_) {
            return ResolveShortLoopByCov(p);
        } else {
            return ResolveShortLoopByPI(p);
        }
    }

    bool InvestigateShortLoop() {
        return investigate_short_loops_ && (use_short_loop_cov_resolver_ || CanInvestigateShortLoop());
    }
protected:
    DECL_LOGGER("LoopDetectingPathExtender")
};

class SimpleExtender: public LoopDetectingPathExtender {

protected:
    shared_ptr<ExtensionChooser> extensionChooser_;
    ShortLoopResolver loop_resolver_;
    double weight_threshold_;

    void FindFollowingEdges(BidirectionalPath& path, ExtensionChooser::EdgeContainer * result) {
        DEBUG("Looking for the following edges")
        result->clear();
        vector<EdgeId> edges;
        DEBUG("Pushing back")
        utils::push_back_all(edges, g_.OutgoingEdges(g_.EdgeEnd(path.Back())));
        result->reserve(edges.size());
        for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
            DEBUG("Adding edge w distance " << g_.int_id(*iter));
            result->push_back(EdgeWithDistance(*iter, 0));
        }
        DEBUG("Following edges found");
    }


public:

    SimpleExtender(const conj_graph_pack &gp,
                   const GraphCoverageMap &cov_map,
                   UsedUniqueStorage &unique,
                   shared_ptr<ExtensionChooser> ec,
                   size_t is,
                   bool investigate_short_loops,
                   bool use_short_loop_cov_resolver,
                   double weight_threshold = 0.0):
        LoopDetectingPathExtender(gp, cov_map, unique, investigate_short_loops, use_short_loop_cov_resolver, is),
        extensionChooser_(ec),
        loop_resolver_(gp.g, make_shared<CombinedLoopEstimator>(gp.g, gp.flanking_cov, extensionChooser_->wc(), weight_threshold)),
        weight_threshold_(weight_threshold) {}

    std::shared_ptr<ExtensionChooser> GetExtensionChooser() const {
        return extensionChooser_;
    }

    bool CanInvestigateShortLoop() const override {
        return extensionChooser_->WeightCounterBased();
    }

    bool ResolveShortLoopByPI(BidirectionalPath& path) override {
        if (extensionChooser_->WeightCounterBased()) {
            LoopDetector loop_detector(&path, cov_map_);
            size_t init_len = path.Length();
            bool result = false;
            while (path.Size() >= 1 && loop_detector.EdgeInShortLoop(path.Back())) {
                loop_resolver_.ResolveShortLoop(path);
                if (init_len == path.Length()) {
                    return result;
                } else {
                    result = true;
                }
                init_len = path.Length();
            }
            return true;
        }
        return false;
    }

    bool MakeSimpleGrowStep(BidirectionalPath& path, PathContainer* paths_storage) override {
        ExtensionChooser::EdgeContainer candidates;
        return FilterCandidates(path, candidates) && AddCandidates(path, paths_storage, candidates);
    }

protected:
    virtual bool FilterCandidates(BidirectionalPath& path, ExtensionChooser::EdgeContainer& candidates) {
        if (path.Size() == 0) {
            return false;
        }
        DEBUG("Simple grow step");
        path.PrintDEBUG();
        FindFollowingEdges(path, &candidates);
        DEBUG("found candidates");
        DEBUG(candidates.size())
        if (candidates.size() == 1) {
            LoopDetector loop_detector(&path, cov_map_);
            if (!investigate_short_loops_ && (loop_detector.EdgeInShortLoop(path.Back()) or loop_detector.EdgeInShortLoop(candidates.back().e_))
                && extensionChooser_->WeightCounterBased()) {
                return false;
            }
        }
        DEBUG("more filtering");
        candidates = extensionChooser_->Filter(path, candidates);
        DEBUG("filtered candidates");
        DEBUG(candidates.size())
        return true;
    }

    virtual bool AddCandidates(BidirectionalPath& path, PathContainer* /*paths_storage*/, ExtensionChooser::EdgeContainer& candidates) {
        if (candidates.size() != 1)
            return false;

        LoopDetector loop_detector(&path, cov_map_);
        DEBUG("loop detecor");
        if (!investigate_short_loops_ &&
            (loop_detector.EdgeInShortLoop(path.Back()) or loop_detector.EdgeInShortLoop(candidates.back().e_))
            && extensionChooser_->WeightCounterBased()) {
            return false;
        }
        DEBUG("push");
        EdgeId eid = candidates.back().e_;
//In 2015 modes when trying to use already used unique edge, it is not added and path growing stops.
//That allows us to avoid overlap removal hacks used earlier.
        Gap gap(candidates.back().d_);
        return TryUseEdge(path, eid, gap);
    }

    DECL_LOGGER("SimpleExtender")
};


class MultiExtender: public SimpleExtender {
public:
    using SimpleExtender::SimpleExtender;

protected:
    bool AddCandidates(BidirectionalPath& path, PathContainer* paths_storage, ExtensionChooser::EdgeContainer& candidates) override {
        if (candidates.size() == 0)
            return false;

        bool res = false;
        LoopDetector loop_detector(&path, cov_map_);
        DEBUG("loop detecor");
        if (!investigate_short_loops_ &&
            (loop_detector.EdgeInShortLoop(path.Back()) or loop_detector.EdgeInShortLoop(candidates.back().e_))
            && extensionChooser_->WeightCounterBased()) {
	    DEBUG("loop deteced");
            return false;
        }
        if (candidates.size() == 1) {
            DEBUG("push");
            EdgeId eid = candidates.back().e_;
            path.PushBack(eid, Gap(candidates.back().d_));
            DEBUG("push done");
            return true;
        }
        else if (candidates.size() == 2) {
             //Check for bulge
            auto v = g_.EdgeStart(candidates.front().e_);
            auto u = g_.EdgeEnd(candidates.front().e_);
            for (auto edge : candidates) {
                if (v != g_.EdgeStart(edge.e_) || u != g_.EdgeEnd(edge.e_))
                    return false;
            }

            //Creating new paths for other than new candidate.
            for (size_t i = 1; i < candidates.size(); ++i) {
                DEBUG("push other candidates " << i);
                BidirectionalPath *p = new BidirectionalPath(path);
                p->PushBack(candidates[i].e_, Gap(candidates[i].d_));
                BidirectionalPath *cp = new BidirectionalPath(p->Conjugate());
                paths_storage->AddPair(p, cp);
            }

            DEBUG("push");
            path.PushBack(candidates.front().e_, Gap(candidates.front().d_));
            DEBUG("push done");
            res = true;

            if (candidates.size() > 1) {
                DEBUG("Found " << candidates.size() << " candidates");
            }
        }

        return res;
    }

protected:
    DECL_LOGGER("MultiExtender")

};


class ScaffoldingPathExtender: public LoopDetectingPathExtender {
    std::shared_ptr<ExtensionChooser> extension_chooser_;
    ExtensionChooser::EdgeContainer sources_;
    std::shared_ptr<GapAnalyzer> gap_analyzer_;
    bool avoid_rc_connections_;

//When check_sink_ set to false we can scaffold not only tips
    bool check_sink_;

    void InitSources() {
        sources_.clear();

        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (g_.IncomingEdgeCount(g_.EdgeStart(*iter)) == 0) {
                sources_.push_back(EdgeWithDistance(*iter, 0));
            }
        }
    }

    bool IsSink(EdgeId e) const {
        return g_.OutgoingEdgeCount(g_.EdgeEnd(e)) == 0;
    }

    Gap ConvertGapDescription(const GapDescription &gap) const {
        if (gap == GapDescription()) {
            return Gap::INVALID();
        }
        return Gap(gap.estimated_dist() + int(g_.k())
                   - int(gap.left_trim()) - int(gap.right_trim()),
                   {uint32_t(gap.left_trim()), uint32_t(gap.right_trim())}, false);
    }

protected:
    virtual bool CheckGap(const Gap &/*gap*/) const {
        return true;
    }

    bool ResolveShortLoopByCov(BidirectionalPath&) override {
        return false;
    }

    bool ResolveShortLoopByPI(BidirectionalPath&) override {
        return false;
    }

    //TODO fix awful design with virtual CheckGap and must_overlap flag!
    bool MakeSimpleGrowStepForChooser(BidirectionalPath& path, std::shared_ptr<ExtensionChooser> ec,
                                      bool must_overlap = false) {
        if (path.Size() < 1 || (check_sink_ && !IsSink(path.Back()))) {
            return false;
        }

        DEBUG("Simple grow step, growing path");
        path.PrintDEBUG();
        ExtensionChooser::EdgeContainer candidates = ec->Filter(path, sources_);
        DEBUG("scaffolding candidates " << candidates.size() << " from sources " << sources_.size());

        DEBUG("Candidate size = " << candidates.size())
        if (candidates.size() != 1) {
            DEBUG("scaffolding end");
            return false;
        }

        EdgeId e = candidates.back().e_;
        if (e == path.Back()
            || (avoid_rc_connections_ && e == g_.conjugate(path.Back()))) {
            return false;
        }

        if (this->DetectCycleScaffolding(path, e)) {
            return false;
        }

        Gap gap;
        //TODO is it ok that we either force joining or ignore its possibility
        if (check_sink_) {
            gap = ConvertGapDescription(gap_analyzer_->FixGap(GapDescription(path.Back(), e,
                                                                             candidates.back().d_ -
                                                                             int(g_.k()))));

            if (gap == Gap::INVALID()) {
                DEBUG("Looks like wrong scaffolding. PathId: "
                              << path.GetId() << " path length: " << path.Length()
                              << ", estimated gap length: " << candidates.back().d_);
                return false;
            }

            DEBUG("Gap after fixing " << gap.gap << " (was " << candidates.back().d_ << ")");

            if (must_overlap && !CheckGap(gap)) {
                DEBUG("Overlap is not large enough")
                return false;
            }
        } else {
            DEBUG("Gap joiners off");
            VERIFY(candidates.back().d_ > int(g_.k()));
            gap = Gap(candidates.back().d_, false);
        }

        return TryUseEdge(path, e, NormalizeGap(gap));
    }

    Gap NormalizeGap(Gap gap) const {
        VERIFY(gap != Gap::INVALID());
        if (gap.overlap_after_trim(g_.k()) > 0)
            gap.trash.current += gap.overlap_after_trim(g_.k());
        return gap;
    }

public:

    ScaffoldingPathExtender(const conj_graph_pack &gp,
                            const GraphCoverageMap &cov_map,
                            UsedUniqueStorage &unique,
                            std::shared_ptr<ExtensionChooser> extension_chooser,
                            std::shared_ptr<GapAnalyzer> gap_analyzer,
                            size_t is,
                            bool investigate_short_loops,
                            bool avoid_rc_connections,
                            bool check_sink = true):
        LoopDetectingPathExtender(gp, cov_map, unique, investigate_short_loops, false, is),
        extension_chooser_(extension_chooser),
        gap_analyzer_(gap_analyzer),
        avoid_rc_connections_(avoid_rc_connections),
        check_sink_(check_sink)
    {
        InitSources();
    }

    bool MakeSimpleGrowStep(BidirectionalPath& path, PathContainer* /*paths_storage*/) override {
        return MakeSimpleGrowStepForChooser(path, extension_chooser_);
    }

    std::shared_ptr<ExtensionChooser> GetExtensionChooser() const {
        return extension_chooser_;
    }

protected:
    DECL_LOGGER("ScaffoldingPathExtender");
};


class RNAScaffoldingPathExtender: public ScaffoldingPathExtender {
    std::shared_ptr<ExtensionChooser> strict_extension_chooser_;

    int min_overlap_;

protected:
    bool CheckGap(const Gap &gap) const override {
        return gap.overlap_after_trim(g_.k()) >= min_overlap_;
    }

public:

    RNAScaffoldingPathExtender(const conj_graph_pack &gp,
                               const GraphCoverageMap &cov_map,
                               UsedUniqueStorage &unique,
                               std::shared_ptr<ExtensionChooser> extension_chooser,
                               std::shared_ptr<ExtensionChooser> strict_extension_chooser,
                               std::shared_ptr<GapAnalyzer> gap_joiner,
                               size_t is,
                               bool investigate_short_loops,
                               int min_overlap = 0):
        ScaffoldingPathExtender(gp, cov_map, unique, extension_chooser, gap_joiner, is, investigate_short_loops, true),
        strict_extension_chooser_(strict_extension_chooser), min_overlap_(min_overlap) {}


    bool MakeSimpleGrowStep(BidirectionalPath& path, PathContainer* /*paths_storage*/) override {
        DEBUG("==== RNA scaffolding ===");
        path.PrintDEBUG();
        bool res = MakeSimpleGrowStepForChooser(path, GetExtensionChooser(), true);
        if (!res) {
            DEBUG("==== Second strategy ====")
            res = MakeSimpleGrowStepForChooser(path, strict_extension_chooser_);
        }
        DEBUG("==== DONE ====")
        return res;
    }

};

}
