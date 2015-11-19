//***************************************************************************
///* Copyright (c) 2011-2014 Saint-Petersburg Academic University
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
#include "path_filter.hpp"
#include <cmath>

namespace path_extend {

class ShortLoopResolver {
public:
    ShortLoopResolver(const Graph& g)
            : g_(g) { }

    virtual ~ShortLoopResolver() { }

    virtual void ResolveShortLoop(BidirectionalPath& path) const = 0;

protected:
    DECL_LOGGER("PathExtender")
    const Graph& g_;

    void UndoCycles(BidirectionalPath& p, EdgeId next_edge) const {
        if (p.Size() <= 2) {
            return;
        }
        EdgeId first_edge = p.Back();
        EdgeId second_edge = next_edge;
        while (p.Size() > 2) {
            if (p.At(p.Size() - 1) == first_edge && p.At(p.Size() - 2) == second_edge) {
                p.PopBack(2);
            } else {
                return;;
            }
        }
    }

    void MakeCycleStep(BidirectionalPath& path, EdgeId e) const {
        if (path.Size() == 0) {
            return;
        }
        EdgeId pathEnd = path.Back();
        path.PushBack(e);
        path.PushBack(pathEnd);
    }
};

class CovShortLoopResolver : public ShortLoopResolver {
public:
    CovShortLoopResolver(const conj_graph_pack& gp)
            : ShortLoopResolver(gp.g), gp_(gp) {

    }

    void ResolveShortLoop(BidirectionalPath& path) const override {
        DEBUG("resolve short loop by coverage");
        path.Print();

        pair<EdgeId, EdgeId> edges;
        if (path.Size() >= 1 && GetLoopAndExit(g_, path.Back(), edges)) {
            DEBUG("Coverage Short Loop Resolver");
            UndoCycles(path, edges.first);
            EdgeId e1 = path.Back();
            EdgeId e2 = edges.first;
            EdgeId e_out = edges.second;
            auto prob_e_in = g_.IncomingEdges(g_.EdgeEnd(e2));
            EdgeId e_in = *prob_e_in.begin();
            size_t count = 0;
            for (auto edge = prob_e_in.begin(); edge != prob_e_in.end(); ++edge) {
                if (*edge != e2)
                    e_in = *edge;
                count++;
            }
            if (count != 2) {
                return;
            }
            double in_cov = gp_.flanking_cov.GetOutCov(e_in); //g_.coverage(e_in);
            double out_cov = gp_.flanking_cov.GetInCov(e_out); //g_.coverage(e_out);
            double cov = (in_cov + out_cov) / 2.0;
            double time1 = math::round(gp_.flanking_cov.GetInCov(e1) / cov);//math::round(gp_.g.coverage(e1) / cov);
            double time2 = math::round(gp_.flanking_cov.GetInCov(e2) / cov);////math::round(gp_.g.coverage(e2) / cov);
            size_t time = (size_t) std::max(0.0, std::min(time1 - 1.0, time2));
            for (size_t i = 0; i < time; ++i) {
                MakeCycleStep(path, edges.first);
            }
            path.PushBack(edges.second);
            DEBUG("loop with start " << g_.int_id(e_in)
            		<<" e1 " << g_.int_id(e1)
            		<< " e2 " << g_.int_id(e2)
            		<< " out " <<g_.int_id(e_out)
            		<< " cov in = " << in_cov
            		<< " cov out " << out_cov
            		<< " cov " << cov
                  << " cov e1 = " << gp_.g.coverage(e1)
                  << " cov e2 = " << gp_.g.coverage(e2)
                  << " time1 = " << time1
                  << " time2 = " << time2
                  << " time = " << time);
        }
    }
private:
    const conj_graph_pack& gp_;
};

class SimpleLoopResolver : public ShortLoopResolver {

public:
    SimpleLoopResolver(Graph& g) : ShortLoopResolver(g) { }

    void ResolveShortLoop(BidirectionalPath& path) const override {
        pair<EdgeId, EdgeId> edges;
        if (path.Size() >= 1 && GetLoopAndExit(g_, path.Back(), edges)) {
            DEBUG("Resolving short loop...");
            EdgeId e = path.Back();
            path.PushBack(edges.first);
            path.PushBack(e);
            path.PushBack(edges.second);
            DEBUG("Resolving short loop done");
        }
    }

protected:
    DECL_LOGGER("PathExtender")
};

class LoopResolver : public ShortLoopResolver {
    static const size_t ITER_COUNT = 10;
    const WeightCounter& wc_;

public:
    LoopResolver(const Graph& g, const WeightCounter& wc)
            : ShortLoopResolver(g),
              wc_(wc) { }

    void MakeBestChoice(BidirectionalPath& path, pair<EdgeId, EdgeId>& edges) const {
        UndoCycles(path, edges.first);
        BidirectionalPath experiment(path);
        double max_weight = wc_.CountWeight(experiment, edges.second);
        double diff = max_weight - wc_.CountWeight(experiment, edges.first);
        size_t maxIter = 0;
        for (size_t i = 1; i <= ITER_COUNT; ++i) {
            double weight = wc_.CountWeight(experiment, edges.first);
            if (weight > 0) {
                MakeCycleStep(experiment, edges.first);
                weight = wc_.CountWeight(experiment, edges.second);
                double weight2 = wc_.CountWeight(experiment, edges.first);
                if (weight > max_weight || (weight == max_weight && weight - weight2 > diff)
                        || (weight == max_weight && weight - weight2 == diff && i == 1)) {
                    max_weight = weight;
                    maxIter = i;
                    diff = weight - weight2;
                }
            }
        }
        for (size_t i = 0; i < maxIter; ++i) {
            MakeCycleStep(path, edges.first);
        }
        path.PushBack(edges.second);
    }

    void ResolveShortLoop(BidirectionalPath& path) const override {
        pair<EdgeId, EdgeId> edges;
        if (path.Size() >=1 && GetLoopAndExit(g_, path.Back(), edges)) {
            DEBUG("Resolving short loop...");
            MakeBestChoice(path, edges);
            DEBUG("Resolving short loop done");
        }
    }
};


class GapJoiner {

public:
    static const int INVALID_GAP = -1000000;
    GapJoiner(const Graph& g)
            : g_(g) { }

    virtual int FixGap(EdgeId sink, EdgeId source, int initial_gap) const = 0;

    virtual ~GapJoiner() { }
protected:
    const Graph& g_;
    DECL_LOGGER("PathExtender")

};

class SimpleGapJoiner : public GapJoiner {

public:
    SimpleGapJoiner(const Graph& g) : GapJoiner(g) { }

    virtual int FixGap(EdgeId sink, EdgeId source, int initial_gap) const {
        if (initial_gap > 2 * (int) g_.k()) {
            return initial_gap;
        }
        for (int l = (int) g_.k(); l > 0; --l) {
            if (g_.EdgeNucls(sink).Subseq(g_.length(sink) + g_.k() - l) == g_.EdgeNucls(source).Subseq(0, l)) {
                DEBUG("Found correct gap length");
                DEBUG("Inintial: " << initial_gap << ", new gap: " << g_.k() - l);
                return (int) g_.k() - l;
            }
        }
        DEBUG("Perfect overlap is not found, inintial: " << initial_gap);
        return initial_gap;
    }
};

class HammingGapJoiner: public GapJoiner {
    static const size_t DEFAULT_PADDING_LENGTH = 10;
    const double min_gap_score_;
    const int must_overlap_threshold_;
    const size_t may_overlap_threshold_;
    const size_t short_overlap_threshold_;
    const size_t basic_overlap_length_;
    const size_t artificial_gap_;
    const bool use_old_score_;

    vector<size_t> DiffPos(const Sequence& s1, const Sequence& s2) const {
        VERIFY(s1.size() == s2.size());
        vector < size_t > answer;
        for (size_t i = 0; i < s1.size(); ++i)
            if (s1[i] != s2[i])
                answer.push_back(i);
        return answer;
    }

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

//    double ScoreGap(const Sequence& s1, const Sequence& s2, int gap, int initial_gap) const {
//        VERIFY(s1.size() == s2.size());
//        return 1.0 - (double) HammingDistance(s1, s2) / (double) s1.size()
//                - (double) abs(gap - initial_gap) / (double) (2 * g_.k());
//    }

    //FIXME use GC content, change match prob and use partition of tip sequence into bad and good part
    double ScoreGap(const Sequence& s1, const Sequence& s2) const {
        static double match_prob = 0.9;
        static double log_match_prob = log2(match_prob);
        static double log_mismatch_prob = log2(1. - match_prob);
        VERIFY(s1.size() == s2.size());
        size_t n = s1.size();
        size_t mismatches = HammingDistance(s1, s2);
        VERIFY(mismatches <= n);
        return 2.*double(n) + double(n - mismatches) * log_match_prob + double(mismatches) * log_mismatch_prob;
    }

    double OldScoreGap(const Sequence& s1, const Sequence& s2) const {
        VERIFY(s1.size() == s2.size());
        return 1.0 - (double) HammingDistance(s1, s2) / (double) s1.size();
    }

public:

    //todo review parameters in usages
    HammingGapJoiner(const Graph& g,
            double min_gap_score,
            int must_overlap_threshold,
            size_t may_overlap_threshold,
            size_t short_overlap_threshold,
            size_t basic_overlap_length,
            size_t artificial_gap = DEFAULT_PADDING_LENGTH,
            bool use_old_score = false):
                GapJoiner(g),
                min_gap_score_(min_gap_score),
                must_overlap_threshold_(must_overlap_threshold),
                may_overlap_threshold_(may_overlap_threshold),
                short_overlap_threshold_(short_overlap_threshold),
                basic_overlap_length_(basic_overlap_length),
                artificial_gap_(artificial_gap),
                use_old_score_(use_old_score)
    {
        DEBUG("HammingGapJoiner params: \n min_gap_score " << min_gap_score_ <<
              "\n must_overlap_threshold " << must_overlap_threshold_ <<
              "\n may_overlap_threshold " << may_overlap_threshold_ <<
              "\n short_overlap_threshold " << short_overlap_threshold_ <<
              "\n basic_overlap_length " << basic_overlap_length_ <<
              "\n artificial_gap " << artificial_gap_);
    }

    //estimated_gap is in k-mers
    virtual int FixGap(EdgeId sink, EdgeId source, int estimated_gap) const {
        DEBUG("Trying to fix estimated gap " << estimated_gap <<
              " between " << g_.str(sink) << " and " << g_.str(source));

        if (estimated_gap > int(g_.k() + may_overlap_threshold_)) {
            DEBUG("Edges are supposed to be too far to check overlaps");
            return estimated_gap;
        }

        size_t corrected_start_overlap = basic_overlap_length_;
        if (estimated_gap < 0) {
            corrected_start_overlap -= estimated_gap;
        }

        corrected_start_overlap = min(corrected_start_overlap,
                                      g_.k() + min(g_.length(sink), g_.length(source)));

        DEBUG("Corrected max overlap " << corrected_start_overlap);

        double best_score = min_gap_score_;
        int fixed_gap = INVALID_GAP;

        double overlap_coeff = 0.3;
        size_t min_overlap = 1ul;
        if (estimated_gap < 0) {
            size_t estimated_overlap = g_.k() - estimated_gap;
            min_overlap = max(size_t(math::round(overlap_coeff * double(estimated_overlap))), 1ul);
        }
        //todo better usage of estimated overlap
        DEBUG("Min overlap " << min_overlap);

        for (size_t l = corrected_start_overlap; l >= min_overlap; --l) {
            //TRACE("Sink: " << g_.EdgeNucls(sink).Subseq(g_.length(sink) + g_.k() - l).str());
            //TRACE("Source: " << g_.EdgeNucls(source).Subseq(0, l));
            double score = 0;
            if(use_old_score_)
                score = OldScoreGap(g_.EdgeNucls(sink).Subseq(g_.length(sink) + g_.k() - l),
                                    g_.EdgeNucls(source).Subseq(0, l));
            else
                score = ScoreGap(g_.EdgeNucls(sink).Subseq(g_.length(sink) + g_.k() - l),
                                    g_.EdgeNucls(source).Subseq(0, l));
            if (math::gr(score, best_score)) {
                TRACE("Curr overlap " << l);
                TRACE("Score: " << score);
                best_score = score;
                fixed_gap = int(g_.k() - l);
            }

            if (l == short_overlap_threshold_ && fixed_gap != INVALID_GAP) {
                //look at "short" overlaps only if long overlaps couldn't be found
                DEBUG("Not looking at short overlaps");
                break;
            }
        }

        if (fixed_gap != INVALID_GAP) {
            DEBUG("Found candidate gap length with score " << best_score);
            DEBUG("Estimated gap: " << estimated_gap <<
                  ", fixed gap: " << fixed_gap << " (overlap " << g_.k() - fixed_gap<< ")");
        } else {
            //couldn't find decent overlap
            if (estimated_gap < must_overlap_threshold_) {
                DEBUG("Estimated gap looks unreliable");
            } else {
                DEBUG("Overlap was not found");
                fixed_gap = max(estimated_gap, int(g_.k() + artificial_gap_));
            }
        }
        DEBUG("Final fixed gap " << fixed_gap);
        return fixed_gap;
    }

private:
    DECL_LOGGER("HammingGapJoiner");
};

//Detects a cycle as a minsuffix > IS present earlier in the path. Overlap is allowed.
class InsertSizeLoopDetector {
protected:
    const Graph& g_;
    const GraphCoverageMap& cov_map_;
    size_t min_cycle_len_;

public:
    InsertSizeLoopDetector(const Graph& g, const GraphCoverageMap& cov_map, size_t is): g_(g), cov_map_(cov_map), min_cycle_len_(is) {
    }

    size_t GetMinCycleLenth() const {
        return min_cycle_len_;
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
};

class RepeatDetector {
public:
    RepeatDetector(const Graph& g, const GraphCoverageMap& cov_map, size_t max_repeat_len)
            : g_(g),
              cov_map_(cov_map),
              used_paths_(),
              repeat_len_(max_repeat_len){
        empty_ = new BidirectionalPath(g_);
    }
    ~RepeatDetector() {
        delete empty_;
    }

    BidirectionalPath* RepeatPath(const BidirectionalPath& p) {
        if (p.Size() == 0) {
            return empty_;
        }
        EdgeId last_e = p.Back();
        BidirectionalPathSet cov_paths = cov_map_.GetCoveringPaths(last_e);
        DEBUG("cov paths for e " << g_.int_id(last_e) << " size " << cov_paths.size());
        size_t max_common_size = 0;
        BidirectionalPath* result_p = empty_;
        for (BidirectionalPath* cov_p : cov_paths) {
            if (used_paths_.find(cov_p) == used_paths_.end() || cov_p == &p || cov_p == p.GetConjPath()) {
                continue;
            }
            size_t common_size = MaxCommonSize(p, *cov_p);
            DEBUG("max comon size with path " << cov_p->GetId() << " is " << common_size);
            if (common_size == 0) {
            	continue;
            }
            VERIFY(common_size <= p.Size());
            if (p.LengthAt(p.Size() - common_size) > repeat_len_) {
                DEBUG("repeat from " << (p.Size() - common_size) << " length " << p.LengthAt(p.Size() - common_size) << " repeat length " << repeat_len_);
                max_common_size = max(common_size, max_common_size);
                result_p = cov_p;
            }
        }
        used_paths_.insert(&p);
        DEBUG("max common size " << max_common_size);
        return result_p;
    }
    size_t MaxCommonSize(const BidirectionalPath& p1, const BidirectionalPath& p2) const {
        DEBUG("max coomon size ")
    	EdgeId last_e = p1.Back();
        vector<size_t> positions2 = p2.FindAll(last_e);
        DEBUG("pos size " << positions2.size())
        size_t max_common_size = 0;
        for (size_t pos2 : positions2) {
            size_t common_size = MaxCommonSize(p1, p1.Size() - 1, p2, pos2);
            DEBUG("max common size from " << pos2 << " is " << common_size);
            max_common_size = max(max_common_size, common_size);
        }
        return max_common_size;
    }
private:
    size_t MaxCommonSize(const BidirectionalPath& p1, size_t pos1, const BidirectionalPath& p2, size_t pos2) const {
        int i1 = (int) pos1;
        int i2 = (int) pos2;
        while (i1 >= 0 && i2 >= 0 &&
        		p1.At((size_t) i1) == p2.At((size_t) i2) &&
        		p1.GapAt((size_t) i1) == p2.GapAt((size_t) i2)) {
            i1--;
            i2--;
        }
        if (i1 >=0 && i2>=0 && p1.At((size_t) i1) == p2.At((size_t) i2)) {
        	i1--;
        	i2--;
        }

        VERIFY(i1 <= (int)pos1);
        return std::max(size_t((int) pos1 - i1), (size_t)1);
    }
    const Graph& g_;
    const GraphCoverageMap& cov_map_;
    set<const BidirectionalPath*> used_paths_;
    size_t repeat_len_;
    BidirectionalPath* empty_;
};

class ContigsMaker {
public:
    ContigsMaker(const Graph & g)
            : g_(g) { }

    virtual ~ContigsMaker() { }

    virtual void GrowPath(BidirectionalPath& path) = 0;

    virtual void GrowPathSimple(BidirectionalPath& path) = 0;

    virtual void GrowAll(PathContainer & paths, PathContainer * result) = 0;

protected:
    const Graph& g_;
    DECL_LOGGER("PathExtender")
};

struct UsedUniqueStorage {
    set<EdgeId> used_;

    shared_ptr<ScaffoldingUniqueEdgeStorage> unique_;
    void insert(EdgeId e) {
        if (unique_->IsUnique(e)) {
            used_.insert(e);
            used_.insert(e->conjugate());
        }
    }
    bool IsUsedAndUnique (EdgeId e) {
        return (unique_->IsUnique(e) && used_.find(e) != used_.end());
    }
    UsedUniqueStorage(  shared_ptr<ScaffoldingUniqueEdgeStorage> unique ):used_(), unique_(unique) {}
};
class PathExtender {
public:
    PathExtender(const Graph & g): g_(g){ }
    virtual ~PathExtender() { }
    virtual bool MakeGrowStep(BidirectionalPath& path) = 0;
    void AddUniqueEdgeStorage(shared_ptr<UsedUniqueStorage> used_storage) {
        used_storage_ = used_storage;
    }
protected:
    const Graph& g_;
    shared_ptr<UsedUniqueStorage> used_storage_;
    DECL_LOGGER("PathExtender")
};

class CompositeExtender : public ContigsMaker {
public:
    CompositeExtender(Graph & g, GraphCoverageMap& cov_map, size_t max_diff_len)
            : ContigsMaker(g),
              cover_map_(cov_map),
              repeat_detector_(g, cover_map_, 2 * cfg::get().max_repeat_length),  //TODO: move to config
              extenders_(),
              max_diff_len_(max_diff_len) {
    }

    CompositeExtender(Graph & g, GraphCoverageMap& cov_map, vector<shared_ptr<PathExtender> > pes, size_t max_diff_len, shared_ptr<ScaffoldingUniqueEdgeStorage> unique)
            : ContigsMaker(g),
              cover_map_(cov_map),
              repeat_detector_(g, cover_map_, 2 * cfg::get().max_repeat_length),  //TODO: move to config
              extenders_(),
              max_diff_len_(max_diff_len) {
        extenders_ = pes;
        used_storage_ = make_shared<UsedUniqueStorage>(UsedUniqueStorage( unique));
        for (auto ex: extenders_) {
            ex->AddUniqueEdgeStorage(used_storage_);
        }
    }

    void AddExtender(shared_ptr<PathExtender> pe) {
        extenders_.push_back(pe);
        pe->AddUniqueEdgeStorage(used_storage_);
    }

    virtual void GrowAll(PathContainer& paths, PathContainer * result) {
        result->clear();
        PathContainer usedPaths;
        GrowAll(paths, usedPaths, result);
        LengthPathFilter filter(g_, 0);
        filter.filter(*result);
    }

    virtual void GrowPath(BidirectionalPath& path) {
        while (MakeGrowStep(path)) { }
    }

    virtual void GrowPathSimple(BidirectionalPath& path) {
        while (MakeGrowStep(path, false)) { }
    }

    bool MakeGrowStep(BidirectionalPath& path, bool detect_repeats_online = true) {
        DEBUG("make grow step composite extender");
        auto sc_mode = cfg::get().pe_params.param_set.sm;
        if (sc_mode == sm_old_pe_2015 || sc_mode == sm_2015 || sc_mode == sm_combined) {
            DEBUG("force switch off online repeats detect, 2015 on");
            detect_repeats_online = false;
        }
        if (detect_repeats_online) {
            BidirectionalPath *repeat_path = repeat_detector_.RepeatPath(path);
            size_t repeat_size = repeat_detector_.MaxCommonSize(path, *repeat_path);

            if (repeat_size > 0) {
                DEBUG("repeat with length " << repeat_size);
                path.Print();
                repeat_path->Print();
                BidirectionalPath repeat = path.SubPath(path.Size() - repeat_size);
                int begin_repeat = repeat_path->FindLast(repeat);
                VERIFY(begin_repeat > -1);
                size_t end_repeat = (size_t) begin_repeat + repeat_size;
                DEBUG("not consistent subpaths ");
                BidirectionalPath begin1 = path.SubPath(0, path.Size() - repeat_size);
                begin1.Print();
                BidirectionalPath begin2 = repeat_path->SubPath(0, begin_repeat);
                begin2.Print();
                int gpa_in_repeat_path = repeat_path->GapAt(begin_repeat);
                BidirectionalPath end2 = repeat_path->SubPath(end_repeat);
                BidirectionalPath begin1_conj = path.SubPath(0, path.Size() - repeat_size + 1).Conjugate();
                BidirectionalPath begin2_conj = repeat_path->SubPath(0, begin_repeat + 1).Conjugate();
                pair<size_t, size_t> last = ComparePaths(0, 0, begin1_conj, begin2_conj, max_diff_len_);
                DEBUG("last " << last.first << " last2 " << last.second);
                path.Clear();
                repeat_path->Clear();
                int gap_len = repeat.GapAt(0);

                if (begin2.Size() == 0 || last.second != 0) { //TODO: incorrect: common edges, but then different ends
                    path.PushBack(begin1);
                    repeat_path->PushBack(begin2);
                } else {
                    gap_len = gpa_in_repeat_path;
                    path.PushBack(begin2);
                    repeat_path->PushBack(begin1);
                }

                path.PushBack(repeat.At(0), gap_len);
                path.PushBack(repeat.SubPath(1));
                path.PushBack(end2);
                DEBUG("new path");
                path.Print();
                return false;
            }
        }

        size_t current = 0;
        while (current < extenders_.size()) {
            DEBUG("step " << current << " from " <<extenders_.size());
            if (extenders_[current]->MakeGrowStep(path)) {
                return true;
            }
           ++current;
        }
        return false;
    }
    
private:
    GraphCoverageMap& cover_map_;
    RepeatDetector repeat_detector_;
    vector<shared_ptr<PathExtender> > extenders_;
    size_t max_diff_len_;
    shared_ptr<UsedUniqueStorage> used_storage_;
    void SubscribeCoverageMap(BidirectionalPath * path) {
        path->Subscribe(&cover_map_);
        for (size_t i = 0; i < path->Size(); ++i) {
            cover_map_.BackEdgeAdded(path->At(i), path, path->GapAt(i));
        }
    }

    void GrowAll(PathContainer& paths, PathContainer& usedPaths, PathContainer * result) {
        cover_map_.Clear();
        for (size_t i = 0; i < paths.size(); ++i) {
            VERBOSE_POWER_T2(i, 100, "Processed " << i << " paths from " << paths.size() << " (" << i * 100 / paths.size() << "%)");
            if (paths.size() > 10 && i % (paths.size() / 10 + 1) == 0) {
                INFO("Processed " << i << " paths from " << paths.size() << " (" << i * 100 / paths.size() << "%)");
            }
//In 2015 modes do not use a seed already used in paths.
            auto sc_mode = cfg::get().pe_params.param_set.sm;
            if (sc_mode == sm_old_pe_2015 || sc_mode == sm_2015 || sc_mode == sm_combined) {
                bool was_used = false;
                for (size_t ind =0; ind < paths.Get(i)->Size(); ind++) {
                    EdgeId eid = paths.Get(i)->At(ind);
                    if (used_storage_->IsUsedAndUnique(eid)) {
                        was_used = true; break;
                    } else {
                        used_storage_->insert(eid);
                    }
                }
                if (was_used) {
                    DEBUG("skipping already used seed");
                    continue;
                }
            }
//TODO: coverage_map should be exterminated
            if (!cover_map_.IsCovered(*paths.Get(i))) {
                usedPaths.AddPair(paths.Get(i), paths.GetConjugate(i));
                BidirectionalPath * path = new BidirectionalPath(*paths.Get(i));
                BidirectionalPath * conjugatePath = new BidirectionalPath(*paths.GetConjugate(i));
                result->AddPair(path, conjugatePath);
                SubscribeCoverageMap(path);
                SubscribeCoverageMap(conjugatePath);
                size_t count_trying = 0;
                size_t current_path_len = 0;
                do {
                    current_path_len = path->Length();
                    count_trying++;
                    GrowPath(*path);
                    GrowPath(*conjugatePath);
                } while (count_trying < 10 && (path->Length() != current_path_len));
                path->CheckConjugateEnd(cfg::get().max_repeat_length);
                DEBUG("result path " << path->GetId());
                path->Print();
            }
        }
    }

};

//All Path-Extenders inherits this one.

class LoopDetectingPathExtender : public PathExtender {

protected:
    size_t maxLoops_;
    bool investigateShortLoops_;
    bool use_short_loop_cov_resolver_;
    CovShortLoopResolver cov_loop_resolver_;

    vector<shared_ptr<BidirectionalPath> > visited_cycles_;
    InsertSizeLoopDetector is_detector_;
    const GraphCoverageMap& cov_map_;

public:
    LoopDetectingPathExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map, size_t max_loops, bool investigateShortLoops,
                              bool use_short_loop_cov_resolver, size_t is)
            : PathExtender(gp.g),
              maxLoops_(max_loops),
              investigateShortLoops_(investigateShortLoops),
              use_short_loop_cov_resolver_(use_short_loop_cov_resolver),
              cov_loop_resolver_(gp),
              is_detector_(gp.g, cov_map, is),
              cov_map_(cov_map) {

    }

    size_t getMaxLoops() const {
        return maxLoops_;
    }

    bool isInvestigateShortLoops() const {
        return investigateShortLoops_;
    }

    void setInvestigateShortLoops(bool investigateShortLoops) {
        this->investigateShortLoops_ = investigateShortLoops;
    }

    void setMaxLoops(size_t maxLoops) {
        if (maxLoops != 0) {
            this->maxLoops_ = maxLoops;
        }
    }
//seems that it is outofdate
    bool InExistingLoop(const BidirectionalPath& path) {
        TRACE("Checking existing loops");
        int j = 0;
        for (auto cycle : visited_cycles_) {
            VERBOSE_POWER2(j++, "checking ");
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
                path.Print();
                DEBUG("last subpath");
                path.SubPath(last_cycle_pos).Print();
                DEBUG("cycle");
                cycle->Print();
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
        visited_cycles_.push_back(std::make_shared<BidirectionalPath>(path.SubPath(pos)));
        DEBUG("add cycle");
        path.SubPath(pos).Print();
    }

    bool DetectCycle(BidirectionalPath& path) {
        DEBUG("detect cycle");
        if (is_detector_.CheckCycled(path)) {
            DEBUG("Checking IS cycle");
            int loop_pos = is_detector_.RemoveCycle(path);
            DEBUG("Removed IS cycle");
            if (loop_pos != -1) {
                AddCycledEdges(path, loop_pos);
                return true;
            }
        }
        return false;
    }

    virtual bool MakeSimpleGrowStep(BidirectionalPath& path) = 0;

    virtual bool ResolveShortLoopByCov(BidirectionalPath& path) = 0;

    virtual bool ResolveShortLoopByPI(BidirectionalPath& path) = 0;

    virtual bool CanInvestigateShortLoop() const {
        return false;
    }

    virtual bool MakeGrowStep(BidirectionalPath& path) {
        if (InExistingLoop(path)) {
            DEBUG("in existing loop");
            return false;
        }
        bool result = false;
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
            result = MakeSimpleGrowStep(path);
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
        return investigateShortLoops_ && (use_short_loop_cov_resolver_ || CanInvestigateShortLoop());
    }
protected:
    DECL_LOGGER("LoopDetectingPathExtender")
};

class SimpleExtender: public LoopDetectingPathExtender {

protected:

    shared_ptr<ExtensionChooser> extensionChooser_;

    void FindFollowingEdges(BidirectionalPath& path, ExtensionChooser::EdgeContainer * result) {
        DEBUG("Looking for the following edges")
        result->clear();
        vector<EdgeId> edges;
        DEBUG("Pushing back")
        push_back_all(edges, g_.OutgoingEdges(g_.EdgeEnd(path.Back())));
        result->reserve(edges.size());
        for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
            DEBUG("Adding edge w distance " << g_.int_id(*iter));
            result->push_back(EdgeWithDistance(*iter, 0));
        }
        DEBUG("Following edges found");
    }


public:

    SimpleExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map, shared_ptr<ExtensionChooser> ec, 
                    size_t is, size_t max_loops, bool investigate_short_loops, bool use_short_loop_cov_resolver):
        LoopDetectingPathExtender(gp, cov_map, max_loops, investigate_short_loops, use_short_loop_cov_resolver, is),
        extensionChooser_(ec) {
    }

    std::shared_ptr<ExtensionChooser> GetExtensionChooser() const {
        return extensionChooser_;
    }

    virtual bool MakeSimpleGrowStep(BidirectionalPath& path) {
        if (path.Size() == 0) {
            return false;
        }
        DEBUG("Simple grow step");
        path.Print();
        ExtensionChooser::EdgeContainer candidates;
        FindFollowingEdges(path, &candidates);
        DEBUG("found candidates");
        DEBUG(candidates.size())
        if (candidates.size() == 1) {
            LoopDetector loop_detector(&path, cov_map_);
            if (!investigateShortLoops_ && (loop_detector.EdgeInShortLoop(path.Back()) or loop_detector.EdgeInShortLoop(candidates.back().e_))
                    && extensionChooser_->WeightCounterBased()) {
                return false;
            }
        }
        DEBUG("more filtering");
        candidates = extensionChooser_->Filter(path, candidates);
        DEBUG("found candidates 2");
        DEBUG(candidates.size())
        if (candidates.size() == 1) {
            LoopDetector loop_detector(&path, cov_map_);
            DEBUG("loop detecor");
            if (!investigateShortLoops_ &&
                    (loop_detector.EdgeInShortLoop(path.Back())  or loop_detector.EdgeInShortLoop(candidates.back().e_))
                    && extensionChooser_->WeightCounterBased()) {
                return false;
            }
            DEBUG("push");
            auto sc_mode = cfg::get().pe_params.param_set.sm;
            EdgeId eid = candidates.back().e_;
//In 2015 modes when trying to use already used unique edge, it is not added and path growing stops.
//That allows us to avoid overlap removal hacks used earlier.
            if (sc_mode == sm_old_pe_2015 || sc_mode == sm_2015 || sc_mode == sm_combined) {
                if (used_storage_->IsUsedAndUnique(eid)) {
                    return false;
                } else {
                    used_storage_->insert(eid);
                }
            }
            path.PushBack(eid, candidates.back().d_);
            DEBUG("push done");
            return true;
        }
        return false;
    }


    bool CanInvestigateShortLoop() const override {
        return extensionChooser_->WeightCounterBased();
    }

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

    virtual bool ResolveShortLoopByPI(BidirectionalPath& path) {
        if (extensionChooser_->WeightCounterBased()) {
            LoopResolver loop_resolver(g_, extensionChooser_->wc());
            LoopDetector loop_detector(&path, cov_map_);
            size_t init_len = path.Length();
            bool result = false;
            while (path.Size() >= 1 && loop_detector.EdgeInShortLoop(path.Back())) {
                loop_resolver.ResolveShortLoop(path);
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

protected:
    DECL_LOGGER("SimpleExtender")

};

class ScaffoldingPathExtender: public LoopDetectingPathExtender {
    std::shared_ptr<ExtensionChooser> extension_chooser_;
    ExtensionChooser::EdgeContainer sources_;
    std::shared_ptr<GapJoiner> gap_joiner_;
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

    bool IsSink(EdgeId e) const	{
		return g_.OutgoingEdgeCount(g_.EdgeEnd(e)) == 0;
	}


public:

    ScaffoldingPathExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map, std::shared_ptr<ExtensionChooser> extension_chooser,
                            std::shared_ptr<GapJoiner> gap_joiner, size_t is, size_t max_loops, bool investigateShortLoops, bool check_sink = true):
        LoopDetectingPathExtender(gp, cov_map, max_loops, investigateShortLoops, false, is),
            extension_chooser_(extension_chooser),
            gap_joiner_(gap_joiner),check_sink_(check_sink)
    {
        InitSources();
    }

    virtual bool MakeSimpleGrowStep(BidirectionalPath& path) {
        if (path.Size() < 1 || (check_sink_ && !IsSink(path.Back())) ) {
            return false;
        }
        DEBUG("scaffolding:");
        DEBUG("Simple grow step, growing path");
        path.Print();
        ExtensionChooser::EdgeContainer candidates = extension_chooser_->Filter(path, sources_);
        DEBUG("scaffolding candidates " << candidates.size() << " from sources " << sources_.size());
        if (candidates.size() == 1) {
            if (candidates[0].e_ == path.Back() || (cfg::get().avoid_rc_connections && candidates[0].e_ == g_.conjugate(path.Back()))) {
                return false;
            }
            int gap = cfg::get().pe_params.param_set.scaffolder_options.fix_gaps ?
                            gap_joiner_->FixGap(path.Back(), candidates.back().e_, candidates.back().d_) : candidates.back().d_;

            if (gap != GapJoiner::INVALID_GAP) {
                DEBUG("Scaffolding. PathId: " << path.GetId() << " path length: " << path.Length() << ", fixed gap length: " << gap);
                auto sc_mode = cfg::get().pe_params.param_set.sm;
                EdgeId eid = candidates.back().e_;
                if (sc_mode == sm_old_pe_2015 || sc_mode == sm_2015 || sc_mode == sm_combined) {
                    if (used_storage_->IsUsedAndUnique(eid)) {
                        return false;
                    } else {
                        used_storage_->insert(eid);
                    }
                }
                path.PushBack(eid, gap);
                return true;


            } else {
                DEBUG("Looks like wrong scaffolding. PathId: " << path.GetId() << " path length: " << path.Length() << ", fixed gap length: " << candidates.back().d_);
                return false;
            }
        }
        DEBUG("scaffolding end");
        return false;
    }

    virtual bool ResolveShortLoopByCov(BidirectionalPath&) {
        return false;
    }

	virtual bool ResolveShortLoopByPI(BidirectionalPath&) {
		return false;
	}

    std::shared_ptr<ExtensionChooser> GetExtensionChooser() const {
        return extension_chooser_;
    }

private:
	DECL_LOGGER("ScaffoldingPathExtender");
};

}
