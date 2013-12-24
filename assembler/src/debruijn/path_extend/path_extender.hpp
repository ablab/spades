//***************************************************************************
///* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * path_extender.hpp
 *
 *  Created on: Mar 5, 2012
 *      Author: andrey
 */

#ifndef PATH_EXTENDER_HPP_
#define PATH_EXTENDER_HPP_

#include "pe_utils.hpp"
#include "extension_chooser.hpp"
#include "path_filter.hpp"

namespace path_extend {

class ShortLoopResolver {
public:
    ShortLoopResolver(const Graph& g)
            : g_(g) { }

    virtual ~ShortLoopResolver() { }

    virtual void ResolveShortLoop(BidirectionalPath& path) = 0;

protected:
    DECL_LOGGER("PathExtender")
    const Graph& g_;

    bool GetLoopAndExit(BidirectionalPath& path, pair<EdgeId, EdgeId>& result) const {
        EdgeId e = path.Head();
        VertexId v = g_.EdgeEnd(e);
        if (g_.OutgoingEdgeCount(v) != 2) {
            return false;
        }
        EdgeId loop;
        EdgeId exit;
        bool loop_found = false;
        bool exit_found = false;
        auto edges = g_.OutgoingEdges(v);
        for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
            if (g_.EdgeEnd(*edge) == g_.EdgeStart(e)) {
                loop = *edge;
                loop_found = true;
            } else {
                exit = *edge;
                exit_found = true;
            }
        }
        result = make_pair(loop, exit);
        return exit_found && loop_found;
    }
};

class SimpleLoopResolver : public ShortLoopResolver {

public:
    SimpleLoopResolver(Graph& g) : ShortLoopResolver(g) { }

    virtual void ResolveShortLoop(BidirectionalPath& path) {
        pair<EdgeId, EdgeId> edges;
        if (GetLoopAndExit(path, edges)) {
            DEBUG("Resolving short loop...");
            EdgeId e = path.Head();
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
    static const size_t iter_ = 10;
    ExtensionChooser& chooser_;

public:
    LoopResolver(const Graph& g, ExtensionChooser& chooser)
            : ShortLoopResolver(g),
              chooser_(chooser) { }

    void MakeCycleStep(BidirectionalPath& path, EdgeId e) {
        EdgeId pathEnd = path.Head();
        path.PushBack(e);
        path.PushBack(pathEnd);
    }

    void MakeBestChoice(BidirectionalPath& path, pair<EdgeId, EdgeId>& edges) {
        EdgeId first_edge = path.Back();
        EdgeId second_edge = edges.first;
        while (path.Size() > 2) {
            if (path.At(path.Size() - 1) == first_edge && path.At(path.Size() - 2) == second_edge) {
                path.PopBack(2);
            } else {
                break;
            }
        }
        chooser_.ClearExcludedEdges();
        BidirectionalPath experiment(path);
        double maxWeight = chooser_.CountWeight(experiment, edges.second);
        double diff = maxWeight - chooser_.CountWeight(experiment, edges.first);
        size_t maxIter = 0;
        for (size_t i = 1; i <= iter_; ++i) {
            double weight = chooser_.CountWeight(experiment, edges.first);
            if (weight > 0) {
                MakeCycleStep(experiment, edges.first);
                weight = chooser_.CountWeight(experiment, edges.second);
                double weight2 = chooser_.CountWeight(experiment, edges.first);
                if (weight > maxWeight || (weight == maxWeight && weight - weight2 > diff)
                        || (weight == maxWeight && weight - weight2 == diff && i == 1)) {
                    maxWeight = weight;
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

    virtual void ResolveShortLoop(BidirectionalPath& path) {
        pair<EdgeId, EdgeId> edges;
        if (GetLoopAndExit(path, edges)) {
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

private:

    double minGapScore_;

    int maxMustHaveOverlap_;

    int maxCanHaveOverlap_;

    int shortOverlap_;

    //int noOverlapGap_;

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


    double ScoreGap(const Sequence& s1, const Sequence& s2, int gap, int initial_gap) const {
        return 1.0 - (double) HammingDistance(s1, s2) / (double) s1.size() - (double) abs(gap - initial_gap) / (double) (2 * g_.k());
    }


public:

    HammingGapJoiner(const Graph& g,
            double minGapScore,
            int mustHaveOverlap,
            int canHaveOverlap,
            int shortOverlap_):
                GapJoiner(g),
                minGapScore_(minGapScore),
                maxMustHaveOverlap_(mustHaveOverlap),
                maxCanHaveOverlap_(canHaveOverlap),
                shortOverlap_(shortOverlap_)
    {
    }

    virtual int FixGap(EdgeId sink, EdgeId source, int initial_gap) const {
        if (initial_gap > (int) g_.k() + maxCanHaveOverlap_) {
            return initial_gap;
        }

        int start = (int) g_.k();
        if (initial_gap < 0) {
            start = (int) g_.k() + min( -initial_gap, (int) min(g_.length(sink), g_.length(source)));
        }

        double max_score = minGapScore_;
        int best_gap = initial_gap;
        bool found = false;
        for (int l = start; l >= shortOverlap_; --l) {
            double score = ScoreGap(g_.EdgeNucls(sink).Subseq((size_t) ((int) g_.length(sink) + (int) g_.k() - l)),
                                    g_.EdgeNucls(source).Subseq(0, (size_t) l),
                                    (int) g_.k() - l,
                                    initial_gap);
            if (score > max_score) {
                max_score = score;
                best_gap = (int) g_.k() - l;
                found = true;
            }
        }

        if (!found) {
            for (int l = shortOverlap_ - 1; l > 0; --l) {
                double score = ScoreGap(g_.EdgeNucls(sink).Subseq((size_t) ((int) g_.length(sink) + (int) g_.k() - l)),
                                        g_.EdgeNucls(source).Subseq(0, (size_t) l),
                                        (int) g_.k() - l,
                                        initial_gap);
                if (score > max_score) {
                    max_score = score;
                    best_gap = (int) g_.k() - l;
                    found = true;
                }
            }
        }

        if (!found) {
            if (initial_gap < maxMustHaveOverlap_) {
                DEBUG("Gap looks like unrealiable: " << initial_gap);
                best_gap = INVALID_GAP;
            }
            else {
                DEBUG("Overlap is not found, initial gap: " << initial_gap << ", not changing.");
                best_gap = initial_gap;
            }
        }
        else {
            DEBUG("Found candidate gap length with score " << max_score);
            DEBUG("Initial: " << initial_gap << ", new gap: " << best_gap);
        }

        return best_gap;
    }

};


class InsertSizeLoopDetector {
protected:
    const Graph& g_;
    size_t min_cycle_len_;

public:
    InsertSizeLoopDetector(const Graph& g, size_t is): g_(g), min_cycle_len_(is) { }

    size_t GetMinCycleLenth() const {
        return min_cycle_len_;
    }

    bool CheckCycled(const BidirectionalPath& path) const {
        return FindCycleStart(path) != -1;
    }
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
        DEBUG("last is pos " << i);
        if (i < 0) return -1;
        BidirectionalPath last = path.SubPath(i);
        last.Print();
        int pos = path.SubPath(0, i).FindFirst(last);
        DEBUG("looking for 1sr IS cycle " << pos);
        return pos;
    }

    int RemoveCycle(BidirectionalPath& path) const {
        int pos = FindCycleStart(path);
        DEBUG("Found IS cycle " << pos);
        if (pos == -1) {
            return -1;
        }

        size_t skip_identical_edges = 0;
        if (path.getLoopDetector().IsCycled(2, skip_identical_edges)) {
            return -1;
            /*int loop_size = (int)path.getLoopDetector().LoopEdges(skip_identical_edges, 2);
            int incorrect_loop_count = (int)path.getLoopDetector().LastLoopCount(loop_size);
            BidirectionalPath small_loop = path.SubPath(path.Size() - loop_size);
            DEBUG("loop size " << loop_size << " incorrect loop count " << incorrect_loop_count);
            small_loop.Print();
            int pos_in_path = (int)path.Size() - loop_size * incorrect_loop_count - 1;
            int pos_in_loop = (int) small_loop.Size() - 1;
            while (pos_in_path >=0 && pos_in_loop >=0 && path.At(pos_in_path) == small_loop.At(pos_in_loop)) {
                pos_in_path--;
                pos_in_loop--;
            }
            pos_in_path++;
            pos_in_loop++;
            DEBUG("pos in path " << pos_in_path << " pos in loop " << pos_in_loop);
            int all_loops_length = (int)path.Size() - pos_in_path;
            DEBUG("all loop len " << all_loops_length);
            int all_loop_count = all_loops_length / loop_size;
            DEBUG("all loop count " << all_loop_count);
            int shift_loop = all_loops_length - all_loop_count * loop_size;
            DEBUG("shift_loop " << shift_loop);
            path.PopBack(shift_loop);
            path.Print();
            DEBUG("Path is cycled after found IS loop, skip identival edges = " << skip_identical_edges);
            path.getLoopDetector().RemoveLoop(skip_identical_edges, false);
            DEBUG("After removing");
            path.Print();
            return pos_in_path;*/
        } else {
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
        set<BidirectionalPath*> cov_paths = cov_map_.GetCoveringPaths(last_e);
        DEBUG("cov paths for e " << g_.int_id(last_e) << " size " << cov_paths.size());
        size_t max_common_size = 0;
        BidirectionalPath* result_p = empty_;
        for (BidirectionalPath* cov_p : cov_paths) {
            if (used_paths_.find(cov_p) == used_paths_.end() || cov_p == &p || cov_p == p.GetConjPath()) {
                continue;
            }
            size_t common_size = MaxCommonSize(p, *cov_p);
            DEBUG("max comon size with path " << cov_p->GetId() << " is " << common_size);
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
        EdgeId last_e = p1.Back();
        vector<size_t> positions2 = p2.FindAll(last_e);
        size_t max_common_size = 0;
        for (size_t pos2 : positions2) {
            size_t common_size = MaxCommonSize(p1, p1.Size() - 1, p2, pos2);
            max_common_size = max(max_common_size, common_size);
        }
        return max_common_size;
    }
private:
    size_t MaxCommonSize(const BidirectionalPath& p1, size_t pos1, const BidirectionalPath& p2, size_t pos2) const {
        int i1 = (int) pos1;
        int i2 = (int) pos2;
        while (i1 >= 0 && i2 >= 0 && p1.At((size_t) i1) == p2.At((size_t) i2)) {
            i1--;
            i2--;
        }
        VERIFY(i1 <= (int)pos1);
        return size_t((int) pos1 - i1);
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

    virtual void GrowAll(PathContainer & paths, PathContainer * result) = 0;

protected:
    const Graph& g_;
    DECL_LOGGER("PathExtender")
};

class PathExtender {
public:
    PathExtender(const Graph & g): g_(g){ }
    virtual ~PathExtender() { }
    virtual bool MakeGrowStep(BidirectionalPath& path) = 0;
protected:
    const Graph& g_;
    DECL_LOGGER("PathExtender")
};

class CompositeExtender : public ContigsMaker {
public:
    CompositeExtender(Graph & g, size_t max_diff_len)
            : ContigsMaker(g),
              coverageMap_(g),
              repeat_detector_(g, coverageMap_, 2 * cfg::get().max_repeat_length),  //TODO: move to config
              extenders_(),
              max_diff_len_(max_diff_len){ }

    CompositeExtender(Graph & g, vector<PathExtender*> pes, size_t max_diff_len)
            : ContigsMaker(g),
              coverageMap_(g),
              repeat_detector_(g, coverageMap_, 2 * cfg::get().max_repeat_length),  //TODO: move to config
              extenders_(),
              max_diff_len_(max_diff_len) {
        extenders_ = pes;
    }

    void AddExender(PathExtender* pe) {
        extenders_.push_back(pe);
    }

    virtual void GrowAll(PathContainer& paths, PathContainer * result) {
        result->clear();
        PathContainer usedPaths;
        for (size_t i = 0; i < paths.size() && !AllPathsCovered(paths); i++) {
            GrowAll(paths, usedPaths, result);
        }

        LengthPathFilter filter(g_, 0);
        filter.filter(*result);
    }

    GraphCoverageMap& GetCoverageMap() {
        return coverageMap_;
    }

    virtual void GrowPath(BidirectionalPath& path) {
        while (MakeGrowStep(path)) { }
    }

    bool MakeGrowStep(BidirectionalPath& path) {
        DEBUG("make grow step composite extender");
//        if (cfg::get().avoid_rc_connections && (path.CameToInterstrandBulge() || path.IsInterstrandBulge())) {
//            DEBUG("Stoping because of interstand bulge");
//            return false;
//        }
        BidirectionalPath* repeat_path = repeat_detector_.RepeatPath(path);
        size_t repeat_size = repeat_detector_.MaxCommonSize(path, *repeat_path);
        if (repeat_size > 0) {
            DEBUG("repeat with length " << repeat_size);
            path.Print();
            repeat_path->Print();
            BidirectionalPath repeat = path.SubPath(path.Size() - repeat_size);
            int begin_repeat = repeat_path->FindLast(repeat);
            VERIFY(begin_repeat > -1);
            size_t end_repeat = (size_t)begin_repeat + repeat_size;
            DEBUG("not consistent subpaths ");
            BidirectionalPath begin1 = path.SubPath(0, path.Size() - repeat_size);
            begin1.Print();
            BidirectionalPath begin2 = repeat_path->SubPath(0, begin_repeat);
            begin2.Print();
            BidirectionalPath end2 = repeat_path->SubPath(end_repeat);
            BidirectionalPath begin1_conj = path.SubPath(0, path.Size() - repeat_size + 1).Conjugate();
            BidirectionalPath begin2_conj = repeat_path->SubPath(0, begin_repeat + 1).Conjugate();
            pair<size_t, size_t> last = ComparePaths(0, 0, begin1_conj, begin2_conj, max_diff_len_);
            DEBUG("last " << last.first << " last2 " << last.second);
            path.Clear();
            repeat_path->Clear();
            if (begin2.Size() == 0 || last.second != 0) { //TODO: incorrect: common edges, but then different ends
               path.PushBack(begin1);
               repeat_path->PushBack(begin2);
            } else {
                path.PushBack(begin2);
                repeat_path->PushBack(begin1);
            }
            path.PushBack(repeat);
            path.PushBack(end2);
            DEBUG("new path");
            path.Print();
            //path.PopBack(repeat_size);
            return false;
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
    GraphCoverageMap coverageMap_;
    RepeatDetector repeat_detector_;
    vector<PathExtender*> extenders_;
    size_t max_diff_len_;

    void SubscribeCoverageMap(BidirectionalPath * path) {
        path->Subscribe(&coverageMap_);
        for (size_t i = 0; i < path->Size(); ++i) {
            coverageMap_.BackEdgeAdded(path->At(i), path, path->GapAt(i));
        }
    }

    bool AllPathsCovered(PathContainer& paths) {
        for (size_t i = 0; i < paths.size(); ++i) {
            if (!coverageMap_.IsCovered(*paths.Get(i)))
                return false;
        }
        return true;
    }

    void GrowAll(PathContainer& paths, PathContainer& usedPaths, PathContainer * result) {
        for (size_t i = 0; i < paths.size(); ++i) {
            if (!coverageMap_.IsCovered(*paths.Get(i))) {
                usedPaths.AddPair(paths.Get(i), paths.GetConjugate(i));
                BidirectionalPath * path = new BidirectionalPath(*paths.Get(i));
                BidirectionalPath * conjugatePath = new BidirectionalPath(*paths.GetConjugate(i));
                result->AddPair(path, conjugatePath);
                SubscribeCoverageMap(path);
                SubscribeCoverageMap(conjugatePath);
                size_t count_trying = 0;
                do {
                    count_trying++;
                    path->CheckGrow();
                    GrowPath(*path);
                    conjugatePath->CheckGrow();
                    GrowPath(*conjugatePath);
                } while (count_trying <10 && (conjugatePath->CheckPrevious() || path->CheckPrevious()));
                path->CheckConjugateEnd();
                DEBUG("result path ");
                path->Print();
            }
        }
    }

    void VerifyMap(PathContainer * result) {
        for (size_t i = 0; i < result->size(); ++i) {
            auto path = result->Get(i);
            for (size_t j = 0; j < path->Size(); ++j) {
                if (coverageMap_.GetCoveringPaths(path->At(j)).count(path) == 0)
                    DEBUG("Inconsistent coverage map");
            }
            path = result->GetConjugate(i);
            for (size_t j = 0; j < path->Size(); ++j) {
                if (coverageMap_.GetCoveringPaths(path->At(j)).count(path) == 0)
                    DEBUG("Inconsistent coverage map");
            }
        }
    }
};


class LoopDetectingPathExtender : public PathExtender {

protected:
    size_t maxLoops_;bool investigateShortLoops_;
    vector<pair<BidirectionalPath*, BidirectionalPath*> > visited_cycles_;
    InsertSizeLoopDetector is_detector_;

public:
    LoopDetectingPathExtender(const Graph & g, size_t max_loops, bool investigateShortLoops, size_t is)
            : PathExtender(g),
              maxLoops_(max_loops),
              investigateShortLoops_(investigateShortLoops),
              is_detector_(g, is) {

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

    bool InExistingLoop(const BidirectionalPath& path) {
        TRACE("Checking existing loops");
        int j = 0;
        for (pair<BidirectionalPath*, BidirectionalPath*> cycle_pair : visited_cycles_) {
            BidirectionalPath* cycle = cycle_pair.first;
            BidirectionalPath* cycle_path = cycle_pair.second;
            VERIFY(!cycle->Empty());
            VERIFY(!cycle_path->Empty());
            VERBOSE_POWER2(j++, "checking ")

            int pos = path.FindLast(*cycle_path);
            if (pos == -1)
                continue;

            int start_cycle_pos = pos + (int) cycle_path->Size();
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
                DEBUG("find cycle " << last_cycle_pos);
                DEBUG("path");
                path.Print();
                DEBUG("last subpath");
                path.SubPath(last_cycle_pos).Print();
                DEBUG("cycle path");
                cycle_path->Print();
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
        int i = (int) pos;
        while (i >= 0 && path.LengthAt(i) < is_detector_.GetMinCycleLenth()) {
            --i;
        }
        if (i < 0)
            i = 0;
        visited_cycles_.push_back(make_pair(new BidirectionalPath(path.SubPath(pos)), new BidirectionalPath(path.SubPath(i))));
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
        DEBUG("check is cycled? " << maxLoops_);
        size_t skip_identical_edges = 0;
        if (path.getLoopDetector().IsCycled(maxLoops_, skip_identical_edges)) {
            size_t loop_size = path.getLoopDetector().LoopEdges(skip_identical_edges, 1);
            DEBUG("Path is Cycled! skip identival edges = " << skip_identical_edges);
            path.Print();
            path.getLoopDetector().RemoveLoop(skip_identical_edges, false);
            DEBUG("After delete");
            path.Print();

            VERIFY(path.Size() >= loop_size);
            AddCycledEdges(path, path.Size() - loop_size);
            return true;
        }
        return false;
    }

    virtual bool MakeSimpleGrowStep(BidirectionalPath& path) = 0;

    virtual bool ResolveShortLoop(BidirectionalPath& path) = 0;

    virtual bool CanInvistigateShortLoop() const {
        return false;
    }

    virtual bool MakeGrowStep(BidirectionalPath& path) {
        if (InExistingLoop(path)) {
            DEBUG("in existing loop");
            return false;
        }
        if (DetectCycle(path)) {
            return false;
        }
        DEBUG("Making step");
        bool result = MakeSimpleGrowStep(path);
        DEBUG("Made step");

        if (DetectCycle(path)) {
            result = false;
        } else if (CanInvistigateShortLoop() && investigateShortLoops_ && path.getLoopDetector().EdgeInShortLoop(path.Back())) {
            DEBUG("Edge in short loop");
            result = ResolveShortLoop(path);
        } else if (CanInvistigateShortLoop() && investigateShortLoops_ && path.getLoopDetector().PrevEdgeInShortLoop()) {
            DEBUG("Prev edge in short loop");
            path.PopBack();
            result = ResolveShortLoop(path);
        }
        return result;
    }

};

class SimpleExtender: public LoopDetectingPathExtender {

protected:

    ExtensionChooser * extensionChooser_;

    LoopResolver loopResolver_;

    void FindFollowingEdges(BidirectionalPath& path, ExtensionChooser::EdgeContainer * result) {
        result->clear();
        vector<EdgeId> edges;
        push_back_all(edges, g_.OutgoingEdges(g_.EdgeEnd(path.Back())));
        result->reserve(edges.size());
        for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
            result->push_back(EdgeWithDistance(*iter, 0));
        }
    }


public:

    SimpleExtender(const Graph& g, ExtensionChooser * ec, size_t is, size_t max_loops, bool investigateShortLoops):
        LoopDetectingPathExtender(g, max_loops, investigateShortLoops, is),
        extensionChooser_(ec),
        loopResolver_(g, *extensionChooser_) {
    }


    virtual bool MakeSimpleGrowStep(BidirectionalPath& path) {
        ExtensionChooser::EdgeContainer candidates;
        FindFollowingEdges(path, &candidates);
        candidates = extensionChooser_->Filter(path, candidates);
        if (candidates.size() == 1) {
            if (!investigateShortLoops_ &&
                    (path.getLoopDetector().EdgeInShortLoop(path.Back())  or path.getLoopDetector().EdgeInShortLoop(candidates.back().e_))
                    && extensionChooser_->WeighConterBased()) {
                return false;
            }
            path.PushBack(candidates.back().e_, candidates.back().d_);
            return true;
        }
        return false;
    }

    virtual bool CanInvistigateShortLoop() const {
        return extensionChooser_->WeighConterBased();
    }

    virtual bool ResolveShortLoop(BidirectionalPath& path) {
        if (extensionChooser_->WeighConterBased()) {
            while (path.getLoopDetector().EdgeInShortLoop(path.Back())) {
                loopResolver_.ResolveShortLoop(path);
            }
            return true;
        }
        return false;
    }

};

class ScaffoldingPathExtender: public LoopDetectingPathExtender {

protected:

    ExtensionChooser * scaffoldingExtensionChooser_;

    //std::vector<int> sizes_;

    ExtensionChooser::EdgeContainer sources_;

    GapJoiner * gapJoiner_;


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

    ScaffoldingPathExtender(const Graph& g, ExtensionChooser * scaffoldingEC, GapJoiner * gapJoiner, size_t is, size_t max_loops, bool investigateShortLoops):
        LoopDetectingPathExtender(g, max_loops, investigateShortLoops, is),
            scaffoldingExtensionChooser_(scaffoldingEC),
            gapJoiner_(gapJoiner)
    {
        InitSources();
    }

    virtual bool MakeSimpleGrowStep(BidirectionalPath& path) {
        DEBUG("scaffolding");
        ExtensionChooser::EdgeContainer candidates;
        bool result = false;

        if (IsSink(path.Back())) {
            candidates = scaffoldingExtensionChooser_->Filter(path, sources_);

            if (candidates.size() == 1) {
                if (candidates[0].e_ == path.Back()) {
                    return false;
                }
                if (cfg::get().avoid_rc_connections && candidates[0].e_ == g_.conjugate(path.Back())) {
                    return false;
                }

                int gap = cfg::get().pe_params.param_set.scaffolder_options.fix_gaps ?
                     gapJoiner_->FixGap(path.Back(), candidates.back().e_, candidates.back().d_) :
                     candidates.back().d_;

                if (gap != GapJoiner::INVALID_GAP) {
                    DEBUG("Scaffolding. PathId: " << path.GetId() << " path length: " << path.Length() << ", fixed gap length: " << gap);
                    path.PushBack(candidates.back().e_, gap);
                    result = true;
                } else {
                    DEBUG("Looks like wrong scaffolding. PathId: " << path.GetId() << " path length: " << path.Length() << ", fixed gap length: " << candidates.back().d_);
                    return false;
                }
            }
        }

        return result;
    }

    virtual bool ResolveShortLoop(BidirectionalPath& /*path*/) {
        return false;
    }

};



}

#endif /* PATH_EXTENDER_HPP_ */
