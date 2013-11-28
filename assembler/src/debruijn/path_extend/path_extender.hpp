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

protected:

    const Graph& g_;

    bool GetLoopAndExit(BidirectionalPath& path,
                        pair<EdgeId, EdgeId>& result) const {
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

public:

    ShortLoopResolver(const Graph& g): g_(g) {

    }

    virtual ~ShortLoopResolver() {

    }

    virtual void ResolveShortLoop(BidirectionalPath& path) = 0;

protected:
    DECL_LOGGER("PathExtender")

};


class SimpleLoopResolver: public ShortLoopResolver {

public:


    SimpleLoopResolver(Graph& g): ShortLoopResolver(g) {

    }

    virtual void ResolveShortLoop(BidirectionalPath& path) {
        pair<EdgeId, EdgeId> edges;

        if (GetLoopAndExit(path, edges)) {
            DEBUG("Resolving short loop...");
            path.Print();

            EdgeId e = path.Head();
            path.PushBack(edges.first);
            path.PushBack(e);
            path.PushBack(edges.second);
            DEBUG("Resolving short loop done");

            path.Print();
        }
    }

protected:
    DECL_LOGGER("PathExtender")

};


class LoopResolver: public ShortLoopResolver {

    static const size_t iter_ = 10;
    ExtensionChooser& chooser_;

public:


    LoopResolver(const Graph& g, ExtensionChooser& chooser): ShortLoopResolver(g), chooser_(chooser) {

    }

    void MakeCycleStep(BidirectionalPath& path, EdgeId e) {
    	EdgeId pathEnd = path.Head();
		path.PushBack(e);
		path.PushBack(pathEnd);
    }

    void MakeBestChoice(BidirectionalPath& path, pair<EdgeId, EdgeId>& edges) {
        DEBUG("Path before deleting");
        path.Print();
        EdgeId first_edge = path.Back();
        EdgeId second_edge = edges.first;
        while (path.Size() > 2) {
            if (path.At(path.Size() - 1) == first_edge
                    && path.At(path.Size() - 2) == second_edge) {
                path.PopBack(2);
            } else {
                break;
            }
        }
        DEBUG("Path after deleting");
        path.Print();
        chooser_.ClearExcludedEdges();
        BidirectionalPath experiment(path);
        double maxWeight = chooser_.CountWeight(experiment, edges.second);
        double diff = maxWeight - chooser_.CountWeight(experiment, edges.first);
        size_t maxIter = 0;
        for (size_t i = 1; i <= iter_; ++i) {
            double weight = chooser_.CountWeight(experiment, edges.first);
            DEBUG("weight " << weight);
            if (weight > 0) {
                MakeCycleStep(experiment, edges.first);
                weight = chooser_.CountWeight(experiment, edges.second);
                double weight2 = chooser_.CountWeight(experiment, edges.first);
                DEBUG("iter " << i << " weight " << weight  << " maxWeight " << maxWeight << " weight 2 " <<  weight2 << " diff " << diff);
                if (weight > maxWeight ||
                        (weight == maxWeight && weight - weight2 > diff) ||
                        (weight == maxWeight && weight - weight2 == diff  && i == 1)) {
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
        DEBUG("path after resolving");
        path.Print();
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

protected:

    const Graph& g_;

public:

    GapJoiner(const Graph& g): g_(g) {
    }

    virtual int FixGap(EdgeId sink, EdgeId source, int initial_gap) const = 0;

    virtual ~GapJoiner() {

    }

    static const int INVALID_GAP = -1000000;

protected:
    DECL_LOGGER("PathExtender")


};


class SimpleGapJoiner: public GapJoiner {

public:

    SimpleGapJoiner(const Graph& g): GapJoiner(g) {
    }

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


class PathExtender {

protected:
    const Graph& g_;

public:
    PathExtender(const Graph & g): g_(g)
    {
    }

    virtual ~PathExtender() {

    }

    virtual bool MakeGrowStep(BidirectionalPath& path) = 0;

    virtual void GrowPath(BidirectionalPath& path) = 0;

    virtual void GrowAll(PathContainer & paths, PathContainer * result) = 0;

protected:
    DECL_LOGGER("PathExtender")
};


class CoveringPathExtender: public PathExtender {

protected:

    GraphCoverageMap coverageMap_;

    void SubscribeCoverageMap(BidirectionalPath * path) {
        path->Subscribe(&coverageMap_);
        for (size_t i = 0; i < path->Size(); ++i) {
            coverageMap_.BackEdgeAdded(path->At(i), path, path->GapAt(i));
        }
    }

    bool AllPathsCovered(PathContainer& paths) {
        for (size_t i = 0; i < paths.size(); ++i) {
            if (!coverageMap_.IsCovered(*paths.Get(i))) {
                return false;
            }
        }
        return true;
    }

    // DEBUG
    void VerifyMap(PathContainer * result) {
        //MAP VERIFICATION
        for (size_t i = 0; i < result->size(); ++i) {
            auto path = result->Get(i);
            for (size_t j = 0; j < path->Size(); ++j) {
                if (coverageMap_.GetCoveringPaths(path->At(j)).count(path) == 0) {
                    DEBUG("Inconsistent coverage map");
                }
            }

            path = result->GetConjugate(i);
            for (size_t j = 0; j < path->Size(); ++j) {
                if (coverageMap_.GetCoveringPaths(path->At(j)).count(path) == 0) {
                    DEBUG("Inconsistent coverage map");
                }
            }
        }
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

                do {
                    path->CheckGrow();
                    GrowPath(*path);
                    conjugatePath->CheckGrow();
                    GrowPath(*conjugatePath);
                }
                while (conjugatePath->CheckPrevious() || path->CheckPrevious());
                path->CheckConjugateEnd();
                DEBUG("result path ");
                path->Print();
            }
        }
    }

public:

    CoveringPathExtender(const Graph& g_): PathExtender(g_), coverageMap_(g_) {
    }

    virtual void GrowPath(BidirectionalPath& path) {
        while (MakeGrowStep(path)) {
        }
    }

    virtual void GrowAll(PathContainer& paths, PathContainer * result) {
        result->clear();
        PathContainer usedPaths;

        for (size_t i = 0; i < paths.size() && !AllPathsCovered(paths); i ++) {
            GrowAll(paths, usedPaths, result);
        }

        LengthPathFilter filter(g_, 0);
        filter.filter(*result);
    }

    GraphCoverageMap& GetCoverageMap() {
        return coverageMap_;
    }

};


class InsertSizeLoopDetector {
protected:
    const Graph& g_;

    size_t min_cycle_len_;

public:
    InsertSizeLoopDetector(const Graph& g, size_t is): g_(g), min_cycle_len_(is) {
    }

    size_t GetMinCycleLenth() const {
        return min_cycle_len_;
    }

    bool CheckCycled(const BidirectionalPath& path) const {
        return FindCycleStart(path) != -1;
    }

    int FindCycleStart(const BidirectionalPath& path) const {

        int i = (int) path.Size() - 1;
        DEBUG("Looking for IS cycle " << min_cycle_len_);
        while (i >= 0 && path.LengthAt(i) < min_cycle_len_) {
            --i;
        }

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
        int last_edge_pos = path.FindLast(path[pos]);
        DEBUG("last edge pos " << last_edge_pos);
        path.Print();
        VERIFY(last_edge_pos > pos);
        for (int i = (int) path.Size() - 1; i >= last_edge_pos; --i) {
            path.PopBack();
        }
        path.Print();
        VERIFY((int) path.Size() == last_edge_pos);

        size_t skip_identical_edges = 0;
        if (path.getLoopDetector().IsCycled(2, skip_identical_edges)) {
            DEBUG("Path is cycled after found IS loop, skip identival edges = " << skip_identical_edges);
            path.Print();
            path.getLoopDetector().RemoveLoop(skip_identical_edges, false);
            DEBUG("After removing");
            path.Print();
        }
        VERIFY(pos < (int) path.Size());
        return pos;
    }
};

class LoopDetectingPathExtender: public CoveringPathExtender {

protected:
    size_t maxLoops_;
    bool investigateShortLoops_;
    vector< pair<BidirectionalPath*, BidirectionalPath*> > visited_cycles_;
    InsertSizeLoopDetector is_detector_;

public:
    LoopDetectingPathExtender(const Graph & g, size_t max_loops,
            bool investigateShortLoops,
            size_t is):
        CoveringPathExtender(g),
        maxLoops_(max_loops),
        investigateShortLoops_(investigateShortLoops),
        is_detector_(g, is)
    {
    }

    size_t getMaxLoops() const
    {
        return maxLoops_;
    }

    bool isInvestigateShortLoops() const
    {
        return investigateShortLoops_;
    }

    void setInvestigateShortLoops(bool investigateShortLoops)
    {
        this->investigateShortLoops_ = investigateShortLoops;
    }

    void setMaxLoops(size_t maxLoops)
    {
        if (maxLoops != 0) {
            this->maxLoops_ = maxLoops;
        }
    }

    bool InExistingLoop(const BidirectionalPath& path) {
        DEBUG("Checking existing loops");
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
            for (int i = start_cycle_pos; i < (int) path.Size() - (int) cycle->Size(); i += (int) cycle->Size()) {
                if (!path.CompareFrom(i, *cycle)) {
                    only_cycles_in_tail = false;
                    break;
                }
                else {
                    last_cycle_pos = i;
                }
            }
            only_cycles_in_tail = only_cycles_in_tail && cycle->CompareFrom(0, path.SubPath(last_cycle_pos + (int) cycle->Size()));
            if (only_cycles_in_tail) {
                return true;
            }
        }
        return false;
    }

    void AddCycledEdges(const BidirectionalPath& path, size_t pos) {
        if (pos >= path.Size()) {
            WARN("Wrong position in IS cycle");
            return;
        }
        int i = (int) pos;
        while (i >= 0 && path.LengthAt(i) < is_detector_.GetMinCycleLenth()) {
            --i;
        }
        if (i < 0)
            i = 0;

        visited_cycles_.push_back(make_pair(new BidirectionalPath(path.SubPath(pos)),
                                            new BidirectionalPath(path.SubPath(i))));
        DEBUG("add cycle");
        path.SubPath(pos).Print();
    }

    bool DetectCycle(BidirectionalPath& path) {
        size_t skip_identical_edges = 0;
        if (is_detector_.CheckCycled(path)) {
            DEBUG("Checking IS cycle");
            int loop_pos = is_detector_.RemoveCycle(path);
            DEBUG("Removed IS cycle");
            VERIFY(loop_pos != -1);
            AddCycledEdges(path, loop_pos);
            DEBUG("Added IS cycle");
            return true;
        } else if (path.getLoopDetector().IsCycled(maxLoops_, skip_identical_edges)) {
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

    virtual bool MakeGrowStep(BidirectionalPath& path) {
        if (InExistingLoop(path)) {
            return false;
        }

        DEBUG("Making step");
        bool result = MakeSimpleGrowStep(path);
        DEBUG("Made step");


        if (DetectCycle(path)) {
            DEBUG("True cycle");
            result = false;
        }
        else if (investigateShortLoops_ && path.getLoopDetector().EdgeInShortLoop(path.Back())) {
            DEBUG("Edge in short loop");
            result = ResolveShortLoop(path);
        }
        else if (investigateShortLoops_ && path.getLoopDetector().PrevEdgeInShortLoop()) {
            DEBUG("Prev edge in short loop");
            path.PopBack();
            result = ResolveShortLoop(path);
        }

        return result;
    }

    virtual void GrowAll(PathContainer& paths, PathContainer * result) {
        result->clear();

        for (size_t i = 0; i < paths.size(); i ++) {
            BidirectionalPath * path = new BidirectionalPath(*paths.Get(i));
            BidirectionalPath * conjugatePath = new BidirectionalPath(*paths.GetConjugate(i));
            result->AddPair(path, conjugatePath);

            do {
                path->CheckGrow();
                GrowPath(*path);
                conjugatePath->CheckGrow();
                GrowPath(*conjugatePath);
            }
            while (conjugatePath->CheckPrevious() || path->CheckPrevious());
        }
    }

};


class CompositeExtender: public CoveringPathExtender {

protected:

    vector<CoveringPathExtender* > extenders_;

public:

    CompositeExtender(Graph & g): CoveringPathExtender(g), extenders_() {
    }

    void AddExender(LoopDetectingPathExtender* pe) {
        extenders_.push_back(pe);
    }

    CompositeExtender(Graph & g, vector<CoveringPathExtender*> pes) :
			CoveringPathExtender(g), extenders_() {
		extenders_ = pes;
	}

    virtual bool MakeGrowStep(BidirectionalPath& path) {
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

    virtual bool ResolveShortLoop(BidirectionalPath& path) {
        if (extensionChooser_->WeighConterBased()) {
            DEBUG("Resolving short loop")
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
