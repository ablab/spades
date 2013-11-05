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

    Graph& g_;

public:

    GapJoiner(Graph& g): g_(g) {
    }

    virtual int FixGap(EdgeId sink, EdgeId source, int initial_gap) const = 0;

    virtual ~GapJoiner() {

    }

    static const int INVALID_GAP = -1000000;

};


class SimpleGapJoiner: public GapJoiner {

public:

    SimpleGapJoiner(Graph& g): GapJoiner(g) {
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

    HammingGapJoiner(Graph& g,
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
//            else if (initial_gap < (int) g_.k()) {
//                best_gap = (int) g_.k() + noOverlapGap_;
//                DEBUG("Overlap is not found, initial gap: " << initial_gap << ", changing to " << best_gap);
//            }
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
};



class LoopDetectingPathExtender: public PathExtender {

protected:
    size_t maxLoops_;

    bool investigateShortLoops_;

public:
    LoopDetectingPathExtender(const Graph & g, size_t max_loops, bool investigateShortLoops): PathExtender(g), maxLoops_(max_loops), investigateShortLoops_(investigateShortLoops)
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

    virtual void GrowPath(BidirectionalPath& path) {
        while (MakeGrowStep(path)) {
            size_t skip_identical_edges = 0;
            if (path.getLoopDetector().IsCycled(maxLoops_, skip_identical_edges)) {
                DEBUG("Path is Cycled!");
                DEBUG("skip identival edges = " << skip_identical_edges);
                path.Print();
                path.getLoopDetector().RemoveLoop(skip_identical_edges, false);
                DEBUG("After delete");
                path.Print();
                return;
            }
        }
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


class CoveringPathExtender: public LoopDetectingPathExtender {

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

//                if (!coverageMap_.IsCovered(*path) || !coverageMap_.IsCovered(*conjugatePath)) {
//                    DEBUG("Paths are not covered after subsciption");
//                }

                do {
					path->CheckGrow();
					GrowPath(*path);
					conjugatePath->CheckGrow();
					GrowPath(*conjugatePath);
                }
                while (conjugatePath->CheckPrevious() || path->CheckPrevious());

//                if (!coverageMap_.IsCovered(*paths.Get(i)) || !coverageMap_.IsCovered(*paths.GetConjugate(i))) {
//                    DEBUG("Seeds are not covered after growing");
//                }
                path->CheckConjugateEnd();
                DEBUG("result path ");
                path->Print();
            }
        }
    }

public:

    CoveringPathExtender(const Graph& g_, size_t max_loops, bool investigateShortLoops): LoopDetectingPathExtender(g_, max_loops, investigateShortLoops), coverageMap_(g_) {
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


class CompositeExtender: public CoveringPathExtender {


protected:

    vector<PathExtender* > extenders_;

public:

    CompositeExtender(Graph & g, size_t max_loops, bool investigateShortLoops): CoveringPathExtender(g, max_loops, investigateShortLoops), extenders_() {
    }

    void AddExender(PathExtender* pe) {
        extenders_.push_back(pe);
    }

    CompositeExtender(Graph & g, size_t max_loops, vector<PathExtender*> pes, bool investigateShortLoops = true) :
			CoveringPathExtender(g, max_loops, investigateShortLoops), extenders_() {
		extenders_ = pes;
	}

    virtual bool MakeGrowStep(BidirectionalPath& path) {
        if (cfg::get().avoid_rc_connections && (path.CameToInterstrandBulge() || path.IsInterstrandBulge())) {
            DEBUG("Stoping because of interstand bulge");
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

};




class SimpleExtender: public CoveringPathExtender {

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

    SimpleExtender(const Graph& g, size_t max_loops, ExtensionChooser * ec,  bool investigateShortLoops = true):
    	CoveringPathExtender(g, max_loops, investigateShortLoops), extensionChooser_(ec), loopResolver_(g, *extensionChooser_) {
    }

    virtual bool MakeGrowStep(BidirectionalPath& path) {
        ExtensionChooser::EdgeContainer candidates;
        bool result = false;
        FindFollowingEdges(path, &candidates);
        candidates = extensionChooser_->Filter(path, candidates);
        if (candidates.size() == 1) {
            if (!investigateShortLoops_
                    && (path.getLoopDetector().EdgeInShortLoop(path.Back())
                            or path.getLoopDetector().EdgeInShortLoop(
                                    candidates.back().e_))
                    && extensionChooser_->WeighConterBased()) {
                return false;
            }
            path.PushBack(candidates.back().e_, candidates.back().d_);
            result = true;
            if (investigateShortLoops_
                    && path.getLoopDetector().EdgeInShortLoop(path.Back())
                    && extensionChooser_->WeighConterBased()) {
                while (path.getLoopDetector().EdgeInShortLoop(path.Back())) {
                    loopResolver_.ResolveShortLoop(path);
                }
            }
        } else if (investigateShortLoops_
                && path.getLoopDetector().PrevEdgeInShortLoop()
                && extensionChooser_->WeighConterBased()) {
            DEBUG("Prev edge in short loop");
            path.PopBack();
            while (path.getLoopDetector().EdgeInShortLoop(path.Back())) {
                loopResolver_.ResolveShortLoop(path);
            }
            result = true;
        } else if (investigateShortLoops_
                && path.getLoopDetector().EdgeInShortLoop(path.Back())
                && extensionChooser_->WeighConterBased()) {
            DEBUG("Edge in short loop");
            while (path.getLoopDetector().EdgeInShortLoop(path.Back())) {
                loopResolver_.ResolveShortLoop(path);
            }
            result = true;
        } else if (candidates.size() >= 1) {
            DEBUG("MORE 1 CANDIDATE");
        }
        return result;
    }

};



class ScaffoldingPathExtender: public CoveringPathExtender {

protected:

    ExtensionChooser * scaffoldingExtensionChooser_;

    //std::vector<int> sizes_;

    ExtensionChooser::EdgeContainer sources_;

    GapJoiner * gapJoiner_;


    void InitSources() {
        sources_.clear();

        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (g_.IncomingEdgeCount(g_.EdgeStart(*iter)) == 0) {
                sources_.push_back(EdgeWithDistance(*iter, 0));
            }
        }
    }

    bool IsSink(EdgeId e)
	{
		return g_.OutgoingEdgeCount(g_.EdgeEnd(e)) == 0;
	}


public:

    ScaffoldingPathExtender(Graph& g, size_t max_loops, ExtensionChooser * scaffoldingEC, GapJoiner * gapJoiner, bool investigateShortLoops = true):
    	CoveringPathExtender(g, max_loops, investigateShortLoops),
            scaffoldingExtensionChooser_(scaffoldingEC),
            gapJoiner_(gapJoiner)
    {
        InitSources();
    }


    virtual bool MakeGrowStep(BidirectionalPath& path) {
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
                DEBUG(candidates.size() << " " << g_.int_id(candidates[0].e_) << " Path id :" << path.GetId()<< "  Edge len : " << g_.length(candidates[0].e_))

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

};



}

#endif /* PATH_EXTENDER_HPP_ */
