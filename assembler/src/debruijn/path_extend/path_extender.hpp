//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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

    Graph& g_;

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
            }
            else {
                exit = *edge;
                exit_found = true;
            }
        }


        result = make_pair(loop, exit);
        return exit_found && loop_found;
    }

public:

    ShortLoopResolver(Graph& g): g_(g) {

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


    LoopResolver(Graph& g, ExtensionChooser& chooser): ShortLoopResolver(g), chooser_(chooser) {

    }

    void MakeCycleStep(BidirectionalPath& path, EdgeId e) {
    	EdgeId pathEnd = path.Head();
		path.PushBack(e);
		path.PushBack(pathEnd);
    }

	void MakeBestChoice(BidirectionalPath& path, pair<EdgeId, EdgeId>& edges) {
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
				INFO("make one more circle step i " << i << " weight " << weight << " weight 2 " << weight2);
				if (weight > maxWeight
						|| (weight == maxWeight && weight - weight2 > diff) || (weight == maxWeight && weight - weight2 == diff && i == 1)) {
					maxWeight = weight;
					maxIter = i;
					diff = weight - weight2;
				}
			}
		}
		INFO("mac iter " << maxIter);
		for (size_t i = 0; i < maxIter; ++i) {
			MakeCycleStep(path, edges.first);
		}
		path.PushBack(edges.second);
	}

    virtual void ResolveShortLoop(BidirectionalPath& path) {
        pair<EdgeId, EdgeId> edges;

        if (GetLoopAndExit(path, edges)) {
            DEBUG("Resolving short loop...");
            //path.Print();
            MakeBestChoice(path, edges);
            DEBUG("Resolving short loop done");
            //path.Print();
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

        for (int l = g_.k() ; l > 0; --l) {
            if (g_.EdgeNucls(sink).Subseq(g_.length(sink) + g_.k() - l) == g_.EdgeNucls(source).Subseq(0, l)) {
                DEBUG("Found correct gap length");
                DEBUG("Inintial: " << initial_gap << ", new gap: " << g_.k() - l);
                DEBUG(g_.EdgeNucls(sink).Subseq(g_.length(sink)).str())
                string s = "";
                for (int i = 0; i < (int) g_.k() - l; ++i) {
                    s += " ";
                }
                DEBUG(s << g_.EdgeNucls(source).Subseq(0, g_.k()).str());
                return g_.k() - l;
            }
        }

        string s = "";
        for (int i = 0; i < initial_gap; ++i) {
            s += " ";
        }

        DEBUG("Perfect overlap is not found, inintial: " << initial_gap);
        DEBUG(g_.EdgeNucls(sink).Subseq(g_.length(sink)).str())
        DEBUG(s << g_.EdgeNucls(source).Subseq(0, g_.k()).str());
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
            //int artificalGap):
                GapJoiner(g),
                minGapScore_(minGapScore),
                maxMustHaveOverlap_(mustHaveOverlap),
                maxCanHaveOverlap_(canHaveOverlap),
                shortOverlap_(shortOverlap_)
                //noOverlapGap_(artificalGap)
    {
    }

    virtual int FixGap(EdgeId sink, EdgeId source, int initial_gap) const {
        if (initial_gap > (int) g_.k() + maxCanHaveOverlap_) {
            return initial_gap;
        }

        int start = g_.k();
        if (initial_gap < 0) {
            start = g_.k() + min( -initial_gap, (int) min(g_.length(sink), g_.length(source)));
        }

        double max_score = minGapScore_;
        int best_gap = initial_gap;
        bool found = false;

        for (int l = start; l >= shortOverlap_; --l) {
            double score = ScoreGap(g_.EdgeNucls(sink).Subseq(g_.length(sink) + g_.k() - l), g_.EdgeNucls(source).Subseq(0, l), g_.k() - l, initial_gap);
            if (score > max_score) {
                max_score = score;
                best_gap = (int) g_.k() - l;
                found = true;
            }
        }

        if (!found) {
            for (int l = shortOverlap_ - 1; l > 0; --l) {
                double score = ScoreGap(g_.EdgeNucls(sink).Subseq(g_.length(sink) + g_.k() - l), g_.EdgeNucls(source).Subseq(0, l), g_.k() - l, initial_gap);
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
                //best_gap = max(initial_gap, (int) g_.k() + noOverlapGap_);
                best_gap = initial_gap;
            }
        }
        else {
            DEBUG("Found candidate gap length with score " << max_score);
            DEBUG("Initial: " << initial_gap << ", new gap: " << best_gap);
        }


        string s = "";
        for (int i = 0; i < best_gap; ++i) {
            s += " ";
        }

        DEBUG(g_.EdgeNucls(sink).Subseq(g_.length(sink)).str())
        DEBUG(s << g_.EdgeNucls(source).Subseq(0, g_.k()).str());
        return best_gap;
    }

};


class PathExtender {

protected:
    Graph& g_;

public:
    PathExtender(Graph & g): g_(g)
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
    LoopDetectingPathExtender(Graph & g, size_t max_loops): PathExtender(g), maxLoops_(max_loops), investigateShortLoops_(true)
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
            if (path.getLoopDetector().IsCycled(maxLoops_)) {
                path.getLoopDetector().RemoveLoop();
                return;
            }
        }
    }

    virtual void GrowAll(PathContainer& paths, PathContainer * result) {
        result->clear();

        for (size_t i = 0; i < paths.size(); i ++) {
            BidirectionalPath * path = new BidirectionalPath(*paths.Get(i));
            path->SetCurrentPathAsSeed();
            BidirectionalPath * conjugatePath = new BidirectionalPath(*paths.GetConjugate(i));
            conjugatePath->SetCurrentPathAsSeed();

            result->AddPair(path, conjugatePath);

            do {
                path->CheckGrow();
                GrowPath(*path);
                conjugatePath->CheckGrow();
                GrowPath(*conjugatePath);
            } while (conjugatePath->CheckPrevious() || path->CheckPrevious());
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
                path->SetCurrentPathAsSeed();
                BidirectionalPath * conjugatePath = new BidirectionalPath(*paths.GetConjugate(i));
                conjugatePath->SetCurrentPathAsSeed();
                result->AddPair(path, conjugatePath);
                SubscribeCoverageMap(path);
                SubscribeCoverageMap(conjugatePath);

                if (!coverageMap_.IsCovered(*path) || !coverageMap_.IsCovered(*conjugatePath)) {
                    DEBUG("Paths are not covered after subsciption");
                }

                do {
					path->CheckGrow();
					GrowPath(*path);
					//verifyMap(result);
					conjugatePath->CheckGrow();
					GrowPath(*conjugatePath);
					//verifyMap(result);
                } while (conjugatePath->CheckPrevious() || path->CheckPrevious());

                if (!coverageMap_.IsCovered(*paths.Get(i)) || !coverageMap_.IsCovered(*paths.GetConjugate(i))) {
                    DEBUG("Seeds are not covered after growing");
                }
                path->CheckConjugateEnd();
            }
        }
    }

    void RemoveSubpaths(PathContainer& usedPaths) {
        for (size_t i = 0; i < usedPaths.size(); ++i) {
            if (coverageMap_.GetUniqueCoverage(*usedPaths.Get(i)) > 1) {

                size_t seedId = usedPaths.Get(i)->GetId();
                size_t seedConjId = usedPaths.GetConjugate(i)->GetId();

                std::set<BidirectionalPath*> coveringPaths = coverageMap_.GetCoveringPaths(*(usedPaths.Get(i)));
                for (auto iter = coveringPaths.begin(); iter != coveringPaths.end(); ++iter) {

                    if ((*iter)->GetId() == seedId) {
                        bool otherSeedFound = false;
                        for (auto it = coveringPaths.begin(); it != coveringPaths.end(); ++it) {
                            if ((*it)->GetId() != seedId && (*it)->GetId() != seedConjId) {
                                otherSeedFound = true;
                                if ((*it)->Length() < (*iter)->Length()) {
                                    DEBUG("Covering path is shorter than the original seed path; seed path is to be removed");
                                }
                                if (!(*it)->Contains(**iter)) {
                                    DEBUG("Not subpaths");
                                    DEBUG((*iter)->Length() << " " << (*it)->Length());
                                    (*iter)->Print();
                                    (*it)->Print();
                                }
                            }
                        }

                        if (otherSeedFound) {
                            (*iter)->Clear();
                        }
                    }
                }
            }
        }
    }

public:

    CoveringPathExtender(Graph& g_, size_t max_loops): LoopDetectingPathExtender(g_, max_loops), coverageMap_(g_) {
    }


    virtual void GrowAll(PathContainer& paths, PathContainer * result) {
        result->clear();
        PathContainer usedPaths;

        for (size_t i = 0; i < paths.size() && !AllPathsCovered(paths); i ++) {
		    GrowAll(paths, usedPaths, result);
		    //RemoveSubpaths(usedPaths);
        }

        LengthPathFilter filter(g_, 0);
        filter.filter(*result);
    }

    GraphCoverageMap& GetCoverageMap() {
        return coverageMap_;
    }

};


class CompositePathExtender: public CoveringPathExtender {


protected:

    vector<PathExtender* > extenders_;

public:

    CompositePathExtender(Graph & g, size_t max_loops): CoveringPathExtender(g, max_loops), extenders_() {
    }

    void AddExender(PathExtender* pe) {
        extenders_.push_back(pe);
    }

    CompositePathExtender(Graph & g, size_t max_loops, PathExtender* pe1): CoveringPathExtender(g, max_loops), extenders_() {
        AddExender(pe1);
    }

    CompositePathExtender(Graph & g, size_t max_loops, PathExtender* pe1, PathExtender* pe2): CoveringPathExtender(g, max_loops), extenders_() {
        AddExender(pe1);
        AddExender(pe2);
    }

    CompositePathExtender(Graph & g, size_t max_loops, PathExtender* pe1, PathExtender* pe2, PathExtender* pe3): CoveringPathExtender(g, max_loops), extenders_() {
        AddExender(pe1);
        AddExender(pe2);
        AddExender(pe3);
    }

    CompositePathExtender(Graph & g, size_t max_loops, vector<PathExtender*> pes) :
			CoveringPathExtender(g, max_loops), extenders_() {
		extenders_ = pes;
	}

    virtual bool MakeGrowStep(BidirectionalPath& path) {
        size_t current = 0;

        while (current < extenders_.size()) {
            if (extenders_[current]->MakeGrowStep(path)) {
                return true;
            }
            ++current;
        }
        return false;
    }

};




class SimplePathExtender: public CoveringPathExtender {

protected:

    ExtensionChooser * extensionChooser_;

    LoopResolver loopResolver_;

    void FindFollowingEdges(BidirectionalPath& path, ExtensionChooser::EdgeContainer * result) {
        result->clear();
        auto edges = g_.OutgoingEdges(g_.EdgeEnd(path.Back()));
        result->reserve(edges.size());
        for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
            result->push_back(EdgeWithDistance(*iter, 0));
        }
    }


public:

    SimplePathExtender(Graph& g, size_t max_loops, ExtensionChooser * ec): CoveringPathExtender(g, max_loops), extensionChooser_(ec), loopResolver_(g, *extensionChooser_) {
    }


    virtual bool MakeGrowStep(BidirectionalPath& path) {
        ExtensionChooser::EdgeContainer candidates;
        bool result = false;
        FindFollowingEdges(path, &candidates);
        candidates = extensionChooser_->Filter(path, candidates);

        if (candidates.size() == 1) {
            path.PushBack(candidates.back().e_, candidates.back().d_);
            result = true;

            if (investigateShortLoops_ && path.getLoopDetector().EdgeInShortLoop() && extensionChooser_->WeighConterBased()) {
                loopResolver_.ResolveShortLoop(path);
            }
        } else if (candidates.size() >= 1){
        	INFO("MORE 1 CANDIDATE");
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

    ScaffoldingPathExtender(Graph& g, size_t max_loops, ExtensionChooser * scaffoldingEC, GapJoiner * gapJoiner): CoveringPathExtender(g, max_loops),
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
                DEBUG(candidates.size() << " " << g_.int_id(candidates[0].e_) << " Path id :" << path.GetId()<< "  Edge len : " << g_.length(candidates[0].e_))

                int gap = params.param_set.scaffolder_options.fix_gaps ?
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
