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

#include "extension_chooser.hpp"
#include "path_filter.hpp"

namespace path_extend {

class GraphCoverageMap: public PathListener {

public:
    typedef std::multiset <BidirectionalPath *> MapDataT;


protected:
    Graph& g_;

    std::map <EdgeId, MapDataT * > edgeCoverage_;

    MapDataT * empty_;

    virtual void EdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
        auto iter = edgeCoverage_.find(e);
        if (iter == edgeCoverage_.end()) {
            edgeCoverage_.insert(std::make_pair(e, new MapDataT()));
        }
        edgeCoverage_[e]->insert(path);
    }

    virtual void EdgeRemoved(EdgeId e, BidirectionalPath * path) {
        auto iter = edgeCoverage_.find(e);
        if (iter != edgeCoverage_.end()) {
            if (iter->second->count(path) == 0) {
                DEBUG("Error erasing path from coverage map");
            } else {
                auto entry = iter->second->find(path);
                iter->second->erase(entry);
            }
        }
    }

public:
    GraphCoverageMap(Graph& g_) : g_(g_), edgeCoverage_() {
        empty_ = new MapDataT();
    }

    virtual void FrontEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
        EdgeAdded(e, path, gap);
    }

    virtual void BackEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
        EdgeAdded(e, path, gap);
    }

    virtual void FrontEdgeRemoved(EdgeId e, BidirectionalPath * path) {
        EdgeRemoved(e, path);
    }

    virtual void BackEdgeRemoved(EdgeId e, BidirectionalPath * path) {
        EdgeRemoved(e, path);
    }

    MapDataT * GetEdgePaths(EdgeId e) const {
        auto iter = edgeCoverage_.find(e);
        if (iter != edgeCoverage_.end()) {
            return iter->second;
        }

        return empty_;
    }


    int GetCoverage(EdgeId e) const {
        return GetEdgePaths(e)->size();
    }


    bool IsCovered(EdgeId e) const {
        return GetCoverage(e) > 0;
    }

    bool IsCovered(const BidirectionalPath& path) const {
        for (size_t i = 0; i < path.Size(); ++i) {
            if (!IsCovered(path[i])) {
                return false;
            }
        }
        return true;
    }

    int GetCoverage(const BidirectionalPath& path) const {
        if (path.Empty()) {
            return 0;
        }

        int cov = GetCoverage(path[0]);
        for (size_t i = 1; i < path.Size(); ++i) {
            int currentCov = GetCoverage(path[i]);
            if (cov > currentCov) {
                cov = currentCov;
            }
        }

        return cov;
    }

    std::set<BidirectionalPath*> GetCoveringPaths(EdgeId e) const {
        auto mapData = GetEdgePaths(e);
        return std::set<BidirectionalPath*>(mapData->begin(), mapData->end());

    }

    std::set<BidirectionalPath*> GetCoveringPaths(const BidirectionalPath& path) const {
        std::set<BidirectionalPath*> result;

        if (!path.Empty()) {
            MapDataT * data;
            data = GetEdgePaths(path.Front());

            result.insert(data->begin(), data->end());

            for (size_t i = 1; i < path.Size(); ++i) {
                data = GetEdgePaths(path[i]);

                std::set<BidirectionalPath*> dataSet;
                dataSet.insert(data->begin(), data->end());

                for (auto iter = result.begin(); iter != result.end(); ) {
                    auto next = iter;
                    ++next;
                    if (dataSet.count(*iter) == 0) {
                        result.erase(iter);
                    }
                    iter = next;
                }
            }
        }

        return result;
    }

    int GetUniqueCoverage(EdgeId e) const {
        return GetCoveringPaths(e).size();
    }

    int GetUniqueCoverage(const BidirectionalPath& path) const {
        return GetCoveringPaths(path).size();
    }

    std::map <EdgeId, MapDataT * >::const_iterator begin() const {
        return edgeCoverage_.begin();
    }

    std::map <EdgeId, MapDataT * >::const_iterator end() const {
        return edgeCoverage_.end();
    }

    // DEBUG

    void PrintUncovered() const {
        DEBUG("Uncovered edges");
        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (!IsCovered(*iter)) {
                DEBUG(g_.int_id(*iter) << " (" << g_.length(*iter) << ") ~ " << g_.int_id(g_.conjugate(*iter)) << " (" << g_.length(g_.conjugate(*iter)) << ")");
            }
        }
    }
};


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
         size_t  maxIter = 0;
         for (size_t i = 1; i <= iter_; ++ i) {
        	 double weight = chooser_.CountWeight(experiment, edges.first);
        	 if (weight > 0) {
				 MakeCycleStep(experiment, edges.first);
				 weight = chooser_.CountWeight(experiment, edges.second);
				 double weight2 = chooser_.CountWeight(experiment, edges.first);
				 //DEBUG("now weight is " << weight << " dif w is: " << weight - weight2)
				 if (weight > maxWeight || (weight == maxWeight && weight - weight2 > diff) ) {
					 maxWeight = weight;
					 maxIter = i;
					 diff = weight - weight2;
				 }
        	 }
         }
         for (size_t i = 0; i < maxIter; ++ i) {
        	 MakeCycleStep(path, edges.first);
         }
         //DEBUG("Max number of iterations: " << maxIter);
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

    size_t maxLoops_;

    bool investigateShortLoops_;

public:
    PathExtender(Graph & g): g_(g), maxLoops_(10), investigateShortLoops_(true)
    {
    }

    virtual ~PathExtender() {

    }

    virtual void GrowAll(PathContainer & paths, PathContainer * result) = 0;

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


    virtual void GrowPath(BidirectionalPath& path) = 0;

public:

    CoveringPathExtender(Graph& g_): PathExtender(g_), coverageMap_(g_) {
    }

    virtual void GrowAll(PathContainer& paths, PathContainer * result) {
        result->clear();
        PathContainer usedPaths;

        for (size_t i = 0; i < paths.size() && !AllPathsCovered(paths); i ++) {
		    GrowAll(paths, usedPaths, result);
		    RemoveSubpaths(usedPaths);
        }

        LengthPathFilter filter(g_, 0);
        filter.filter(*result);
    }

    GraphCoverageMap& GetCoverageMap() {
        return coverageMap_;
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

    virtual void GrowPath(BidirectionalPath& path) {
        ExtensionChooser::EdgeContainer candidates;
        do {
            FindFollowingEdges(path, &candidates);
            candidates = extensionChooser_->Filter(path, candidates);

            if (candidates.size() == 1) {
                path.PushBack(candidates.back().e_, candidates.back().d_);

                if (investigateShortLoops_ && path.getLoopDetector().EdgeInShortLoop()) {
                    loopResolver_.ResolveShortLoop(path);
                }
            }

            if (path.getLoopDetector().IsCycled(maxLoops_)) {
                path.getLoopDetector().RemoveLoop();
                break;
            }

        } while (candidates.size() == 1);
    }

public:

    SimplePathExtender(Graph& g, ExtensionChooser * ec): CoveringPathExtender(g), extensionChooser_(ec), loopResolver_(g, *extensionChooser_) {
    }

};



class ScaffoldingPathExtender: public SimplePathExtender {

protected:

    ExtensionChooser * scaffoldingExtensionChooser_;

    std::vector<int> sizes_;

    ExtensionChooser::EdgeContainer sources_;

    HammingGapJoiner gapJoiner_;


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

    virtual void GrowPath(BidirectionalPath& path)
    {
        ExtensionChooser::EdgeContainer candidates;
        do {
        	FindFollowingEdges(path, &candidates);
			candidates = extensionChooser_->Filter(path, candidates);

			if (candidates.size() == 1) {
				path.PushBack(candidates.back().e_, candidates.back().d_);

	            if (investigateShortLoops_ && path.getLoopDetector().EdgeInShortLoop()) {
	                loopResolver_.ResolveShortLoop(path);
	            }
			}
			else if (IsSink(path.Back())) {
            	candidates = scaffoldingExtensionChooser_->Filter(path, sources_);

            	if (candidates.size() < sizes_.size()) {
            	    sizes_[candidates.size()] ++;
            	} else {
                    sizes_.resize(candidates.size() + 1, 0);
                    sizes_[candidates.size()] ++;
            	}
				if (candidates.size() == 1) {
					 DEBUG(candidates.size() << " " << g_.int_id(candidates[0].e_) << " Path id :" << path.GetId()<< "  Edge len : " << g_.length(candidates[0].e_))

                     int gap = cfg::get().pe_params.param_set.scaffolder_options.fix_gaps ?
                             gapJoiner_.FixGap(path.Back(), candidates.back().e_, candidates.back().d_) :
                             candidates.back().d_;

					 if (gap != GapJoiner::INVALID_GAP) {
					     DEBUG("Scaffolding. PathId: " << path.GetId() << " path length: " << path.Length() << ", fixed gap length: " << gap);
					     path.PushBack(candidates.back().e_, gap);
					 } else {
					     DEBUG("Looks like wrong scaffolding. PathId: " << path.GetId() << " path length: " << path.Length() << ", fixed gap length: " << candidates.back().d_);
					     break;
					 }
					 //path.Print();
			    }
            }

			if (path.getLoopDetector().IsCycled(maxLoops_)) {
				path.getLoopDetector().RemoveLoop();
				break;
			}
        }
        while (candidates.size() == 1);
    }


public:

    ScaffoldingPathExtender(Graph& g, ExtensionChooser * usualEC, ExtensionChooser * scaffoldingEC): SimplePathExtender(g, usualEC),
            scaffoldingExtensionChooser_(scaffoldingEC),
            gapJoiner_(g, cfg::get().pe_params.param_set.scaffolder_options.min_gap_score,
                    (int) (cfg::get().pe_params.param_set.scaffolder_options.max_must_overlap * g.k()),
                    (int) (cfg::get().pe_params.param_set.scaffolder_options.max_can_overlap * g.k()),
                    cfg::get().pe_params.param_set.scaffolder_options.short_overlap)
                    //cfg::get().pe_paramsparam_set.scaffolder_options.artificial_gap)
    {
        InitSources();
    }

    virtual ~ScaffoldingPathExtender() {

    }

};




class ScaffoldingOnlyPathExtender: public SimplePathExtender {

protected:

    ExtensionChooser * scaffoldingExtensionChooser_;

    std::vector<int> sizes_;

    ExtensionChooser::EdgeContainer sources_;


    void InitSources() {
        sources_.clear();

        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (g_.IncomingEdgeCount(g_.EdgeStart(*iter)) == 0) {
                sources_.push_back(EdgeWithDistance(*iter, 0));
            }
        }
        DEBUG("Found " << sources_.size() << " source edges");
    }

    bool IsSink(EdgeId e)
    {
        return g_.OutgoingEdgeCount(g_.EdgeEnd(e)) == 0;
    }

    virtual void GrowPath(BidirectionalPath& path)
    {
        ExtensionChooser::EdgeContainer candidates;
        do {
            candidates.clear();

            if (IsSink(path.Back())) {
                candidates = scaffoldingExtensionChooser_->Filter(path, sources_);

                if (candidates.size() < sizes_.size()) {
                    sizes_[candidates.size()] ++;
                } else {
                    sizes_.resize(candidates.size() + 1, 0);
                    sizes_[candidates.size()] ++;
                }
                if (candidates.size() == 1) {
                     DEBUG(candidates.size() << " " << g_.int_id(candidates[0].e_) << " Path id :" << path.GetId()<< "  Edge len : " << g_.length(candidates[0].e_))
                     path.PushBack(candidates.back().e_, candidates.back().d_);
                     //path.Print();
                }
            }

            if (path.getLoopDetector().IsCycled(maxLoops_)) {
                path.getLoopDetector().RemoveLoop();
                break;
            }
        }
        while (candidates.size() == 1);
    }


public:

    ScaffoldingOnlyPathExtender(Graph& g_,  ExtensionChooser * scaffoldingEC): SimplePathExtender(g_, 0),
            scaffoldingExtensionChooser_(scaffoldingEC)  {
        InitSources();
    }

    virtual ~ScaffoldingOnlyPathExtender() {

    }

};

}

#endif /* PATH_EXTENDER_HPP_ */
