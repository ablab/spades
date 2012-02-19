/*
 * pe_resolver.hpp
 *
 *  Created on: Mar 12, 2012
 *      Author: andrey
 */

#ifndef PE_RESOLVER_HPP_
#define PE_RESOLVER_HPP_

#include "path_extender.hpp"

namespace path_extend {


class OverlapRemover {

protected:
    Graph& g;

    GraphCoverageMap& coverageMap;

    std::set<EdgeId> multiCoveredEdges;


    void findPathPairs(BidirectionalPath& sample, std::vector<BidirectionalPath*> * ending, std::vector<BidirectionalPath*> * starting) {
        starting->clear();
        ending->clear();

        if (coverageMap.GetCoverage(sample) <= 1) {
            return;
        }

        auto coveringPaths = coverageMap.GetCoveringPaths(sample);
        for (auto iter = coveringPaths.begin(); iter != coveringPaths.end(); ++iter) {
            if ((*iter)->StartsWith(sample)) {
                starting->push_back(*iter);
            }
            else if ((*iter)->EndsWith(sample)) {
                ending->push_back(*iter);
            }
        }
    }

    void cropOverlaps(std::vector<BidirectionalPath*> * ending, std::vector<BidirectionalPath*> * starting, size_t length) {
        if (ending->size() <= starting->size()) {
//            if (ending->size() > 0) {
//                INFO("Removing overlap of length " << length);
//                INFO("Starting paths")
//                for (auto iter = starting->begin(); iter != starting->end(); ++iter) {
//                    (*iter)->print();
//                }
//                INFO("Ending paths")
//                for (auto iter = ending->begin(); iter != ending->end(); ++iter) {
//                    (*iter)->print();
//                }
//            }

            for (auto iter = ending->begin(); iter != ending->end(); ++iter) {
                for (size_t i = 0; i < length; ++i) {
                    if ((*iter)->Back() != g.conjugate((*iter)->Back())) {
                        (*iter)->PopBack();
                    }
                }
            }
        }
        //Removing overlaps from starting paths will be performed at RC strand
    }


    void removeOverlapsOfLength(PathContainer& paths, size_t length) {
        std::vector<BidirectionalPath*> ending;
        std::vector<BidirectionalPath*> starting;

        for (auto iter = multiCoveredEdges.begin(); iter != multiCoveredEdges.end(); ++iter) {
            auto coveringPaths = coverageMap.GetCoveringPaths(*iter);

            for (auto pathIter = coveringPaths.begin(); pathIter != coveringPaths.end(); ++pathIter) {
                if ((*pathIter)->Size() > length && (*pathIter)->At(length - 1) == *iter) {
                    BidirectionalPath sample = (*pathIter)->SubPath(0, length);
                    findPathPairs(sample, &starting, &ending);
                    cropOverlaps(&starting, &ending, length);
                }

                if ((*pathIter)->Size() > length && (*pathIter)->ReverseAt(length - 1) == *iter) {
                    BidirectionalPath sample = (*pathIter)->SubPath((*pathIter)->Size() - length);
                    findPathPairs(sample, &starting, &ending);
                    cropOverlaps(&starting, &ending, length);
                }
            }
        }
    }

    void updateMultiCoveredEdges() {
        for (auto iter = multiCoveredEdges.begin(); iter != multiCoveredEdges.end(); ) {
            auto next = iter;
            ++next;
            if (coverageMap.GetCoverage(*iter) <= 1) {
                multiCoveredEdges.erase(iter);
            }
            iter = next;
        }
    }


public:
    OverlapRemover(Graph& g_, GraphCoverageMap& cm): g(g_), coverageMap(cm), multiCoveredEdges() {
        for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (coverageMap.GetCoverage(*iter) > 1) {
                multiCoveredEdges.insert(*iter);
            }
        }
    }

    void removeOverlaps(PathContainer& paths) {
        if (paths.size() == 0) {
            return;
        }

        size_t maxLen = paths.Get(0)->Size();
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if (iter.get()->Size() > maxLen) {
                maxLen = iter.get()->Size();
            }
        }

        removeOverlapsOfLength(paths, 1);
        updateMultiCoveredEdges();
//        for (size_t i = 1; i < maxLen; ++i) {
//        }
    }

};


class PathExtendResolver {

protected:
    Graph& g;
    size_t k_;

    bool isEdgeSuitable(EdgeId e, double minEdgeCoverage = 0.0) {
        size_t len = g.length(e);
        size_t delta = params.param_set.seed_selection.chimeric_delta;
        delta = delta > k_ ? k_ : delta;

        return math::ge(g.coverage(e), minEdgeCoverage) &&
                (!params.param_set.seed_selection.exclude_chimeric
                        || (len < k_ - delta
                        || len > k_ + params.param_set.seed_selection.chimeric_delta));
    }

public:
    PathExtendResolver(Graph& g_, size_t k): g(g_), k_(k) {
    }

    PathContainer ReturnEdges(double minEdgeCoverage){
        std::set<EdgeId> included;
        PathContainer edges;
        size_t count = 0;
        for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (included.count(*iter) == 0 && isEdgeSuitable(*iter, minEdgeCoverage)) {
                count ++;
                edges.AddPair(new BidirectionalPath(g, *iter), new BidirectionalPath(g, g.conjugate(*iter)));
                included.insert(*iter);
                included.insert(g.conjugate(*iter));
            }
        }
        return edges;
    }

    PathContainer returnFilteredEdges(double minEdgeCoverage){
		std::set<EdgeId> included;
		PathContainer edges;
		size_t count = 0;
		for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			if (included.count(*iter) == 0 && isEdgeSuitable(*iter, minEdgeCoverage)) {
				count ++;
				if (!InCycle(*iter)) {
					edges.AddPair(new BidirectionalPath(g, *iter), new BidirectionalPath(g, g.conjugate(*iter)));
					included.insert(*iter);
					included.insert(g.conjugate(*iter));
				}
			}
		}
		return edges;
	}

    bool InCycle(EdgeId e)
    {
    	auto edges = g.OutgoingEdges(g.EdgeEnd(e));
    	if (edges.size() > 1) {
			for (auto it = edges.begin(); it != edges.end();  ++ it) {
				if (g.EdgeStart(e) == g.EdgeEnd(*it)) {
				   return true;
				}
			}
    	}
    	return false;
    }

    PathContainer makeSeeds(PathExtender& seedExtender, double minEdgeCoverage = 0.0) {
    	PathContainer edges = returnFilteredEdges(minEdgeCoverage);
        PathContainer seeds;
        seedExtender.GrowAll(edges, &seeds);
        return seeds;
    }

    PathContainer extendSeeds(PathContainer& seeds, PathExtender& pathExtender) {
        PathContainer paths;
        pathExtender.GrowAll(seeds, &paths);
        return paths;
    }

    void removeOverlaps(PathContainer& paths, GraphCoverageMap& coverageMap) {
        OverlapRemover overlapRemover(g, coverageMap);
        overlapRemover.removeOverlaps(paths);
    }

    void addUncoveredEdges(PathContainer& paths, GraphCoverageMap& coverageMap) {
        std::set<EdgeId> included;
        for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (included.count(*iter) == 0 && !coverageMap.IsCovered(*iter)) {
                paths.AddPair(new BidirectionalPath(g, *iter), new BidirectionalPath(g, g.conjugate(*iter)));
                included.insert(*iter);
                included.insert(g.conjugate(*iter));
            }
        }
    }

};

} /* PE_RESOLVER_HPP_ */

#endif
