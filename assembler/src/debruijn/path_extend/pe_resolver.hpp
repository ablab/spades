//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * pe_resolver.hpp
 *
 *  Created on: Mar 12, 2012
 *      Author: andrey
 */

#ifndef PE_RESOLVER_HPP_
#define PE_RESOLVER_HPP_

#include "path_extender.hpp"
#include "pe_io.hpp"

namespace path_extend {


class SimpleOverlapRemover {

protected:
	Graph& g_;

	GraphCoverageMap& coverageMap;

private:


	pair<int, int> ComparePaths(int startPos1, int startPos2,
			BidirectionalPath& path1, BidirectionalPath& path2, size_t maxOverlap){
		int curPos = startPos1;
		int lastPos2 = startPos2;
		int lastPos1 = curPos;
		curPos++;
		size_t diff_len = 0;
		while (curPos < (int)path1.Size()) {
			EdgeId currentEdge = path1[curPos];
			vector<size_t> poses2 = path2.FindAll(currentEdge);
			bool found = false;
			for (size_t pos2 = 0; pos2 < poses2.size(); ++pos2) {
				if ((int) poses2[pos2] > lastPos2) {
					if (diff_len > maxOverlap) {
						return make_pair(-1, -1);
					}
					lastPos2 = poses2[pos2];
					lastPos1 = curPos;
					found = true;
					break;
				}
			}
			if (!found) {
				diff_len += g_.length(currentEdge);
			} else {
				diff_len = 0;
			}
			curPos++;
		}
		return make_pair(lastPos1, lastPos2);
	}


	bool cutOverlaps(BidirectionalPath* path1, int path1_first, int path1_last,
			int size1, BidirectionalPath* path2, int path2_first,
			int path2_last, int size2) {
		if (path1_first == 0 && path1_last == size1 - 1) {
			DEBUG("delete path 1");
			path1->Clear();
		} else if (path2_first == 0 && path2_last == size2 - 1) {
			DEBUG("delete path 2");
			path2->Clear();
		} else if (path2_first == 0 && path1_first == 0) {
			if (size1 < size2) {
				DEBUG("delete begin path 1");
				path1->getConjPath()->PopBack(path1_last + 1);
			} else {
				DEBUG("delete begin path 2");
				path2->getConjPath()->PopBack(path2_last + 1);
			}
		} else if ((path1_last == size1 - 1 && path2_last == size2 - 1)) {
			if (size1 < size2) {
				DEBUG("delete end path 1");
				path1->PopBack(path1_last + 1 - path1_first);
			} else {
				DEBUG("delete end path 2");
				path2->PopBack(path2_last + 1 - path2_first);
			}
		} else if (path2_first == 0) {
			DEBUG("delete path 2 begin");
			path2->getConjPath()->PopBack(path2_last + 1);
		} else if (path2_last == size2 - 1) {
			DEBUG("delete path 2 end");
			path2->PopBack(path1_last + 1 - path1_first);
		} else if (path1_first == 0) {
			DEBUG("delete path1 begin");
			path1->getConjPath()->PopBack(path1_last + 1);
		} else if (path1_last == size1 - 1) {
			path1->PopBack(path1_last + 1 - path1_first);
			DEBUG("delete path1 end")
		} else {
			DEBUG("nothing delete");
			return false;
		}
		return true;
	}

	class PathIdComparator {
		public:

		    bool operator() (const BidirectionalPath* p1, const BidirectionalPath* p2) const {
		        if (p1->GetId() == p2->GetId()){
		        	DEBUG("WRONG ID !!!!!!!!!!!!!!!!!");
		        	return p1->Length() < p2->Length();
		        }
		    	return p1->GetId() < p2->GetId();
		    }
	};
	class EdgeIdComparator {
		Graph& g_;
	public:
		EdgeIdComparator(Graph& g): g_(g){

		}

		bool operator()(const EdgeId& e1, const EdgeId& e2) const {
			if (g_.length(e1) < g_.length(e2)) {
				return true;
			}
			if (g_.length(e2) < g_.length(e1)) {
				return false;
			}
			return e1.int_id() < e2.int_id();
		}
	};


	std::vector<EdgeId> getSortedEdges() {
		std::set<EdgeId> edges_set;
		for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			edges_set.insert(*iter);
			edges_set.insert(g_.conjugate(*iter));
		}
		std::vector<EdgeId> edges(edges_set.begin(), edges_set.end());
		std::sort(edges.begin(), edges.end(), EdgeIdComparator(g_));
		return edges;
	}

public:
	SimpleOverlapRemover(Graph& g, GraphCoverageMap& cm) :
			g_(g), coverageMap(cm) {
	}

	void removeSimilarPaths(PathContainer& paths, size_t maxOverlapLength) {
		std::vector<EdgeId> edges = getSortedEdges();
		for (size_t edgeId = 0; edgeId < edges.size(); ++edgeId) {
			EdgeId edge = edges.at(edgeId);
			std::set<BidirectionalPath *> covPaths = coverageMap.GetCoveringPaths(edge);
			std::vector<BidirectionalPath*> covVector(covPaths.begin(), covPaths.end());
			std::sort(covVector.begin(), covVector.end(), PathIdComparator());
			for (size_t vect_i = 0; vect_i < covVector.size(); ++vect_i) {
				BidirectionalPath* path1 = covVector.at(vect_i);
				if (covPaths.find(path1) == covPaths.end()) {
					continue;
				}
				for (size_t vect_i1 = vect_i + 1; vect_i1 < covVector.size(); ++vect_i1) {
					BidirectionalPath* path2 = covVector.at(vect_i1);
					if (covPaths.find(path2) == covPaths.end()) {
						continue;
					}
					if ((*path1) == (*path2)) {
						path2->Clear();
						covPaths = coverageMap.GetCoveringPaths(edge);
						DEBUG("delete path2 it is equal");
						continue;
					}
					if (g_.length(edge) <= maxOverlapLength) {
						DEBUG("small edge");
						continue;
					}
					vector<size_t> poses1 = path1->FindAll(edge);
					for (size_t index_in_posis = 0; index_in_posis < poses1.size(); ++index_in_posis) {
						vector<size_t> poses2 = path2->FindAll(edge);
						for (size_t index_in_posis2 = 0; index_in_posis2 < poses2.size(); ++index_in_posis2) {
							int path2_last = poses2[index_in_posis2];
							int path1_last = poses1[index_in_posis];
							if (path1_last >= (int)path1->Size() || path2_last >= (int)path2->Size()){
								continue;
							}

							vector<int> other_path_end;
							pair<int, int > posRes = ComparePaths(path1_last, path2_last, *path1, *path2, maxOverlapLength);
							path1_last = posRes.first;
							path2_last = posRes.second;
							if (path1_last == -1 || path2_last == -1){
								continue;
							}
							BidirectionalPath* conj_current = path1->getConjPath();
							BidirectionalPath* conj_second = path2->getConjPath();
							int path1_first = conj_current->Size() - poses1[index_in_posis] - 1;
							int path2_first = conj_second->Size() - poses2[index_in_posis2] - 1;
							posRes = ComparePaths(path1_first, path2_first, *conj_current, *conj_second, maxOverlapLength);
							path1_first = posRes.first;
							path2_first = posRes.second;
							if (path1_first == -1 || path2_first == -1) {
								continue;
							}
							DEBUG("pos " << path1_last << " pos " << path2_last );
							DEBUG("try to delete smth ");
							path1->Print();
							DEBUG("second path");
							path2->Print();
							path2_first = conj_second->Size() - path2_first - 1;
							path1_first = conj_current->Size() - path1_first - 1;
							DEBUG("path1 begin " << path1_first << " path1 end "
									<< path1_last << " path2_begin "
									<< path2_first << " path2_end "
									<< path2_last);
							cutOverlaps(path1, path1_first, path1_last,
									path1->Size(), path2, path2_first,
									path2_last, path2->Size());
							covPaths = coverageMap.GetCoveringPaths(edge);
						}
					}
				}
			}

		}
		DEBUG("END ALL CUT")
	}


	void removeOverlaps(PathContainer& paths) {
		for (size_t i = 0; i < paths.size(); i++) {
			removePathOverlap(paths, paths.Get(i));
			removePathOverlap(paths, paths.GetConjugate(i));
		}
	}

	bool hasAlreadyOverlapedEnd(BidirectionalPath * path){
		return !path->isOverlap() and path->getOverlapedEnd().size() > 0;
	}

	bool hasAlreadyOverlapedBegin(BidirectionalPath * path){
		return !path->isOverlap() and path->getOverlapedBegin().size() > 0;
	}

	bool isSamePath(BidirectionalPath * path1, BidirectionalPath * path2){
		return path2->GetId() == path1->GetId() or path2->GetId() == path1->getConjPath()->GetId();
	}


	void removePathOverlap(PathContainer& pathsContainer, BidirectionalPath* currentPath) {
        int last_index = currentPath->Size() - 1;
		if (last_index <= 0
				or coverageMap.GetCoverage(currentPath->At(last_index)) <= 1
				or hasAlreadyOverlapedEnd(currentPath)){
			return;
		}

		std::set<BidirectionalPath *> paths = coverageMap.GetCoveringPaths(currentPath->At(last_index));

		BidirectionalPath* overlapedPath = NULL;

		size_t overlapedSize = 0;

		for (auto path_iter = paths.begin(); path_iter != paths.end(); ++path_iter) {
			if (isSamePath(*path_iter, currentPath)
					|| hasAlreadyOverlapedBegin(*path_iter)) {
				continue;
			}
			size_t over_size = currentPath->OverlapEndSize(*path_iter);
			if (over_size > overlapedSize) {
				overlapedSize = over_size;
				overlapedPath = *path_iter;
			} else if (over_size == overlapedSize
					&& (overlapedPath == NULL
							|| (*path_iter)->GetId() < overlapedPath->GetId())) {
				overlapedPath = *path_iter;
			}
		}
        if (overlapedPath == NULL){
            return;
        }

		BidirectionalPath* conj_path2 = overlapedPath->getConjPath();
		if (overlapedSize > 0) {
			DEBUG("remove overlaps, change paths ");
			if (currentPath->isOverlap() && overlapedSize == currentPath->Size()) {
				conj_path2->PopBack(overlapedSize);
				DEBUG("change second path");
				overlapedPath->changeOverlapedBeginTo(currentPath);
			}else if (overlapedPath->isOverlap() && overlapedPath->Size() == overlapedSize) {
				currentPath->PopBack(overlapedSize);
				DEBUG("change first path");
				currentPath->changeOverlapedEndTo(overlapedPath);
			}else if (overlapedSize < overlapedPath->Size() && overlapedSize < currentPath->Size()) {
				BidirectionalPath* overlap = new BidirectionalPath(g_, currentPath->Head());
				BidirectionalPath* conj_overlap = new BidirectionalPath(g_, g_.conjugate(currentPath->Head()));
				pathsContainer.AddPair(overlap, conj_overlap);
				currentPath->PopBack();
				conj_path2->PopBack();
				for (size_t i = 1; i < overlapedSize; ++i) {
					conj_overlap->Push(g_.conjugate(currentPath->Head()));
					currentPath->PopBack();
					conj_path2->PopBack();
				}
				coverageMap.Subscribe(overlap);
				overlap->setOverlap(true);
				coverageMap.Subscribe(conj_overlap);
				currentPath->changeOverlapedEndTo(overlap);
				overlapedPath->changeOverlapedBeginTo(overlap);
				DEBUG("add new overlap");

			}
		}
	}

};


// Old one

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


    void JoinByOverlaps(BidirectionalPath * ending, BidirectionalPath * starting, size_t length) {
        for (size_t i = length; i < starting->Size(); ++i) {
            ending->PushBack(starting->At(i), starting->GapAt(i));
        }

        starting->Clear();
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

        for (size_t i = 1; i < maxLen; ++i) {
            removeOverlapsOfLength(paths, i);
            updateMultiCoveredEdges();
        }
    }

};


class PathExtendResolver {

protected:
    Graph& g_;
    size_t k_;

    bool isEdgeSuitable(EdgeId e, double minEdgeCoverage = 0.0) {
//        int len = g.length(e);
//
//        int max_len = (int) cfg::get().K; /*cfg::get().pe_params.param_set.seed_selection.max_len >= 0 ?
//                cfg::get().pe_params.param_set.seed_selection.max_len :
//                -cfg::get().pe_params.param_set.seed_selection.max_len * (int) cfg::get().K;*/

        return math::ge(g_.coverage(e), minEdgeCoverage) &&
                !(cfg::get().pe_params.param_set.seed_selection.exclude_chimeric /*&& len >= cfg::get().pe_params.param_set.seed_selection.min_len && len <= max_len*/);
    }

public:
    PathExtendResolver(Graph& g): g_(g), k_(g.k()) {
    }

    PathContainer ReturnEdges(double minEdgeCoverage){
        std::set<EdgeId> included;
        PathContainer edges;
        size_t count = 0;
        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (included.count(*iter) == 0 && isEdgeSuitable(*iter, minEdgeCoverage)) {
                count ++;
                edges.AddPair(new BidirectionalPath(g_, *iter), new BidirectionalPath(g_, g_.conjugate(*iter)));
                included.insert(*iter);
                included.insert(g_.conjugate(*iter));
            }
        }
        return edges;
    }

    PathContainer returnFilteredEdges(double minEdgeCoverage){
		std::set<EdgeId> included;
		PathContainer edges;
		size_t count = 0;
		for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			if (included.count(*iter) == 0 && isEdgeSuitable(*iter, minEdgeCoverage)) {
				count ++;
				if (!InCycle(*iter) and !InBuble(*iter)) {
					edges.AddPair(new BidirectionalPath(g_, *iter), new BidirectionalPath(g_, g_.conjugate(*iter)));
					included.insert(*iter);
					included.insert(g_.conjugate(*iter));
				}
			}
		}
		return edges;
	}

    bool InCycle(EdgeId e)
    {
    	auto edges = g_.OutgoingEdges(g_.EdgeEnd(e));
    	if (edges.size() >= 1) {
			for (auto it = edges.begin(); it != edges.end();  ++ it) {
				if (g_.EdgeStart(e) == g_.EdgeEnd(*it)) {
				   return true;
				}
			}
    	}
    	return false;
    }

    bool InBuble(EdgeId e){
    	auto edges = g_.OutgoingEdges(g_.EdgeStart(e));
    	auto endVertex = g_.EdgeEnd(e);
    	for (auto it = edges.begin(); it != edges.end(); ++it){
    		if ((g_.EdgeEnd(*it) == endVertex) and (*it != e)){
    			return true;
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

    PathContainer makeSimpleSeeds() {
		std::set<EdgeId> included;
		PathContainer edges;
		for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			if (g_.int_id(*iter) <= 0 or InCycle(*iter)) {
				continue;
			}
			if (included.count(*iter) == 0) {
				edges.AddPair(new BidirectionalPath(g_, *iter),
						new BidirectionalPath(g_, g_.conjugate(*iter)));
				included.insert(*iter);
				included.insert(g_.conjugate(*iter));
			}
		}
		return edges;
	}

    PathContainer extendSeeds(PathContainer& seeds, PathExtender& pathExtender) {
        PathContainer paths;
        pathExtender.GrowAll(seeds, &paths);
        return paths;
    }

    void removeOverlaps(PathContainer& paths, GraphCoverageMap& coverageMap, size_t maxOverlapedSize, ContigWriter& writer, string output_dir) {
        SimpleOverlapRemover overlapRemover(g_, coverageMap);
        paths.SortByLength();
        DEBUG("Removing overlaps");
        paths.SortByLength();
        overlapRemover.removeOverlaps(paths);
        paths.SortByLength();
        DEBUG("remove equals paths");
        paths.SortByLength();
		overlapRemover.removeSimilarPaths(paths, 0);
		paths.SortByLength();
		overlapRemover.removeSimilarPaths(paths, maxOverlapedSize);
		paths.SortByLength();
		DEBUG("end removing");
    }

    void addUncoveredEdges(PathContainer& paths, GraphCoverageMap& coverageMap) {
        std::set<EdgeId> included;
        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (included.count(*iter) == 0 && !coverageMap.IsCovered(*iter)) {
                paths.AddPair(new BidirectionalPath(g_, *iter), new BidirectionalPath(g_, g_.conjugate(*iter)));
                included.insert(*iter);
                included.insert(g_.conjugate(*iter));
            }
        }
    }

};

} /* PE_RESOLVER_HPP_ */

#endif
