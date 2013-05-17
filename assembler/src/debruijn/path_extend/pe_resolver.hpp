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

namespace path_extend {


class SimpleOverlapRemover {

protected:
	Graph& g;

	GraphCoverageMap& coverageMap;

private:
	void ComparePaths(int startPosition1, int startPosition2, int& other_path_last_pos, int& this_path_last_pos, BidirectionalPath& first, BidirectionalPath& second, vector<int>& result){
		int end_position = startPosition1;
		other_path_last_pos = startPosition2;
		this_path_last_pos = end_position;
		result.push_back(startPosition2);
		end_position++;
		while (end_position < (int)first.Size() && other_path_last_pos < (int)second.Size()) {
			EdgeId currentEdge = first[end_position];
			vector<size_t> other_path_poses = second.FindAll(currentEdge);
			bool found = false;
			for (size_t other_path_pose = 0; other_path_pose < other_path_poses.size();++other_path_pose) {
				if ((int)other_path_poses[other_path_pose] > other_path_last_pos) {
					result.push_back(other_path_poses[other_path_pose]);
					other_path_last_pos = other_path_poses[other_path_pose];
					this_path_last_pos = end_position;
					found = true;
					break;
				}
			}
			if (!found) {
				result.push_back(-1);
			}
			end_position++;
		}
	}

	bool IsPathDifferent(BidirectionalPath* currentPath, vector<int> second_path_poses, int this_path_first_pos, int this_path_last_pos, int maxOverlapLength) {
		int diffLen = 0;
		for (int i = this_path_first_pos; i <= this_path_last_pos; ++i) {
			if (second_path_poses[i] == -1) {
				diffLen = g.length(currentPath->At(i));
				if (diffLen > maxOverlapLength) {
					return true;
				}
			} else {
				diffLen = 0;
			}
		}
		return false;
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


public:
	SimpleOverlapRemover(Graph& g_, GraphCoverageMap& cm) :
			g(g_), coverageMap(cm) {
	}

	void removeLastOverlaps(PathContainer& paths, double maxOverlapLength) {
		for (auto edge = g.SmartEdgeBegin(); !edge.IsEnd(); ++edge) {
			if (coverageMap.GetCoverage(*edge) > 1) {
				bool change = true;
				std::set<BidirectionalPath *> coveredPathes = coverageMap.GetCoveringPaths(*edge);
				while (coveredPathes.size() > 1 and change) {
					change = false;
					BidirectionalPath* currentPath = *(coveredPathes.begin());
					BidirectionalPath* secondPath = *(++coveredPathes.begin());
					DEBUG("delete similar paths");
					currentPath->Print();
					DEBUG("second path");
					secondPath->Print();
					if ((*currentPath) == (*secondPath)){
						secondPath->Clear();
						change = true;
						coveredPathes = coverageMap.GetCoveringPaths(*edge);
						DEBUG("paths equals, delete second");
						continue;
					}
					if (g.length(*edge) <= maxOverlapLength) {
						 continue;
					}
					vector<size_t> positions = currentPath->FindAll(*edge);
					for (size_t index_in_posis = 0; index_in_posis < positions.size(); ++index_in_posis) {
						int path2_last = secondPath->FindFirst(*edge);
						int path1_last = positions[index_in_posis];
						vector<int> other_path_end;
						ComparePaths(path1_last, path2_last, path2_last, path1_last, *currentPath, *secondPath, other_path_end);
						BidirectionalPath* conj_current = paths.FindConjugate(currentPath);
						BidirectionalPath* conj_second = paths.FindConjugate(secondPath);
						int path1_first = conj_current->Size() - positions[index_in_posis] - 1;
						int path2_first  = secondPath->Size() - secondPath->FindFirst(*edge) - 1;
						vector<int> other_path_start_temp;
						ComparePaths(path1_first, path2_first, path2_first, path1_first, *conj_current, *conj_second, other_path_start_temp);
						vector<int> other_path_start;
						for(size_t i = 0; i < other_path_start_temp.size(); ++i){
							other_path_start.push_back( other_path_start_temp[i] < 0 ?  other_path_start_temp[i] : conj_second->Size() - other_path_start_temp[i] - 1);
						}
						path2_first = conj_second->Size() - path2_first - 1;
						path1_first =  conj_current->Size() - 1 - path1_first ;
						std::reverse(other_path_start.begin(), other_path_start.end());
						other_path_start.pop_back();
						other_path_start.insert(other_path_start.end(), other_path_end.begin(), other_path_end.end());
						DEBUG("path1 begin "<< path1_first<< " path1 end "<<path1_last << " path2_begin " << path2_first << " path2_end " << path2_last);
						bool different = IsPathDifferent(currentPath, other_path_start, path1_first, path1_last, maxOverlapLength);
						if (different){
							DEBUG("different");
							continue;
						}
						change = cutOverlaps(currentPath, path1_first, path1_last, currentPath->Size(),
								secondPath, path2_first, path2_last, secondPath->Size());
						coveredPathes = coverageMap.GetCoveringPaths(*edge);
					}
				}
			}
		}
	}

	void removeOverlaps(PathContainer& paths) {
		for (size_t i = 0; i < paths.size(); i++) {
			removeEndOverlapedSomeBegin(paths, paths.Get(i), paths.GetConjugate(i));
			removeEndOverlapedSomeBegin(paths, paths.GetConjugate(i), paths.Get(i));
		}
	}

	void removeEndOverlapedSomeBegin(PathContainer& pathsContainer, BidirectionalPath * currentPath, BidirectionalPath * currentConjPath) {
        int last_index = currentPath->Size() - 1;
		if (last_index <= 0 or coverageMap.GetCoverage(currentPath->At(last_index)) <= 1 or (!currentPath->isOverlap() and currentPath->getOverlapedEnd().size() > 0)){
			return;
		}
		std::set<BidirectionalPath *> paths = coverageMap.GetCoveringPaths(currentPath->At(last_index));
        size_t max_overlaped_size = 0;
		BidirectionalPath* max_overlaped_path = NULL;
		for (auto path_iter = paths.begin(); path_iter != paths.end(); ++path_iter) {
			BidirectionalPath* path = *path_iter;
			if (path->GetId() == currentPath->GetId() or (path->GetId() == currentConjPath->GetId()) or (!path->isOverlap() and path->getOverlapedBegin().size() > 0)) {
                continue;
			}
			size_t  overlaped_size = currentPath->OverlapEndSize(path);
            if (overlaped_size > max_overlaped_size) {
				max_overlaped_size = overlaped_size;
				max_overlaped_path = path;
			} else if (overlaped_size == max_overlaped_size && (max_overlaped_path == NULL or path->GetId() < max_overlaped_path->GetId())){
				max_overlaped_path = path;
			}
		}
        if (max_overlaped_path == NULL){
            return;
        }
		BidirectionalPath* conj_second_path = pathsContainer.FindConjugate(max_overlaped_path);
		if (max_overlaped_size > 0) {
			DEBUG("remove overlaps, change paths ");
//			DEBUG("before first");
//			currentPath->Print();
//			DEBUG("before second");
//			max_overlaped_path->Print();
			if (currentPath->isOverlap() && max_overlaped_size == currentPath->Size()) {
				conj_second_path->PopBack(max_overlaped_size);
				DEBUG("change second path");
				//max_overlaped_path->Print();
				max_overlaped_path->changeOverlapedBeginTo(currentPath);
			}else if (max_overlaped_path->isOverlap() && max_overlaped_path->Size() == max_overlaped_size) {
				currentPath->PopBack(max_overlaped_size);
				DEBUG("change first path");
				//currentPath->Print();
				currentPath->changeOverlapedEndTo(max_overlaped_path);
			}else if (max_overlaped_size < max_overlaped_path->Size() && max_overlaped_size < currentPath->Size()) {
				BidirectionalPath* overlap = new BidirectionalPath(g, currentPath->Head());
				BidirectionalPath* conj_overlap = new BidirectionalPath(g, g.conjugate(currentPath->Head()));
				pathsContainer.AddPair(overlap, conj_overlap);
				currentPath->PopBack();
				conj_second_path->PopBack();
				for (size_t i = 1; i < max_overlaped_size; ++i) {
					conj_overlap->Push(g.conjugate(currentPath->Head()));
					currentPath->PopBack();
					conj_second_path->PopBack();
				}
				coverageMap.Subscribe(overlap);
				overlap->setOverlap(true);
				coverageMap.Subscribe(conj_overlap);
				currentPath->changeOverlapedEndTo(overlap);
				max_overlaped_path->changeOverlapedBeginTo(overlap);
				DEBUG("add new overlap");
//				DEBUG("now first");
//				currentPath->Print();
//				DEBUG("now second");
//				max_overlaped_path->Print();
//				DEBUG("overlap");
//				overlap->Print();

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

    void removeOverlaps(PathContainer& paths, GraphCoverageMap& coverageMap, size_t maxOverlapedSize) {
        SimpleOverlapRemover overlapRemover(g_, coverageMap);
        INFO("Removing overlaps");
        overlapRemover.removeOverlaps(paths);
        DEBUG("remove equals paths");
		overlapRemover.removeLastOverlaps(paths, maxOverlapedSize);
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
