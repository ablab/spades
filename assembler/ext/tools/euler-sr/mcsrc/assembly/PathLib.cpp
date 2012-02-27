/***************************************************************************
 * Title:          PathLib.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/21/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "PathLib.h"
#include "PathBranch.h"
#include <sstream>


ssize_t FindPathCoverage(PathInterval *path, ssize_t pathLength,
										 PathBranch &pathTree) {
	ssize_t pi;
	PathBranch::BranchMap::iterator branchIt;
	PathBranch *curBranch = &pathTree;
	ssize_t edgeIndex;
	for (pi = 0; pi < pathLength; pi++) {
		edgeIndex = path[pi].edge;
		// REMOVE THIS WHEN PATHS ARE CONDENSED!
		if (pi > 0 and edgeIndex == path[pi-1].edge)
			continue;
		branchIt  = curBranch->branches.find(edgeIndex);
		if (branchIt == curBranch->branches.end())
			return -1;
		else
			curBranch = branchIt->second;
	}
	return curBranch->count;
}

void CollectPathTreeOnPath(PathInterval* path, ssize_t pathLength,
													 PathBranch &pathTree) {
	PathBranch *curBranch, *nextBranch;
	curBranch = &pathTree;
	PathBranch::BranchMap::iterator branchIt;
	ssize_t pi, edgeIndex;
	ssize_t level = 0;
	for (pi = 0; pi < pathLength; pi++ ){ 
		while (pi < pathLength and path[pi].edge == -1) {
			++pi;
		}
		if (pi == pathLength)
			break;

		edgeIndex = path[pi].edge;
			

		// This path is from a multi-segment read on the same edge, don't bother combinging it here.

		
		// REMOVE THIS WHEN MULTI-SEGMENT PATHS 
		// Advance while the intervals are on the same edge.
		while (pi < pathLength - 1 and path[pi].edge == path[pi+1].edge)
			pi++;


		branchIt  = curBranch->branches.find(edgeIndex);
		// each iteration goes down a level in the tree.
		if (branchIt == curBranch->branches.end()) {
			//			cout << "created " << level << " " << edgeIndex << endl;
			nextBranch = new PathBranch;
			curBranch->branches[edgeIndex] = nextBranch;
			curBranch = nextBranch;
			curBranch->count = 0;
		}
		else {
			//			cout << "contains: " << level << " " << edgeIndex << endl;
			curBranch = branchIt->second;
		}
		++ level;
	}
	if (pathLength > 0) 
		curBranch->count++;
}

void CollectPathTree(PathIntervalList &paths,
										 PathLengthList &pathLengths,
										 PathBranch &pathTree) {
	size_t p;
	for (p = 0; p < paths.size(); p++) { 
		CollectPathTreeOnPath(paths[p], pathLengths[p], pathTree);
	}
}

void PrintPathTree(PathBranch &pathTree) {
	std::list<ssize_t> curPath;
	PrintPathTree(pathTree, curPath);
}

void PrintPathTree(PathBranch &pathTree, std::list<ssize_t> &curPath) {
	// print the path tree in a DFS manner.
	PathBranch::BranchMap::iterator branchIt;
	std::list<ssize_t>::iterator listIt;
	if (pathTree.count > 0 and pathTree.branches.size() > 0) {
		for (listIt = curPath.begin(); listIt != curPath.end(); ++listIt) {
			std::cout << *listIt << " ";
		}
		cout << " (" << pathTree.count << ")" << endl;
	}		
	for (branchIt = pathTree.branches.begin();
			 branchIt != pathTree.branches.end(); 
			 ++branchIt) {
		//		std::cout << branchIt->first << " " << flush;
		curPath.push_back(branchIt->first);
		if (branchIt->second->branches.size() == 0) {
			for (listIt = curPath.begin(); listIt != curPath.end(); ++listIt) {
				std::cout << *listIt << " ";
			}
			std::cout << " (";
			std::cout << branchIt->second->count << ")" << std::endl;
		}
		else 
			PrintPathTree(*branchIt->second, curPath);
		curPath.pop_back();
	}
}

void PrintPath(PathTrace &pathTrace) {
	size_t e;
	for (e = 0; e < (*pathTrace.edges).size(); e++) {
		std::cout << (*pathTrace.edges)[e] << " ";
	}
	std::cout << std::endl;
}
	
void PrintPathTraceList(PathTraceList &pathTraces) {
	size_t p;
	for (p = 0; p < pathTraces.size(); p++) {
		PrintPath(pathTraces[p]);
	}
}

void PathTreeToPathList(PathBranch &pathTree, PathTraceList &pathTraces) {
	std::list<ssize_t> curPath;
	PathTreeToPathList(pathTree, curPath, pathTraces);
}



void PathTreeToPathList(PathBranch &pathTree, std::list<ssize_t> &curPath, 
												PathTraceList &pathTraces) {
	PathBranch::BranchMap::iterator branchIt;
	for (branchIt = pathTree.branches.begin();
			 branchIt != pathTree.branches.end(); 
			 ++branchIt) {
		//		std::cout << branchIt->first << " ";
		curPath.push_back(branchIt->first);
		if (branchIt->second->branches.size() == 0) {
			ssize_t curTrace = pathTraces.size();
			pathTraces.push_back(PathTrace());
			std::list<ssize_t>::iterator listIt;
			assert(curPath.size() > 0);
			for (listIt = curPath.begin(); listIt != curPath.end(); ++listIt) {
				pathTraces[curTrace].edges->push_back(*listIt);
			}
			pathTraces[curTrace].count = branchIt->second->count;
			pathTraces[curTrace].parent = curTrace;
		}
		else 
			PathTreeToPathList(*branchIt->second, curPath, pathTraces);
		curPath.pop_back();
	}
} 



void RemoveLowCountPaths(IntervalGraph    &g,
												 PathIntervalList &paths,
												 PathLengthList   &pathLengths,
												 PathBranch       &pathTree, ssize_t minCount) {
	ssize_t p;
	ssize_t numPaths = paths.size();
	ssize_t pathCoverage;
	for (p = 0; p < numPaths; p++) {
		if (pathLengths[p] > 1) {
			pathCoverage = FindPathCoverage(paths[p], pathLengths[p], pathTree);
			if ( pathCoverage < minCount) {
				// Eventually remove the interval from the graph.
				ssize_t pi;
				for (pi = 0; pi < pathLengths[p]; pi++) {
					g.edges[paths[p][pi].edge].MarkIntervalForRemoval(paths[p][pi].index);
				}
				// 
				delete[] paths[p];
				pathLengths[p] = 0;
			}
		}
	}
	g.RemoveMarkedIntervalsNoPaths();
	
}

ssize_t TrimLowCoverageBranches(PathBranch &pathTree, ssize_t minCount) {
	// If this branch is a terminal one.
	if (pathTree.branches.size() == 0) {
		return pathTree.count >= minCount;
	}

	// Ohterwise, check to see if the branches from here work.
	PathBranch::BranchMap::iterator branchIt, deletedIt;
	branchIt = pathTree.branches.begin();
	ssize_t branchIsHighCoverage;
	ssize_t branchContainsHighCoverage = 0;
	while (branchIt != pathTree.branches.end()) { 
		//		curPath.push_back((*branchIt).first);
		branchIsHighCoverage = TrimLowCoverageBranches(*((*branchIt).second), minCount);
		if (!branchIsHighCoverage) {
			delete (*branchIt).second;
			deletedIt = branchIt;
			++branchIt;
			pathTree.branches.erase(deletedIt);
			/*			std::cout << "deleting path: ";
				std::list<ssize_t>::iterator listIt, endList;
				endList = curPath.end();
				for (listIt = curPath.begin(); listIt != endList; ++listIt) {
				std::cout << *listIt << " ";
				}
				std::cout << std::endl;
			*/
		}
		else {
			++branchIt;
			branchContainsHighCoverage = 1;
		}
		//		curPath.pop_back();
	}
	if (branchContainsHighCoverage || pathTree.count >= minCount)
		return 1;
	else
		return 0;
}

void RemoveLowCountTraces(PathTraceList &pathTraces, ssize_t minCount, 
													PathBranch &removedPathTree) {
	size_t p, curPath;
	ssize_t numUnder = 0;
	for (p = 0; p < pathTraces.size(); p++) { 
		if (pathTraces[p].count < minCount)
			++numUnder;
	}
	//	std::cout << "the number under count: " << numUnder << std::endl;
	p = 0; curPath = 0;
	ssize_t numPathTraces= pathTraces.size();
	while ((numPathTraces > p ) and (pathTraces[p].count >= minCount)) {
		p++; curPath++;
	}
	 
	for (;p < pathTraces.size(); p++) {
		if (pathTraces[p].count >= minCount) {
			pathTraces[curPath] = pathTraces[p];
			++curPath;
		}
		else {
			// add this path to a tree
			/*			PathInterval *removedPath = new PathInterval[pathTraces[p].edges->size()];
				ssize_t pi;
				for (pi = 0; pi < (*pathTraces[p].edges).size(); pi++) {
				removedPath[pi].edge = (*pathTraces[p].edges)[pi];
				}

				CollectPathTreeOnPath(removedPath, (*pathTraces[p].edges).size(), removedPathTree);
			*/
			//			delete pathTraces[p].edges;
			//			delete[] removedPath;
		}
	}
	pathTraces.resize(curPath);
}


ssize_t LocatePathStart(PathTraceList &pathTraces, ssize_t startEdge) {
	ssize_t low= 0;
	ssize_t high = pathTraces.size();
	ssize_t cur = (high + low) / 2;

	while (cur >= low and cur < high) {
		assert(pathTraces[cur].edges->size() > 0);
		ssize_t curEdge = (*pathTraces[cur].edges)[0];
		if ( curEdge == startEdge)
			break;

		else if (curEdge < startEdge ){
			low = cur + 1;
		}
		else {
			high = cur;
		}
		cur = (high + low) / 2;
	}
	if (cur < high) {
		// Advance to the first
		while (cur > 0 and (*pathTraces[cur-1].edges)[0] == startEdge) 
			cur--;
		return cur;
	}
	else {
		return -1;
	}
}

void MarkEnclosedPaths(PathTraceList &pathTraces) {
	ssize_t p;
	ssize_t s;
	//UNUSED// ssize_t pathContained;
	for (p = 0; p < pathTraces.size(); p++ ){
		// Look to see if path 'p' includes any subpaths.
		if (pathTraces[p].edges->size() >= 2) {
			ssize_t pe;
			//			std::cout << "checking path: " << p << std::endl;
			for (pe = 1; pe < pathTraces[p].edges->size(); pe++) {
				ssize_t internalEdge = (*pathTraces[p].edges)[pe];
				ssize_t internalEdgeIndex = LocatePathStart(pathTraces, internalEdge);
				if (internalEdgeIndex >= 0) {
					// make sure this path may be contained by the cur path
					while (internalEdgeIndex < pathTraces.size()) {

						// Make sure this path starts on the same as the internal edge
						if ((*pathTraces[internalEdgeIndex].edges)[0] != internalEdge)
							break;

						// This path is already deleted, ignore it.
						if (pathTraces[internalEdgeIndex].parent != internalEdgeIndex) {
							++internalEdgeIndex;
							continue;
						}
							
						// If this path is smaller than the internal edge, it can't 
						// contain it.
						if (pathTraces[internalEdgeIndex].edges->size() >
								pathTraces[p].edges->size() - pe ) {
							++internalEdgeIndex;
							continue;
						}

						assert(pathTraces[internalEdgeIndex].edges->size() > 0);
						
						// Look to see if the rest of the path contains internalEdgeIndex as a subpath
						ssize_t pathsMatch = 1;
						for (s = 0; s < pathTraces[internalEdgeIndex].edges->size(); s++) {
							if ((*pathTraces[internalEdgeIndex].edges)[s] !=
									(*pathTraces[p].edges)[s + pe]) {
								pathsMatch = 0;
								break;
							}
						}
						if (pathsMatch) {
							/*							std::cout << "path: " << p << " includes path: " << internalEdgeIndex << std::endl;
								ssize_t i;
								for (i = 0 ; i < pathTraces[p].edges->size(); i++) 
								std::cout << (*pathTraces[p].edges)[i] <<" ";
								std::cout << std::endl;
								for (i = 0 ; i < pathTraces[internalEdgeIndex].edges->size(); i++)
								std::cout << (*pathTraces[internalEdgeIndex].edges)[i] <<" ";
								std::cout << std::endl;
							*/
							pathTraces[internalEdgeIndex].parent = p;
						}

						// Advance to the next path that possibly starts
						// witn edge 'internalEdge'.
						internalEdgeIndex++;
					}
				} // done looking to see if any other paths start with this edge.
			} // done searching subpaths of this path.
		}
	}
}

void RemoveEnclosedPaths(PathTraceList &pathTraces) {
	ssize_t curTraceIndex = 0;
	ssize_t p;
	p = 0;
	// advance to the first bad case
	if (pathTraces.size() == 0)
		return;

	while (p < pathTraces.size() and pathTraces[p].parent == p) p++;
	curTraceIndex = p;
	for (; p < pathTraces.size(); p++ ) {
		if (pathTraces[p].parent == p) {
			/*			std::cout  << curTraceIndex << " <-- " << p << std::endl;*/
			pathTraces[curTraceIndex] = pathTraces[p];
			pathTraces[p].edges = NULL;
			curTraceIndex++;
		}
		else {
			/*			std::cout << "skipping: " << p << ": ";
				ssize_t i;
				for (i = 0; i < pathTraces[p].edges->size(); ++i) {
				std::cout << (*pathTraces[p].edges)[i] << " ";
				}
				std::cout << std::endl;
			*/
		}
	}
	//	std::cout << "the last path is at: " << (*pathTraces[curTraceIndex-1].edges)[0] << std::endl;
	pathTraces.resize(curTraceIndex);
}


void MarkExInternal(PathTraceList &pathTraces, 
										std::vector<ssize_t> &startList,
										std::vector<ssize_t> &endList,
										std::vector<ssize_t> &intList) {
	ssize_t p = 0; 
	ssize_t pi;
	ssize_t traceLength;
	for (p = 0; p < pathTraces.size(); p++ ){ 
		traceLength = pathTraces[p].edges->size();
		if (traceLength > 1) {
			// Mark the boundaries as extList (unresovled)
			startList[(*pathTraces[p].edges)[0]]++;
			endList[(*pathTraces[p].edges)[traceLength - 1]]++;
			
			// Mark the int part as intList (resolved).
			for (pi = 1; pi < traceLength - 1; pi++ ) {
				intList[(*pathTraces[p].edges)[pi]]++;
			}
		}
	}
}

void PrintPathTraceResolution(TEdgeList &edges,
															PathTraceList &pathTraces, 
															TraceMapMatrix &traceMaps) {
	std::cout << "path resolution: " << std::endl;
	ssize_t p;
	for (p = 0; p < pathTraces.size(); p++ ){ 
		ssize_t pi;
		ssize_t traceLength = (*pathTraces[p].edges).size();
		//UNUSED// int firstEdge = (*pathTraces[p].edges)[0];
		//UNUSED// int lastEdge  = (*pathTraces[p].edges)[traceLength-1];
		std::cout << "path: " << p << " count: " << pathTraces[p].count << " ";
		for (pi = 0; pi < traceLength; pi++ ) { 
			std::cout << (*pathTraces[p].edges)[pi] << " ("
								<< traceMaps[(*pathTraces[p].edges)[pi]].numStart << " "
								<< traceMaps[(*pathTraces[p].edges)[pi]].numInternal << " "
								<< traceMaps[(*pathTraces[p].edges)[pi]].numEnd << " "
								<< traceMaps[(*pathTraces[p].edges)[pi]].inResolved << " "
								<< traceMaps[(*pathTraces[p].edges)[pi]].outResolved << " "
								<< traceMaps[(*pathTraces[p].edges)[pi]].inAdjacent << " "
								<< traceMaps[(*pathTraces[p].edges)[pi]].outAdjacent << " "
								<< edges[(*pathTraces[p].edges)[pi]].length << ") ";
		}
		std::cout << std::endl;
	}
}


void PrintPathTraces(PathTraceList &pathTraces,
										 IntervalGraph &graph,
										 std::ostream &out) {
	
	ssize_t p;
	ssize_t pi;
	for (p = 0; p < pathTraces.size(); p++ ){
		out << "trace: " << p << " ";
		for (pi = 0; pi < pathTraces[p].edges->size(); pi++ ) {
			out << (*pathTraces[p].edges)[pi] << " " << " [" << graph.edges[(*pathTraces[p].edges)[pi]].index << "] "; 
		}
		out << "(" << pathTraces[p].count << ")";
		out << std::endl;
	}
}
	
void PrintEdgeTraces(PathTraceList &pathTraces, TraceMapMatrix &traceMaps, ssize_t edge) {
	ssize_t t;
	for (t = 0 ; t < traceMaps[edge].traces.size(); t++ ) {
		ssize_t e;
		ssize_t path = traceMaps[edge].traces[t].trace;
		for (e = 0; e < (*pathTraces[path].edges).size(); e++ ){ 
			cout << (*pathTraces[path].edges)[e] << " ";
		}
		cout << ", " << pathTraces[path].count << endl;
	}
}

void StoreTraceMaps(PathTraceList &pathTraces, 
										TraceMapMatrix &traceMaps) {
	ssize_t p, pi;
	for (p = 0; p < pathTraces.size(); p++ ) {
		for (pi = 0 ; pi < pathTraces[p].edges->size(); pi++ ){
			traceMaps[(*pathTraces[p].edges)[pi]].traces.push_back(TraceMap(p, pi));
			if (pi == 0) {
				traceMaps[(*pathTraces[p].edges)[pi]].numStart++;
				traceMaps[(*pathTraces[p].edges)[pi]].startEdge = (*pathTraces[p].edges)[0];
			}
			else if (pi > 0 and pi < pathTraces[p].edges->size() - 1)
				traceMaps[(*pathTraces[p].edges)[pi]].numInternal++;
			else {
				traceMaps[(*pathTraces[p].edges)[pi]].numEnd++;
				traceMaps[(*pathTraces[p].edges)[pi]].endEdge = (*pathTraces[p].edges)[pi];
			}
		}
	}
}


void CountAdjacencies(PathTraceList &pathTraces, TraceMapMatrix &traceMap) {
	//UNUSED// ssize_t p, pi;
	//UNUSED// ssize_t numPathTraces = pathTraces.size();
	ssize_t e, t;
	ssize_t path, pos;
	for (e = 0; e < traceMap.size(); e++) {
		std::set<ssize_t> in, out;
		for (t = 0; t < traceMap[e].traces.size(); t++) {
			path = traceMap[e].traces[t].trace;
			pos  = traceMap[e].traces[t].pos;
			// skip on error code.
			if (path < 0)
				continue;

			if (pos > 0) {
				in.insert((*pathTraces[path].edges)[pos-1]);
			}
			if (pos < pathTraces[path].edges->size() - 1) {
				out.insert((*pathTraces[path].edges)[pos+1]);
			}
		}
		traceMap[e].inAdjacent = in.size();
		traceMap[e].outAdjacent = out.size();
	}
}

ssize_t CheckForwardConsistency(PathTrace& path1, ssize_t pos1, PathTrace &path2, ssize_t pos2,
														ssize_t minMatchLength) {
	ssize_t end1 = path1.edges->size();
	ssize_t end2 = path2.edges->size();
	//UNUSED// ssize_t match = 1;
	ssize_t start1 = pos1;
	//UNUSED// ssize_t start2 = pos2;
	while (pos1 < end1 and pos2 < end2) {
		if ((*path1.edges)[pos1] != (*path2.edges)[pos2]) {
			return 0;
		}
		++pos1;
		++pos2;
	}
	if (pos1 - start1 >= minMatchLength)
		return 1;
	return 0; // TODO: check
}

ssize_t CheckReverseConsistency(PathTrace &path1, ssize_t pos1, PathTrace &path2, ssize_t pos2,
														ssize_t minMatchLength) {
	ssize_t start1 = pos1;
	//UNUSED// ssize_t start2 = pos2;
	while(pos1 >= 0 and pos2 >= 0) {
		if ((*path1.edges)[pos1] != (*path2.edges)[pos2])
			return 0;
		--pos1;
		--pos2;
	}
	if (start1 - pos1 >= minMatchLength)
		return 1;
	return 0; // TODO: check
}

// Count the paths that are equal for the first N positions
// then consistent afterwards

ssize_t CountForwardConsistent(PathTrace &refPath, ssize_t refPos, ssize_t matchLength,
													 PathTraceList &pathTraces,
													 TraceMapMatrix &traceMap,
													 ssize_t edgeIndex, ssize_t &forwardConsistentTrace) {
	// some sanity checks
	assert(edgeIndex < traceMap.size());
	assert(refPos < refPath.edges->size());
	ssize_t traceIndex;
	ssize_t trace, tracePos;
	ssize_t numConsistentTraces = 0;
	for (traceIndex = 0; traceIndex < traceMap[edgeIndex].traces.size(); traceIndex++) { 
		trace = traceMap[edgeIndex].traces[traceIndex].trace;
		tracePos = traceMap[edgeIndex].traces[traceIndex].pos;
		if (tracePos == 0) {
			if (CheckForwardConsistency(refPath, refPos, pathTraces[trace], tracePos)) {
				numConsistentTraces++;
				forwardConsistentTrace = trace;
			}
		}
	}
	return numConsistentTraces;
}

ssize_t CountReverseConsistent(PathTrace &refPath, ssize_t refPos,
													 PathTraceList &pathTraces,
													 TraceMapMatrix &traceMap,
													 ssize_t edgeIndex, ssize_t &reverseConsistentTrace) {
	assert(edgeIndex < traceMap.size());
	assert(refPos < refPath.edges->size());
	ssize_t traceIndex;
	ssize_t trace, tracePos;
	ssize_t numConsistentTraces = 0;
	for (traceIndex = 0; traceIndex < traceMap[edgeIndex].traces.size(); traceIndex++) { 
		trace = traceMap[edgeIndex].traces[traceIndex].trace;
		tracePos = traceMap[edgeIndex].traces[traceIndex].pos;
		if (tracePos == pathTraces[trace].edges->size() - 1) {
			if (CheckReverseConsistency(refPath, refPos, pathTraces[trace], tracePos)) {
				numConsistentTraces++;
				reverseConsistentTrace = trace;
			}
		}
	}
	return numConsistentTraces;
}


ssize_t IsEdgePathContained(PathTraceList &pathTraces,
												TraceMapMatrix &traceMap,
												ssize_t edgeIndex) {
	ssize_t t;
	ssize_t traceLength;
	for (t = 0; t < traceMap[edgeIndex].traces.size(); t++ ){
		traceLength = pathTraces[traceMap[edgeIndex].traces[t].trace].size();
		if (traceMap[edgeIndex].traces[t].pos == 0 or
				traceMap[edgeIndex].traces[t].pos == traceLength - 1) 
			return 0;
	}
	return 1;
}
									
ssize_t AreEntranceAndExitEdgesPaired(PathTraceList &pathTraces,
																	TraceMapMatrix &traceMap,
																	TVertexList &vertices,
																	TEdgeList &edges,
																	ssize_t centerEdgeIndex, std::set<InOutEdgePair> &edgePairs) {

	ssize_t src, dest;
	src  = edges[centerEdgeIndex].src;
	dest = edges[centerEdgeIndex].dest;
	// Make sure we have a clean slate.
	edgePairs.clear();
	if (vertices[src].InDegree() != vertices[dest].OutDegree())
		return 0;

	//UNUSED// ssize_t inMapIndex = 0;
	ssize_t centerTrace;
	for (centerTrace = 0; centerTrace < traceMap[centerEdgeIndex].traces.size(); centerTrace++) {
		// Do not process this path if it starts in the center
		if (traceMap[centerEdgeIndex].traces[centerTrace].pos == 0)
			continue;
				
		ssize_t traceIndex  = traceMap[centerEdgeIndex].traces[centerTrace].trace;
		ssize_t tracePos    = traceMap[centerEdgeIndex].traces[centerTrace].pos;
		ssize_t traceLength = pathTraces[traceIndex].size();

		// This trace ends at the center edge, don't look 
		// to see if it continues to an out edge.
		if (tracePos > traceLength - 1) 
			continue;

		// Check to see if the edge before the center on the trace is the in edge
		edgePairs.insert(InOutEdgePair((*pathTraces[traceIndex].edges)[tracePos-1],
																	 (*pathTraces[traceIndex].edges)[tracePos+1]));
	}
		
	// Now check to see if edges are properly paired.  A necessary
	// condition is that there are as many pairs as in degrees (and out 
	// because of the first check)
	if (edgePairs.size() != vertices[src].InDegree())
		return 0;

	// The sufficient condition is that one of every in edge and out edge is present
	
	ssize_t inEdgeIndex, inEdge;
	ssize_t outEdgeIndex, outEdge;
	std::set<InOutEdgePair>::iterator setIt;
	for (inEdgeIndex = vertices[src].FirstIn();
			 inEdgeIndex != vertices[src].EndIn();
			 inEdgeIndex = vertices[src].NextIn(inEdgeIndex)) {
		inEdge = vertices[src].in[inEdgeIndex];
		
		for (setIt = edgePairs.begin(); setIt != edgePairs.end(); ++setIt) {
			if ((*setIt).first == inEdge)
				break;
		}
		if (setIt == edgePairs.end()) {
			// One of the in edges is not represented, bail out.
			edgePairs.clear();
			return 0;
		}
	}

	for (outEdgeIndex = vertices[dest].FirstIn();
			 outEdgeIndex != vertices[dest].EndOut();
			 outEdgeIndex = vertices[dest].NextOut(outEdgeIndex)) {
		outEdge = vertices[dest].out[outEdgeIndex];
		
		for (setIt = edgePairs.begin(); setIt != edgePairs.end(); ++setIt) {
			if ((*setIt).second == outEdge)
				break;
		}
		if (setIt == edgePairs.end()) {
			// One of the in edges is not represented, bail out.
			edgePairs.clear();
			return 0;
		}
	}
				
	return 1;
}


void DetachPath(PathTraceList &pathTraces,
								TraceMapMatrix &traceMap,
								IntervalGraph &graph,
								TVertexList &vertices, TEdgeList &edges,
								PathIntervalList &paths, PathLengthList &pathLengths,
								PathTrace &trace, ssize_t traceIndex, ssize_t removeInconsistentPathIntervals, ssize_t doPrint) {
	ReadMateList emptyList;
	DetachPath(pathTraces, traceMap, 
						 graph,vertices, edges,
						 paths, pathLengths,
						 trace, traceIndex,
						 0, emptyList, removeInconsistentPathIntervals, doPrint);
}

void DetachPath(PathTraceList &pathTraces,
								TraceMapMatrix &traceMap,
								IntervalGraph &graph,
								TVertexList &vertices, TEdgeList &edges,
								PathIntervalList &paths, PathLengthList &pathLengths,
								PathTrace &trace, ssize_t traceIndex, 
								ReadMateList &mateList, ssize_t removeInconsistentPathIntervals, ssize_t doPrint) {
	DetachPath(pathTraces, traceMap, 
						 graph,vertices, edges,
						 paths, pathLengths,
						 trace, traceIndex,
						 1, mateList, removeInconsistentPathIntervals, doPrint);
}




void TrimPathEnds(IntervalGraph &graph, ssize_t path, ssize_t revPathBegin, ssize_t forPathEnd) {
	ssize_t pi;
	for (pi = 0; pi < revPathBegin; pi++ ){
		if (graph.paths[path][pi].edge != -1 and 
				graph.paths[path][pi].index != -1) {
			graph.edges[graph.paths[path][pi].edge].MarkIntervalForRemoval(graph.paths[path][pi].index);
			graph.paths[path][pi].edge = -1;
			graph.paths[path][pi].index = -1;
		}
	}
	for (pi = forPathEnd; pi < graph.pathLengths[path]; pi++) {
		if (graph.paths[path][pi].edge != -1 and
				graph.paths[path][pi].index != -1 ) {
			graph.edges[graph.paths[path][pi].edge].MarkIntervalForRemoval(graph.paths[path][pi].index);
			graph.paths[path][pi].edge = -1;
			graph.paths[path][pi].index = -1;
		}
	}
}


void DetachPath(PathTraceList &pathTraces,
								TraceMapMatrix &traceMap,
								IntervalGraph &graph,
								TVertexList &vertices, TEdgeList &edges,
								PathIntervalList &paths, PathLengthList &pathLengths,
								PathTrace &trace, ssize_t traceIndex, 
								ssize_t isMatePath, ReadMateList &mateList, 
								ssize_t removeInconsistentPathIntervals,  // remove paths from source if they are not con-
								// sistent wthe trace
								ssize_t doPrint) {
	ssize_t srcVertex;
	ssize_t firstDest;
	ssize_t lastSrc;
	ssize_t destVertex;
	ssize_t traceSize = trace.size();

	ssize_t firstEdge = (*trace.edges)[0];
	ssize_t lastEdge  = (*trace.edges)[traceSize - 1];

	ssize_t traceContainsDuplications;
	
	traceContainsDuplications = TraceContainsDuplications(trace);

	assert(trace.size() > 0);
	
	// Although doing this on a trace of size 1 should be a no-op,
	// I'm asserting here since it doesn't make sense to be called
	// on single edges.
	assert(trace.size() > 1);

	// ******************** 
	// Step 1. 
	// Change the connectivity so that the first edge
	// connects to the dest of the last edge.

	srcVertex = edges[(*trace.edges)[0]].src;
	destVertex = edges[(*trace.edges)[traceSize-1]].dest;

	firstDest = edges[(*trace.edges)[0]].dest;
	lastSrc   = edges[(*trace.edges)[traceSize-1]].src;

	// ********************
	// Step 2. 
	// Remove this trace from all edge maps that reference it.
	// Each edge has a list of traces that it is part of.  This trace
	// is resolved, so it should be removed from the lists in these edges.
	//
	ssize_t tracePos;
	ssize_t traceEdge;
	ssize_t traceMapIndex;
	for (tracePos = 0; tracePos < traceSize; tracePos++ ) {
		traceEdge = (*trace.edges)[tracePos];
		traceMapIndex = 0; 
		while ( traceMapIndex < traceMap[traceEdge].traces.size()) { 
			if (traceMap[traceEdge].traces[traceMapIndex].trace == traceIndex) {
				traceMap[traceEdge].traces.erase(traceMap[traceEdge].traces.begin() + 
																				 traceMapIndex);
			}
			else {
				traceMapIndex++;
			}
		}
	}
	
	// ********************
	// Step 3. 
	// The last edge becomes merged with the first edge of the trace.
	// Update traces that use the last edge to reference the first
	// edge instead.
	//
	ssize_t lastEdgeTrace;
	for (lastEdgeTrace = 0; lastEdgeTrace < traceMap[lastEdge].traces.size(); 
			 lastEdgeTrace++ ) {
		ssize_t lastEdgePathTrace = traceMap[lastEdge].traces[lastEdgeTrace].trace;
		ssize_t lastEdgePathPos   = traceMap[lastEdge].traces[lastEdgeTrace].pos;
		if (lastEdgePathPos == 0) {
			(*pathTraces[lastEdgePathTrace].edges)[0] = firstEdge;
		}
		traceMap[lastEdge].startEdge = firstEdge;

		// Add the 
		if (lastEdgePathTrace != traceIndex)
			traceMap[firstEdge].traces.push_back(TraceMap(lastEdgePathTrace, lastEdgePathPos));
	}
	

	// *******************
	// Step 4.
	// 
	// Move intervals from the edges along the trace into the new edge, and
	// do some cleaning of inconsistent intervals if that flag is set.
	//
	//UNUSED// ssize_t intvIndex, intvEdge;

	ssize_t path, pathPos, pathLength;

	// Grow the intervals from the beginning of the trace into the trace edge
	//UNUSED// ssize_t i;
	//UNUSED// ssize_t curPathEnd;
	//UNUSED// ssize_t numGrew = 0;
	//UNUSED// ssize_t numAdded = 0;
	//UNUSED// ssize_t curPathBegin;

	// Add the intervals from the trace to the new edge.
	// Transform intervals that start in the first edge

	ssize_t numMarkedForRemoval;


	//
	// Make sure none of the intervals along the path 
	// have been designated as detached yet.
	//	

	for (tracePos = 0; tracePos < traceSize; tracePos++ ){
		traceEdge = (*trace.edges)[tracePos];
		ReadIntervalList *traceEdgeIntervals;
		ssize_t numTraceEdgeIntervals;
		traceEdgeIntervals = edges[traceEdge].intervals;
		numTraceEdgeIntervals = edges[traceEdge].intervals->size();
		ssize_t i;
		for (i = 0; i < numTraceEdgeIntervals; i++) {
			(*traceEdgeIntervals)[i].detached = 0;
		}
	}


	for (tracePos = 0; tracePos < traceSize; tracePos++ ){ 
		traceEdge = (*trace.edges)[tracePos];
		numMarkedForRemoval = 0;

		//UNUSED// ssize_t numMoved = 0;
		ssize_t traceEdgeIntv;
		ssize_t numTraceEdgeIntervals = edges[traceEdge].intervals->size();
	
		for (traceEdgeIntv = 0; traceEdgeIntv < numTraceEdgeIntervals; traceEdgeIntv++) { 
			// Don't process this interval if it is set to be removed.
			// it's possible that this edge is traversed multiple times 
			// by a trace, in which case they will be processed when starting 
			// on the first interval in this edge.

			path    = (*edges[traceEdge].intervals)[traceEdgeIntv].read;
			pathPos = (*edges[traceEdge].intervals)[traceEdgeIntv].pathPos;

			if ((*edges[traceEdge].intervals)[traceEdgeIntv].IsMarkedForDeletion())
				continue;

			// This should be a valid path.
			assert(path != -1);
			assert(pathPos != -1);
			// remember where this path started so we can update it
			// in the path list
			//UNUSED// ssize_t pathStartPos = pathPos;
			//UNUSED// ssize_t p, pi;
			
			pathLength   = pathLengths[path];

			// First check to see if the path here is consistent with the thread that 
			// starts at this edge. 
			// We should only be separating edges that have consistent paths.
			ssize_t forPathEnd = -1, revPathBegin = -1;
			ssize_t forConsistent, revConsistent;

			//
			//  By default, try and detach any path.  If detaching mate-paths, only
			//  detach paths that are linked to the trace by mate-pairs, or that overlap
			//  the start or the end of the trace.
			//

			ssize_t readPathIsPaired = 0;

			if (isMatePath) {

				/*
				 * If detaching a trace that is generated from mate-ends, there 
				 * may be paths on the trace that do not span the entire trace.  Of 
				 * these paths, some belong to the trace (are part of the repeat copy 
				 * being resolved), and others are part of another repeat copy.  Pull
				 * out the paths that belong to the copy being resolved if they
				 * have a mate-pair linking them outside this repeat.
				 */


				ssize_t mateReadIndex;
				ssize_t readPathPos;
				if (PathContainsEdge(paths[path], pathLength, lastEdge, readPathPos) or
						PathContainsEdge(paths[path], pathLength, firstEdge, readPathPos)) {
					// This guarantees the path should be detached into the
					// new trace, since it overlaps with the resolved edges of the 
					// trace (the first or last).
					readPathIsPaired = 1;
				}
				//
				// Otherwise, it is necessary to look at the mate paths to see if 
				// they land in the first or last edges that are being transformed.
				//
				else if (mateList[path / 2].mateIndex == -1) {
					// No op here. This read cannot be placed on a mate-path 
					// because it does not start in the first edge, does 
					// not end in the last edge, and there is no mate
					// linking it to either edge.
					readPathIsPaired = 0;
				}
				else {
					ssize_t forwardRead;
					if (path % 2 == 0) {
						// This is a forward read, look to the end edge to see if
						// the read-path is mapped to this mate-path
						mateReadIndex = mateList[path / 2].mateIndex * 2 + 1;
						forwardRead = 1;
					}
					else {
						//
						// This is a reverse read, look on the forward stand 
						// for the mate-read.
						mateReadIndex = mateList[path / 2].mateIndex * 2;
						forwardRead = 0;
					}


					//
					//  Check to see if the mate contains one of the resolved
					//  (first or last) edges of the trace.  If so, the path 
					//  should be detached into the trace.
					//
					ssize_t mateFirstEdgePos = -1, mateLastEdgePos  = -1;
						

					// 
					//  If the read is forward, the read must have a mate overlapping the 
					//  last edge (case 1), since the path traversed the mate after traversing
					//  the repeat.  
					//
					//  (case1)                     FORWARDREAD-----------MATEREAD
					//            firstEdge----repeat1--repeat2----------lastEdge
					//  (case2)     MATEREAD---------REVERSEREAD
					//
					//  In (case2), the read is a reverse read, so the mate must be on the path 
					//  before the repeat.
					if ((forwardRead and PathContainsEdge(paths[mateReadIndex], pathLengths[mateReadIndex], 
																								lastEdge, mateLastEdgePos)) or 
							(!forwardRead and PathContainsEdge(paths[mateReadIndex], pathLengths[mateReadIndex], 
																								 firstEdge, mateFirstEdgePos)))  {
						//
						// The read path is mapped to the ends of the trace that are 
						// resolved, so it should be pulled out into the detached trace.
						//
						if (mateFirstEdgePos != -1 or mateLastEdgePos != -1) {
							readPathIsPaired = 1;
						}
					
						//
						// Now do some error checking on the mated path that links 
						// the read path into the trace. 
						//

						// 
						// It is possible the path of the mate-read is not consistent
						// with the trace-path.  If this is the case, remove the mate-path.
						//
						ssize_t mateForConsistent = 1;
						ssize_t mateRevConsistent = 1;
						ssize_t mateForPathEnd, mateRevPathBegin;
						mateForPathEnd = mateRevPathBegin = -1;
						if (mateFirstEdgePos != -1) {
							if (!traceContainsDuplications) {
								mateForConsistent = 
									ArePathAndTraceForwardConsistent(paths[mateReadIndex], 
																									 mateFirstEdgePos, pathLengths[mateReadIndex], 
																									 trace, 0, mateForPathEnd);
							}
							else {
								mateForConsistent =
									FindLongestPathTraceForwardConsistency(paths[mateReadIndex], 
																												 mateFirstEdgePos, pathLengths[mateReadIndex], 
																												 trace, 0, mateForPathEnd);
							}

							// 
							// The path of the mate begins in this trace, but contains
							// sequence that is not consistent with the trace.  Remove that inconsistent
							// sequence.
							if (!mateForConsistent and removeInconsistentPathIntervals) {
								TrimPathEnds(graph, mateReadIndex, 0, mateForPathEnd);
							}

						}
						if (mateLastEdgePos != -1) {
							if (!traceContainsDuplications) {
								mateRevConsistent = 
									ArePathAndTraceReverseConsistent(paths[mateReadIndex], 
																									 mateLastEdgePos, pathLengths[mateReadIndex],
																									 trace, trace.edges->size() - 1, mateRevPathBegin);
							}
							else {
								mateRevConsistent =
									FindLongestPathTraceReverseConsistency(paths[mateReadIndex], 
																												 mateLastEdgePos, pathLengths[mateReadIndex],
																												 trace, trace.edges->size() - 1, mateRevPathBegin);
							}

							if (!mateRevConsistent and removeInconsistentPathIntervals) {
								TrimPathEnds(graph, mateReadIndex, mateRevPathBegin, pathLengths[mateReadIndex]);
							}
						}
					}
				}
			}

			//
			// Look to see if this read path is consistent with the trace.
			// From the previous step, the read path does not overlap
			// with the beginning of the trace, so it must be consistent
			// with the end of the trace.  
			// 

			forConsistent = revConsistent = 0;
			if (!isMatePath or readPathIsPaired) {
				if (!traceContainsDuplications) {
					forConsistent = ArePathAndTraceForwardConsistent(paths[path], 
																													 pathPos, pathLength, 
																													 trace, tracePos, forPathEnd);
					
					revConsistent = ArePathAndTraceReverseConsistent(paths[path], 
																													 pathPos, pathLength, 
																													 trace, tracePos, revPathBegin);
				}
				else {
					ssize_t longestOverlap;
					ssize_t longestPathBegin, longestPathEnd, longestTraceBegin, longestTraceEnd;
					longestOverlap = FindLongestPathTraceOverlap(paths[path], pathLength, 
																											 longestPathBegin, longestPathEnd,
																											 trace, 
																											 longestTraceBegin, longestTraceEnd, (path % 2)==0);
					
					// If there are duplications in the repeat being resolved (tandem repeat)
					// the path may map to the second duplication. 
					ssize_t overlapTracePos = tracePos;
					if (longestOverlap > 0) {
						pathPos = longestPathBegin;
						overlapTracePos = longestTraceBegin;
					}

					forConsistent = FindLongestPathTraceForwardConsistency(paths[path], 
																																 pathPos, pathLength, 
																																 trace, overlapTracePos, forPathEnd);
					
					revConsistent = FindLongestPathTraceReverseConsistency(paths[path], 
																																 pathPos, pathLength, 
																																 trace, overlapTracePos, revPathBegin);
				}
			}
			

			//
			// Some mate-paths are mapped to this path, but stray from the trace due 
			// to sequencing erros.  When this is the case, and the flag is set
			// to remove the intervals that are inconsistent, trim a path until it is
			// consistent.
			if (isMatePath and 
					readPathIsPaired and
					removeInconsistentPathIntervals and 
					(!forConsistent || !revConsistent)) {
				
				TrimPathEnds(graph, path, revPathBegin, forPathEnd);
				// 
				// This path has now been trimmed so that it is both
				// forward consistent and reverse consistent with the trace.
				// Set these flags to 1 so that the path may be detached
				// to the new trace edge.
				//
				forConsistent = revConsistent = 1;
			}

			//
			// The path is either consistent up to the beginning of the 
			// trace, or part of it was removed after pathPos.
			//
			if (revPathBegin < pathPos)
				revPathBegin = pathPos;

			if (tracePos > 0 and   // Only move intervals after the 0'th trace edge, since 
					                   // the intervals are moving to that edge.
					forConsistent and  // only move intervals in consistent paths.
					revConsistent and
					revPathBegin < forPathEnd and // A NULL path (r >= f) is not worth detaching
					revPathBegin == pathPos // Only detach if this path was not trimmed at pathPos
					                        // if it was trimmed, the intervals will be moved later.
					)	{
					
				//
				// Move all intervals to the new edge.
				//
				for (pathPos = revPathBegin; pathPos < forPathEnd; ++pathPos) {

					//
					// Add the interval at path,pathPos to the new edge.
					//
					ssize_t pathEdge, pathIntv;
					pathEdge = paths[path][pathPos].edge;
					pathIntv = paths[path][pathPos].index;
					
					// This interval is already on the first edge, no need to put it there 
					// again.
					if (pathEdge == firstEdge)
						continue;
					

					//
					// Mark this interval as being detached.  This will be moved
					// to the first edge in a later phase.
					//
					(*edges[pathEdge].intervals)[pathIntv].detached = 1;
				}
			} // End checking to see if this path should be removed.
		} // Done looking through all intervals on the current edge

	} // Done looking through all edges on the trace.


	ssize_t firstEdgeLength = edges[firstEdge].length;
	ssize_t traceEdgeOffset = firstEdgeLength - vertices[firstDest].vertexSize;

	// 
	// Move all detached intervals to the first edge.
	//
	for (tracePos = 1; tracePos < traceSize; tracePos++ ){
		traceEdge = (*trace.edges)[tracePos];
		ReadIntervalList *traceEdgeIntervals;
		ssize_t numTraceEdgeIntervals;
		traceEdgeIntervals    = edges[traceEdge].intervals;
		numTraceEdgeIntervals = edges[traceEdge].intervals->size();
		ssize_t i;
		for (i = 0; i < numTraceEdgeIntervals; i++) {

			if ((*traceEdgeIntervals)[i].detached == 1 and 
					(*traceEdgeIntervals)[i].markedForDeletion == 0) {
				//UNUSED// ssize_t pathEdge, pathIndex;
				ssize_t path, pathPos;
				path = (*traceEdgeIntervals)[i].read;
				pathPos = (*traceEdgeIntervals)[i].pathPos;
				
				ssize_t newIntvIndex = edges[firstEdge].intervals->size();
				(*edges[firstEdge].intervals).push_back((*traceEdgeIntervals)[i]);

				//
				// Now fix the fields that should not be the same.
				// The new interval is offset into the new edge according to 
				// the path length before the interval (traceEdgeOffset).
				//
				// Also, the old interval was both marked for deltion and to be
				// detached.  The new interval should not be deleted and should
				// stay put.

				(*edges[firstEdge].intervals)[newIntvIndex].edgePos = traceEdgeOffset + 
					(*traceEdgeIntervals)[i].edgePos;
				(*edges[firstEdge].intervals)[newIntvIndex].markedForDeletion = 0;
				(*edges[firstEdge].intervals)[newIntvIndex].detached = 0;

				// Mark the old interval for removal.
				(*traceEdgeIntervals)[i].markedForDeletion = 1;
				
				// update the path reference to the edge.
				paths[path][pathPos].edge  = firstEdge;
				paths[path][pathPos].index = newIntvIndex;
			}
		}
		traceEdgeOffset += edges[traceEdge].length - vertices[edges[traceEdge].dest].vertexSize;
	}

	// The only edge that is guaranteed to be removed is the 
	// last edge, which will have a different balanced edge
	// after the transformation.
	edges[firstEdge].balancedEdge = edges[lastEdge].balancedEdge;
	if (lastEdge != firstEdge)
		edges[lastEdge].balancedEdge = -1;
	
	// 1.0  Update the sequence to contain the sequences
	//      of all edges.
	SimpleSequence traceSequences;
	CollectPathSequence(vertices, edges, trace, traceSequences);
	delete[] edges[firstEdge].seq.seq;
	edges[firstEdge].seq.seq = traceSequences.seq;
	edges[firstEdge].seq.length = traceSequences.length;
	edges[firstEdge].length = traceSequences.length;
	
	// 1.1 Update the connectivity.
	// link the new dest.
	if (lastEdge != firstEdge) {
		ssize_t destInIndex = vertices[edges[lastEdge].dest].LookupInIndex(lastEdge);
		vertices[edges[lastEdge].dest].in[destInIndex] = firstEdge;

		// unlink the first edge from the first src
		ssize_t firstDestInIndex = vertices[firstDest].LookupInIndex(firstEdge);
		vertices[firstDest].in[firstDestInIndex] = -1;
		
		ssize_t lastSrcOutIndex = vertices[lastSrc].LookupOutIndex(lastEdge);
		vertices[lastSrc].out[lastSrcOutIndex] = -1;

		// Re-point the first edges
		edges[firstEdge].dest = edges[lastEdge].dest;

		// unlink the last edge

		edges[lastEdge].src  = -1;
		edges[lastEdge].dest = -1;
	}
	/*	
		SortReadIntervalsByReadPos(*edges[firstEdge].intervals);
		graph.UpdatePathIndices(firstEdge);
	*/
}


ssize_t AreTracesForwardConsistent(PathTrace &traceA, ssize_t traceAPos, 
															 PathTrace &traceB, ssize_t traceBPos) {

	ssize_t a, b;
	for (a = traceAPos, b = traceBPos; 
			 a < traceA.edges->size() and b < traceB.edges->size();
			 ++a, ++b) {
		if ((*traceA.edges)[a] != (*traceB.edges)[b]) {
			return 0;
		}
	}
	return 1;
}

ssize_t AreTracesReverseConsistent(PathTrace &traceA, ssize_t traceAPos, 
															 PathTrace &traceB, ssize_t traceBPos) {

	ssize_t a, b;
	for (a = traceAPos, b = traceBPos; 
			 a >= 0 and b >= 0;
			 --a, --b) {
		if ((*traceA.edges)[a] != (*traceB.edges)[b]) {
			return 0;
		}
	}
	return 1;
}


ssize_t AreTracesConsistent(PathTrace &traceA, ssize_t traceAPos, 
												PathTrace &traceB, ssize_t traceBPos) {

	return (AreTracesForwardConsistent(traceA, traceAPos, traceB, traceBPos) and
					AreTracesReverseConsistent(traceA, traceAPos, traceB, traceBPos));
}


ssize_t ArePathAndTraceReverseConsistent(PathInterval *path,
																		 ssize_t pathPos, ssize_t pathLength,
																		 PathTrace &trace, ssize_t tracePos,
																		 ssize_t &pathBegin) {
	ssize_t t, p;
	pathBegin = pathPos;
	for (t = tracePos, p = pathPos;
			 t >= 0 and p >= 0 and path[p].edge == (*trace.edges)[t];
			 t--) {

		// REMOVE THIS WHEN PATHS ARE CONDENSED
		// rewind while the path is the same.
		while (p > 0 and path[p].edge != -1 and path[p].edge == path[p-1].edge ) {
			p--;
			pathBegin--;
		}
		// Move p back.
		p--;
		pathBegin--;
	}
	// Always moved pathBegin past the beginning of the path, or 
	// to one before where the path and trace match.
	// Fix that here.
	pathBegin++;
	if (t < 0 or p < 0) {
		return 1;
	}
	else {
		return 0;
	}
}


ssize_t TraceContainsDuplications(PathTrace &trace ){ 
	std::vector<ssize_t> traceEdges;
	traceEdges = *trace.edges;
	std::sort(traceEdges.begin(), traceEdges.end());
	ssize_t i;
	for (i = 0; i < traceEdges.size() - 1; i++) {
		if (traceEdges[i] == traceEdges[i+1])
			return 1;
	}
	return 0;
}
	

ssize_t  FindLongestPathTraceForwardConsistency(PathInterval *path, 
																						ssize_t pathPos, ssize_t pathLength,
																						PathTrace &trace, ssize_t tracePos, 
																						ssize_t &pathEnd) {

	//
	// Step 1. Look to see if the trace has multiple edges that are the same.
	//
	ssize_t curPathEnd;
	ssize_t maxConsistentPathEnd = -1;
	ssize_t maxPathEnd = -1;
	while (tracePos < trace.size()) {
		if (ArePathAndTraceForwardConsistent(path, pathPos, pathLength, trace, tracePos, curPathEnd)) {
			if (maxConsistentPathEnd == -1 or curPathEnd > maxConsistentPathEnd) 
				maxConsistentPathEnd = curPathEnd;
		}
		else {
			// The trace is not consistent, but it fits in a long portion of
			// the subpath.
			if (maxPathEnd == -1 or curPathEnd > maxPathEnd) {
				maxPathEnd = curPathEnd;
			}
		}
		++tracePos;
	}
	if (maxConsistentPathEnd != -1) {
		pathEnd = maxConsistentPathEnd;
		return 1;
	}
	else {
		pathEnd = maxPathEnd;
		return 0;
	}
}

ssize_t FindLongestPathTraceOverlap(PathInterval *path, ssize_t pathLength, ssize_t &pathBegin, ssize_t &pathEnd,
																PathTrace &trace, 
																ssize_t &traceBegin, ssize_t &traceEnd,
																ssize_t findFirst) {
	// There are faster algorithms to do this, but for now do an exhaustive search.

	ssize_t p, t;

	pathBegin = pathEnd = -1;
	traceBegin = traceEnd = -1;
	std::vector<ssize_t>* traceEdges = trace.edges;
	ssize_t longestMatch = 0, longestPathMatch = 0;
	
	ssize_t inAMatch;
	ssize_t traceLength = trace.size();
	for (t = 0; t < trace.size() and longestMatch < pathLength; t++ ) {
		ssize_t pathSearchEnd = pathLength;
		if (pathSearchEnd + t> trace.size())
			pathSearchEnd -= (pathSearchEnd + t - trace.size());
		
		ssize_t pathMatchBegin = -1, pathMatchEnd = -1;
		//UNUSED// ssize_t matchLength;
		inAMatch = 0;
		ssize_t traceMatchBegin = -1, traceMatchEnd = -1;
		ssize_t ti = t;
		for (p = 0; p < pathLength && ti < traceLength; p++, ti++ ){ 
			if (path[p].edge != (*traceEdges)[ti]) {
				if (inAMatch) {
					// 
					// Closing a match here.
					//
					if (pathMatchEnd - pathMatchBegin > longestPathMatch) {
						traceBegin = traceMatchBegin;
						traceEnd   = traceMatchEnd;
						pathBegin  = pathMatchBegin;
						pathEnd    = pathMatchEnd;

						// Use the trace as the distance measure since
						// the path may contain multiple itervals on the same edge.
						longestMatch = traceMatchEnd - traceMatchBegin;
						longestPathMatch = pathMatchEnd - pathMatchBegin;
					}
				}
				inAMatch = 0;
			}
			else {
				if (inAMatch == 1) {
					pathMatchEnd ++;
					traceMatchEnd++;
				}
				else {
					pathMatchBegin = p;
					pathMatchEnd   = p+1;
					traceMatchBegin = ti;
					traceMatchEnd = ti + 1;
					inAMatch = 1;
				}
			}
			while (inAMatch and 
						 (p < pathLength - 1) and 
						 path[p].edge == path[p+1].edge) {
				p++;
				pathMatchEnd++;
			}
		}
		if (inAMatch) {
			if (pathMatchEnd - pathMatchBegin > longestPathMatch) {
				/*			if ((findFirst and 
					 traceMatchEnd - traceMatchBegin > longestMatch and
					 (pathBegin == -1 or pathBegin > pathMatchBegin)) or 
					(!findFirst and 
					 traceMatchEnd - traceMatchBegin >= longestMatch and
					 (pathBegin == -1 or pathBegin < pathMatchBegin))) {
				*/
				longestPathMatch = pathMatchEnd - pathMatchBegin;
				longestMatch = traceMatchEnd - traceMatchBegin;
				traceBegin = traceMatchBegin;
				traceEnd   = traceMatchEnd;
				pathEnd = pathMatchEnd;
				pathBegin = pathMatchBegin;
			}
		}
	}
	return longestMatch;
}



ssize_t FindLongestPathTraceReverseConsistency(PathInterval *path,
																					 ssize_t pathPos, ssize_t pathLength,
																					 PathTrace &trace, ssize_t tracePos,
																					 ssize_t &pathBegin) {
	ssize_t curPathBegin;
	ssize_t minConsistentPathBegin = trace.size(), minPathBegin = trace.size();
	
	while(tracePos >= 0) {
		if (ArePathAndTraceReverseConsistent(path, pathPos, pathLength, trace, tracePos, curPathBegin)) {
			if (minConsistentPathBegin == -1 or curPathBegin < minConsistentPathBegin) {
				minConsistentPathBegin = curPathBegin;
			}
		}
		else {
			if (minPathBegin == -1 or curPathBegin < minPathBegin) {
				minPathBegin = curPathBegin;
			}
		}
		--tracePos;
	}
	if (minConsistentPathBegin < trace.size()) {
		pathBegin = minConsistentPathBegin;
		return 1;
	}
	else {
		pathBegin = minPathBegin;
		return 0;
	}
}

ssize_t  ArePathAndTraceForwardConsistent(PathInterval *path, 
																			ssize_t pathPos, ssize_t pathLength,
																			PathTrace &trace, ssize_t tracePos, 
																			ssize_t &pathEnd) {
	ssize_t t, p;
	pathEnd = pathPos;
	for (t = tracePos, p = pathPos;
			 t < trace.size() and p < pathLength and path[p].edge == (*trace.edges)[t];
			 t++, p++, pathEnd++ ) {

		// REMOVE THIS WHEN PATHS ARE CONDENSED
		// Advance intervals in the path while they are on the same edge.
		while (p < pathLength - 1 and path[p].edge != -1 and path[p].edge == path[p+1].edge) {
			p++;
			pathEnd++;
		}
	}
	
	if (t == trace.size() or p == pathLength)
		return 1;
	else 
		return 0;
}

void CollectPathSequence(TVertexList &vertices, 
												 TEdgeList &edges,
												 PathTrace &path, 
												 SimpleSequence &seq) {
	ssize_t e;
	ssize_t pathSize = path.size();
	ssize_t edgeIndex;
	ssize_t seqLength = 0;
	for (e = 0; e < pathSize - 1; e++ ) {
		edgeIndex = (*path.edges)[e];
		seqLength += edges[edgeIndex].length - vertices[edges[edgeIndex].dest].vertexSize;
	}
	edgeIndex = (*path.edges)[e];
	seqLength += edges[edgeIndex].length;
	seq.seq = new unsigned char[seqLength];
	seq.length = seqLength;
	ssize_t curPos = 0;
	for (e = 0; e < pathSize - 1; e++ ) {
		edgeIndex = (*path.edges)[e];
		ssize_t curEdgeLength = edges[edgeIndex].length - vertices[edges[edgeIndex].dest].vertexSize;
		memcpy(&seq.seq[curPos], edges[edgeIndex].seq.seq, curEdgeLength);
		curPos += curEdgeLength;
	}

	edgeIndex = (*path.edges)[e];
	memcpy(&seq.seq[curPos], edges[edgeIndex].seq.seq, edges[edgeIndex].length);	
}


ssize_t IsPathEndResolved(PathTrace &trace, TraceMapMatrix &traceMap) {

	ssize_t firstEdge = (*trace.edges)[0];
	ssize_t lastEdge  = (*trace.edges)[trace.size()-1];

	return (traceMap[firstEdge].numStart <= 1 and
					traceMap[firstEdge].numInternal == 0 and
					traceMap[lastEdge].numEnd <= 1 and
					traceMap[lastEdge].numInternal == 0);
}

ssize_t AreInternalPathEdgesResolved(PathTraceList &pathTraces,
																 TraceMapMatrix &traceMap,
																 ssize_t pathIndex) {
	ssize_t pe;
	for (pe = 1; pe < pathTraces[pathIndex].size() - 1; pe++ ) {
		if (!IsTangleEdgeResolved(pathTraces, traceMap, (*pathTraces[pathIndex].edges)[pe]))
			return 0;
	}
	return 1;
} 


ssize_t IsTangleEdgeResolved(PathTraceList &pathTraces,
												 TraceMapMatrix &traceMap,
												 ssize_t edgeIndex) {
	ssize_t edgeTrace;
	for (edgeTrace = 0; edgeTrace < traceMap[edgeIndex].traces.size(); edgeTrace++ ) {
		ssize_t pathIndex = traceMap[edgeIndex].traces[edgeTrace].trace;
		if (!IsPathEndResolved(pathTraces[pathIndex], traceMap))
			return 0;
	}
	return 1;
}



ssize_t FindCompatibleTrace(PathTraceList &pathTraces,
												TraceMapMatrix &traceMap,
												PathTrace &trace, ssize_t traceIndex, ssize_t tracePos,
												ssize_t edgeIndex,
												ssize_t &compatibleTrace, ssize_t &compatibleTracePos) {
	ssize_t t;
	ssize_t numCompatible = 0;
	for (t = 0; t < traceMap[edgeIndex].traces.size(); t++ ) {
		if (traceMap[edgeIndex].traces[t].trace != traceIndex and
				AreTracesConsistent(trace, tracePos,
														pathTraces[traceMap[edgeIndex].traces[t].trace],
														traceMap[edgeIndex].traces[t].pos)) {
			compatibleTrace = traceMap[edgeIndex].traces[t].trace;
			compatibleTracePos = traceMap[edgeIndex].traces[t].pos;
			++numCompatible;
		}
	}
	return numCompatible;
}

ssize_t ExtendPath(PathTraceList &pathTraces,
							 TraceMapMatrix &traceMaps,
							 PathTrace &trace,
							 ssize_t traceIndex, ssize_t tracePos,
							 PathTrace &newTrace) {
	/* 
	 * This is experimental code designed to merge overlapping paths.
	 *
	 */
	
	ssize_t startEdge, endEdge;
	ssize_t traceLength = trace.edges->size();
	startEdge = (*trace.edges)[0];
	endEdge   = (*trace.edges)[traceLength - 1];

	// We should only be considering paths that make sense to us.
	//	assert(traceMaps[startEdge].resolved || traceMaps[endEdge].resolved);
	
	ssize_t pi;
	/*
		std::cout << "extending path: " << traceIndex << " ";
		for (pi = tracePos + 1; pi < traceLength; pi++ ){ 
		std::cout << (*trace.edges)[pi] << " (" << traceMaps[(*trace.edges)[pi]].numStart << " " 
		<< traceMaps[(*trace.edges)[pi]].numInternal << " "
		<< traceMaps[(*trace.edges)[pi]].numEnd << "), ";
		}
		std::cout << std::endl;
	*/
	ssize_t compatibleTrace = -1;
	ssize_t compatibleTracePos = -1;
	ssize_t compatiblePathFound = 0;
	for (pi = tracePos + 1; pi < traceLength and !compatiblePathFound; pi++ ) {
		ssize_t traceEdge = (*trace.edges)[pi];
		/*		if (traceMaps[traceEdge].numStart <= 1 and
			traceMaps[traceEdge].numEnd   <= 1) {
		*/
		if (traceMaps[traceEdge].numEnd   <= 1 and 
				traceMaps[traceEdge].inResolved == 0 and 
				traceMaps[traceEdge].outResolved == 0 ) {
			// Found an edge that does not have too many paths that start/end in it.
			// That means it may be compatible with another path.
			if (FindCompatibleTrace(pathTraces, traceMaps, trace, traceIndex, pi,
															traceEdge, compatibleTrace, compatibleTracePos) == 1) {
				std::cout << "path: " << compatibleTrace << " " << compatibleTracePos << " is consistent."
									<< std::endl;
				ssize_t cp;
				for (cp = 0; cp < pathTraces[compatibleTrace].edges->size(); cp++ ){
					std::cout << (*pathTraces[compatibleTrace].edges)[cp] << " ";
				}
				std::cout << std::endl;
				// only one trace is compatible with this one.
				ssize_t prevLength = trace.edges->size();
				ssize_t compatibleTraceLength = pathTraces[compatibleTrace].edges->size();
				
				// Combine the two traces into a new one.
				ssize_t newLength = pi + (compatibleTraceLength - compatibleTracePos);
				
				if (newLength >= prevLength) {
					// Stop searching for compatible paths after this... there
					// should only be one compatible path though.
					compatiblePathFound = 1;

					// For now I'm not filtering out matching a path with itself (that just
					// needs a new parameter to this function, so nothing big).

					newTrace.edges = new std::vector<ssize_t>;
					ssize_t ni;
					// Copy over the old path up until and including the resolved edge.
					for (ni = 0; ni <= pi; ni++) {
						newTrace.edges->push_back((*trace.edges)[ni]);
					}
					
					// Copy the new path.
					for (ni = compatibleTracePos+1; ni < compatibleTraceLength; ni++ ) {
						newTrace.edges->push_back((*pathTraces[compatibleTrace].edges)[ni]);
					}

					/*					// Now add this path trace to all paths.
						pathTraces.push_back(newTrace);
					
						ssize_t newTraceIndex = pathTraces.size() - 1;
						// Now add this path trace to all maps that index it/
						for (ni = 0; ni < newLength; ni++) { 
						traceMaps[(*newTrace.edges)[ni]].traces.push_back(TraceMap(newTraceIndex, ni));
						}
					*/
				}
			}
		}
	}
	return compatiblePathFound;
}


ssize_t MarkResolvedPaths(PathTraceList &pathTraces, TraceMapMatrix &traceMaps, ssize_t notStrict ){
	ssize_t p;
	ssize_t numResolved = 0;
	for (p = 0; p < pathTraces.size(); p++ ) {
		ssize_t traceEnd  = pathTraces[p].edges->size() - 1;
		ssize_t firstEdge = (*pathTraces[p].edges)[0];
		ssize_t lastEdge  = (*pathTraces[p].edges)[traceEnd];
		if (pathTraces[p].size() < 2)
			continue;

		if (traceMaps[firstEdge].numStart == 1 and 
				traceMaps[firstEdge].numInternal == 0 and 
				traceMaps[lastEdge].numEnd == 1 and
				traceMaps[lastEdge].numInternal == 0) {
			//			std::cout << "possibly resolved" << std::endl;
			if (notStrict or AreInternalPathEdgesResolved(pathTraces, traceMaps, p)) {
				traceMaps[firstEdge].outResolved = 1;
				traceMaps[lastEdge].inResolved   = 1;
				++numResolved;
			}
		}
	}
	return numResolved;
}

ssize_t CheckPathTraceListBalance(TEdgeList &edges, PathTraceList &pathTraces) {
	ssize_t pa1, pa2;
	for (pa1 = 0; pa1 < pathTraces.size(); pa1++) {
		ssize_t balFound = 0;
		for (pa2 = 0; pa2 < pathTraces.size() and !balFound; pa2++) {
			ssize_t pos;
			if (pathTraces[pa1].edges->size() == pathTraces[pa2].edges->size()) {
				ssize_t pathLength = pathTraces[pa1].edges->size();
				for (pos = 0; pos < pathLength; pos++ ){
					if ((*pathTraces[pa1].edges)[pos] != 
							edges[(*pathTraces[pa2].edges)[pathLength - pos - 1]].balancedEdge)
						break;
				}
				if (pos == pathLength)
					balFound = 1;
			}
		}
		if (!balFound) {
			ssize_t pos;
			
			for (pos = 0; pos < pathTraces[pa1].edges->size(); pos++) {
				std::cout << (*pathTraces[pa1].edges)[pos] << " (" 
									<< edges[(*pathTraces[pa1].edges)[pos]].balancedEdge << ") ";
			}
			std::cout << std::endl;
			return 0;
		}
	}
	return 1;
}


ssize_t CheckPathBalance(TEdgeList &edges, PathTraceList &pathTraces, TraceMapMatrix &traceMaps) {
	// Check the balance of all traces.
	ssize_t p;
	for (p = 0; p < pathTraces.size(); p++ ){ 
		ssize_t traceEnd  = pathTraces[p].edges->size() - 1;
		//UNUSED//		int firstEdge = (*pathTraces[p].edges)[0];
		ssize_t lastEdge  = (*pathTraces[p].edges)[traceEnd];
		ssize_t traceLength = (*pathTraces[p].edges).size();
		ssize_t balEdge   = edges[lastEdge].balancedEdge;
		
		// There must exist a path that starts in the balance of the last edge and is the exact
		// balance of this path.
		ssize_t balPathFound = 0;
		ssize_t t;
		for (t = 0; t < traceMaps[balEdge].traces.size() and !balPathFound; t++) {
			ssize_t bpi, bp;
			
			bp = traceMaps[balEdge].traces[t].trace;
			if (pathTraces[bp].edges->size() == traceLength) {
				ssize_t pathsAgree = 1;
				for (bpi = 0; bpi < pathTraces[bp].edges->size() and pathsAgree; bpi++ ) {
					if (edges[(*pathTraces[bp].edges)[bpi]].balancedEdge != (*pathTraces[p].edges)[traceLength - bpi - 1])
						pathsAgree = 0;
				}
				if (pathsAgree)
					balPathFound = 1;
			}
		}
		assert(balPathFound);
	}
	return 1;
}

ssize_t PrintCandidatePaths(PathTraceList &pathTraces,
												TraceMapMatrix &traceMaps,
												TEdgeList &edges) {
	ssize_t p, pi;
	ssize_t numCandidates = 0;
	for (p = 0; p < pathTraces.size(); p++ ) {
		ssize_t traceLength = (*pathTraces[p].edges).size();
		if (traceLength <= 0)
			continue;
		ssize_t firstEdge = (*pathTraces[p].edges)[0];
		//UNUSED// int lastEdge  = (*pathTraces[p].edges)[traceLength-1];
		if (traceMaps[firstEdge].numStart == 1 and
				traceMaps[firstEdge].numInternal == 0 and
				traceMaps[firstEdge].numEnd <= 1) {
			ssize_t candidate = 0;
			for (pi = 2; pi < traceLength; pi++ ){ 
				if (traceMaps[(*pathTraces[p].edges)[pi]].numEnd == 1 and
						traceMaps[(*pathTraces[p].edges)[pi]].numStart == 0) {
					candidate = 1;
				}
			}
			if (candidate) { 
				++numCandidates;
				std::cout << "path: " << p << " is a candidate path." << std::endl;
				for (pi = 0; pi < traceLength; pi++ ) { 
					std::cout << (*pathTraces[p].edges)[pi] << " ("
										<< traceMaps[(*pathTraces[p].edges)[pi]].numStart << " "
										<< traceMaps[(*pathTraces[p].edges)[pi]].numInternal << " "
										<< traceMaps[(*pathTraces[p].edges)[pi]].numEnd << " "
										<< edges[(*pathTraces[p].edges)[pi]].length << ") ";
				}
				std::cout << std::endl;
			}
		}
	}
	return numCandidates;
}

void DeletePathTraceList(PathTraceList &list) {
	ssize_t i;
	for (i =0 ;i< list.size(); i++) {
		delete list[i].edges;
	}
	list.clear();
}

ssize_t PathContainsEdge(PathInterval *path, ssize_t pathLength,
										 ssize_t edgeIndex, ssize_t &pathIndex) {
	ssize_t pi;
	for (pi = 0; pi < pathLength; pi++ ){
		if (path[pi].edge == edgeIndex) {
			pathIndex = pi;
			return 1;
		}
	}
	return 0;
}


void PrintTracesAsReads(TVertexList &vertices, TEdgeList &edges, 
												PathTraceList &traces, ssize_t endEdgeLength,
												std::string pathSeqName,
												std::ostream &report) {

	std::ofstream pathSeqOut;
	openck(pathSeqName, pathSeqOut, std::ios::out, report);

	// Print the paths as reads.
	ssize_t t, e;
	for (t = 0; t < traces.size(); t++ ) {
		// print a long version of the suffix of the last edge.
		ssize_t traceLength = traces[t].edges->size();
		
		if (traceLength == 0) 
			continue;

		pathSeqOut << ">" << t << " " << traceLength << std::endl;
		ssize_t firstDest = edges[(*traces[t].edges)[0]].dest;
		ssize_t lastSrc   = edges[(*traces[t].edges)[traceLength-1]].src;
		ssize_t firstDestVertexSize = vertices[firstDest].vertexSize;
		ssize_t lastSrcVertexSize = vertices[lastSrc].vertexSize;
		ssize_t firstEdgeEndLength = 
			std::min(edges[(*traces[t].edges)[0]].length - firstDestVertexSize, 
							 endEdgeLength);
		
		ssize_t lastEdgeEndLength  = 
			std::min(edges[(*traces[t].edges)[traceLength-1]].length - lastSrcVertexSize, 
							 endEdgeLength);

		
		unsigned char *seqPtr;
		ssize_t s;
		ssize_t firstEdgeLength = edges[(*traces[t].edges)[0]].length;
		seqPtr =  &edges[(*traces[t].edges)[0]].seq.seq[firstEdgeLength - firstEdgeEndLength];

		// output the first edge
		ssize_t pathSeqLength = 0;
		for (s = 0; s < firstEdgeEndLength; s++ ) {
			pathSeqOut << seqPtr[s];
			++pathSeqLength;
			if (pathSeqLength and pathSeqLength % 50 == 0)
				pathSeqOut << std::endl;
		}
		
		// output all intermediate edges.
		ssize_t edgeLength;
		for (e = 1; e < traceLength; e++) {
			ssize_t dest = edges[(*traces[t].edges)[e]].dest;
			int vertexSize = vertices[dest].vertexSize;
											 
			if (e < traceLength - 1)
				edgeLength = edges[(*traces[t].edges)[e]].length - vertexSize;
			else 
				edgeLength = lastEdgeEndLength;

			seqPtr = &edges[(*traces[t].edges)[e]].seq.seq[vertexSize];
			for (s = 0; s < edgeLength; s++ ){ 
				pathSeqOut << seqPtr[s];
				pathSeqLength++;
				if (pathSeqLength and pathSeqLength % 50 == 0)
					pathSeqOut << std::endl;
			}
		}

		if (pathSeqLength % 50 != 0)
			pathSeqOut << std::endl;
	}
}
/*
	ssize_t FindBranch(PathTree &tree, ssize_t branch) {
	ssize_t b;
	for (b = 0; b < tree.size(); b++) {
	if (trunk == tree[b].edge) {
	break;
	}
	}
	assert(b < tree.size());
	return b;
	}

	void AddBranch(PathTree &tree, ssize_t trunk, ssize_t branch) {
	ssize_t b = FindBranch(tree, trunk);
	ssize_t treeSize = tree.size();
	tree[trunk].branches.push_back(treeSize - 1);
	tree.push_back(PathBranch(branch));
	}

	void InitTree(PathTree &tree, ssize_t trunk) {
	tree.insert(PathBranch(trunk));
	}
*/


void BalancedDetachPaths(PathTraceList &pathTraces,
												 TraceMapMatrix &traceMap,
												 IntervalGraph &graph, 
												 TVertexList &vertices, TEdgeList &edges,
												 PathIntervalList &paths, PathLengthList &pathLengths,
												 PathTrace &trace, ssize_t traceIndex,
												 PathTrace &balTrace, ssize_t balTraceIndex) {
	ReadMateList mateList;
	BalancedDetachPaths(pathTraces, traceMap, graph, vertices, edges,
											paths, pathLengths,
											trace, traceIndex, balTrace, balTraceIndex, mateList, 0,0);
}


void BalancedDetachPaths(PathTraceList &pathTraces,
												 TraceMapMatrix &traceMap,
												 IntervalGraph &graph, 
												 TVertexList &vertices, TEdgeList &edges,
												 PathIntervalList &paths, PathLengthList &pathLengths,
												 PathTrace &trace, ssize_t traceIndex,
												 PathTrace &balTrace, ssize_t balTraceIndex,
												 ReadMateList &mateList, ssize_t removeInconsistentPathIntervals, ssize_t isMatePath) {

	ssize_t srcVertex;
	ssize_t firstDest;
	ssize_t lastSrc;
	ssize_t destVertex;
	ssize_t traceSize = trace.size();
	ssize_t balTraceSize = balTrace.size();
	assert(traceSize == balTraceSize);

	ssize_t firstEdge = (*trace.edges)[0];
	ssize_t lastEdge  = (*trace.edges)[traceSize - 1];

	ssize_t balFirstEdge = (*balTrace.edges)[0];
	ssize_t balLastEdge = (*balTrace.edges)[balTraceSize - 1];

	ssize_t traceContainsDuplications;
	ssize_t balTraceContainsDuplications;
	

	//UNUSED// ssize_t e;
	
	traceContainsDuplications = TraceContainsDuplications(trace);
	balTraceContainsDuplications = TraceContainsDuplications(balTrace);
	assert(traceContainsDuplications == balTraceContainsDuplications);

	assert(trace.size() > 0);
	
	// Although doing this on a trace of size 1 should be a no-op,
	// I'm asserting here since it doesn't make sense to be called
	// on single edges.
	assert(trace.size() > 1);

	// ******************** 
	// Step 1. 
	// Change the connectivity so that the first edge
	// connects to the dest of the last edge.

	srcVertex  = edges[(*trace.edges)[0]].src;
	destVertex = edges[(*trace.edges)[traceSize-1]].dest;

	firstDest  = edges[(*trace.edges)[0]].dest;
	lastSrc    = edges[(*trace.edges)[traceSize-1]].src;

	ssize_t balFirstDest;
	balFirstDest = edges[(*balTrace.edges)[0]].dest;
	ssize_t balLastSrc;
	balLastSrc = edges[(*balTrace.edges)[traceSize-1]].src;

	map<ssize_t, IntVector> traceIntervalDest;

	// 
	// Initialize storage for where the trace intervals go.
	//
	ssize_t tracePos, balTracePos;
	ssize_t traceEdge, balTraceEdge;
	for (tracePos = 0; tracePos < traceSize; tracePos++ ) {
		traceEdge = (*trace.edges)[tracePos];
		if (traceIntervalDest.find(traceEdge) == traceIntervalDest.end()) {
			traceIntervalDest[traceEdge].resize(edges[traceEdge].intervals->size());
			std::fill(traceIntervalDest[traceEdge].begin(),
								traceIntervalDest[traceEdge].end(), -1);
		}
	}

	// 
	// This will create all new vectors... UNLESS the path has palindromic
	// edges on it.  Then the vector is reused.
	//
	for (balTracePos = 0; balTracePos < traceSize; ++balTracePos) {
		balTraceEdge = (*balTrace.edges)[balTracePos];
		if (traceIntervalDest.find(balTraceEdge) == traceIntervalDest.end()) {
			traceIntervalDest[balTraceEdge].resize(edges[balTraceEdge].intervals->size());
			std::fill(traceIntervalDest[balTraceEdge].begin(),
								traceIntervalDest[balTraceEdge].end(), -1);
		}		
	}


	// ********************
	// Step 2. 
	// Remove this trace from all edge maps that reference it.
	// Each edge has a list of traces that contain the edge.  Since this
	// trace is resolved, it will no longer exist, therefore should be
	// removed from the map of edges->trace.
	//


	ssize_t traceMapIndex, balTraceMapIndex;
	for (tracePos = 0; tracePos < traceSize; tracePos++ ) {
		traceEdge = (*trace.edges)[tracePos];
		traceMapIndex = 0; 
		while ( traceMapIndex < traceMap[traceEdge].traces.size()) { 
			if (traceMap[traceEdge].traces[traceMapIndex].trace == traceIndex) {
				traceMap[traceEdge].traces.erase(traceMap[traceEdge].traces.begin() + 
																				 traceMapIndex);
			}
			else {
				traceMapIndex++;
			}
		}
		balTraceEdge = (*balTrace.edges)[tracePos];
		balTraceMapIndex = 0;
		while (balTraceMapIndex < traceMap[balTraceEdge].traces.size()) {
			if (traceMap[balTraceEdge].traces[balTraceMapIndex].trace == balTraceIndex) {
				traceMap[balTraceEdge].traces.erase(traceMap[balTraceEdge].traces.begin() + 
																						balTraceMapIndex);
			}
			else {
				balTraceMapIndex++;
			}
		}
	}
	
	// ********************
	// Step 3. 
	// The last edge becomes merged with the first edge of the trace.
	// Update traces that use the last edge to reference the first
	// edge instead.
	//
	ssize_t lastEdgeTrace;
	for (lastEdgeTrace = 0; 
			 lastEdgeTrace < traceMap[lastEdge].traces.size(); 
			 lastEdgeTrace++ ) {
		ssize_t lastEdgePathTrace = traceMap[lastEdge].traces[lastEdgeTrace].trace;
		ssize_t lastEdgePathPos   = traceMap[lastEdge].traces[lastEdgeTrace].pos;
		if (lastEdgePathPos == 0) {
			(*pathTraces[lastEdgePathTrace].edges)[0] = firstEdge;
		}
		traceMap[lastEdge].startEdge = firstEdge;

		// Add the reference to the trace that started in the last edge
		// to the first edge.
		if (lastEdgePathTrace != traceIndex)
			traceMap[firstEdge].traces.push_back(TraceMap(lastEdgePathTrace, lastEdgePathPos));
	}

	//
	// Perform the balanced operation to the balanced path.
	// Balancd paths are merged into the last edge instead
	// of the first since the last edge of the balanced 
	// trace is the balance of the first edge of the original trace
	// So, any path trace that reference the first edge should 
	// be changed to reference the last edge.
	ssize_t balFirstEdgeTrace;
	for (balFirstEdgeTrace = 0; 
			 balFirstEdgeTrace < traceMap[balFirstEdge].traces.size(); 
			 balFirstEdgeTrace++ ) {
		ssize_t balFirstEdgePathTrace = traceMap[balFirstEdge].traces[balFirstEdgeTrace].trace;
		ssize_t balFirstEdgePathPos   = traceMap[balFirstEdge].traces[balFirstEdgeTrace].pos;
		if (balFirstEdgePathPos == pathTraces[balFirstEdgePathTrace].edges->size()-1) {
			(*pathTraces[balFirstEdgePathTrace].edges)[balFirstEdgePathPos] = balLastEdge;
		}

		traceMap[balFirstEdge].endEdge = balLastEdge;

		//
		// Add the reference to the trace that ends in the first edge
		// to the map of traces in the last edge since traces that
		// ended in the bal first edge now end in the bal last edge.
		//
		if (balFirstEdgePathTrace != balTraceIndex)
			traceMap[balLastEdge].traces.push_back(TraceMap(balFirstEdgePathTrace, balFirstEdgePathPos));
	}


	

	// *******************
	// Step 4.
	// 
	// Move intervals from the edges along the trace into the new edge, and
	// do some cleaning of inconsistent intervals if that flag is set.
	//
	//UNUSED// ssize_t intvIndex, intvEdge;

	//UNUSED+// ssize_t pathLength;
	ssize_t path, pathPos ;

	// Grow the intervals from the beginning of the trace into the trace edge
	ssize_t i;
	//UNUSED// ssize_t curPathEnd;
	//UNUSED// ssize_t numGrew = 0;
	//UNUSED// ssize_t numAdded = 0;
	//UNUSED// ssize_t curPathBegin;

	// Add the intervals from the trace to the new edge.
	// Transform intervals that start in the first edge

	ssize_t numMarkedForRemoval;

	//
	// NOW, only process the path information for the forward trace,
	// and 

	for (tracePos = 0; tracePos < traceSize; tracePos++ ){ 
		traceEdge = (*trace.edges)[tracePos];
		numMarkedForRemoval = 0;

		//UNUSED// ssize_t numMoved = 0;
		ssize_t traceEdgeIntv;
		ssize_t numTraceEdgeIntervals = edges[traceEdge].intervals->size();
		ssize_t balPath, balPathPos;
		ssize_t pathLength;
		ssize_t balEdge, balEdgeIntv;
		for (traceEdgeIntv = 0; traceEdgeIntv < numTraceEdgeIntervals; traceEdgeIntv++) { 
			// Don't process this interval if it is set to be removed.
			// it's possible that this edge is traversed multiple times 
			// by a trace, in which case they will be processed when starting 
			// on the first interval in this edge.

			path    = (*edges[traceEdge].intervals)[traceEdgeIntv].read;
			pathPos = (*edges[traceEdge].intervals)[traceEdgeIntv].pathPos;
			pathLength = pathLengths[path];

			if (pathLength == 0)
				continue;
			//
			// Determine how to index into the balanced path.
			//
			if (path % 2 == 0) {
				balPath = path + 1;
			}
			else {
				balPath = path - 1;
			}
			balPathPos = pathLength - pathPos - 1;
			balEdge = paths[balPath][balPathPos].edge;
			balEdgeIntv = paths[balPath][balPathPos].index;

			if ((*edges[traceEdge].intervals)[traceEdgeIntv].IsMarkedForDeletion()) {
				continue;
			}

			// This should be a valid path.
			assert(path != -1);
			assert(pathPos != -1);
			// remember where this path started so we can update it
			// in the path list
			//UNUSED// ssize_t pathStartPos = pathPos;
			//UNUSED+// ssize_t  p;
			ssize_t pi;
			
			// First check to see if the path here is consistent with the thread that 
			// starts at this edge. 
			// We should only be separating edges that have consistent paths.
			ssize_t forPathEnd = -1, revPathBegin = -1;
			ssize_t forConsistent, revConsistent;

			//
			//  By default, try and detach any path.  If detaching mate-paths, only
			//  detach paths that are linked to the trace by mate-pairs, or that overlap
			//  the start or the end of the trace.
			//

			ssize_t readPathIsPaired = 0;

			if (isMatePath) {

				/*
				 * If detaching a trace that is generated from mate-ends, there 
				 * may be paths on the trace that do not span the entire trace.  Of 
				 * these paths, some belong to the trace (are part of the repeat copy 
				 * being resolved), and others are part of another repeat copy.  Pull
				 * out the paths that belong to the copy being resolved if they
				 * have a mate-pair linking them outside this repeat.
				 */


				ssize_t mateReadIndex;
				ssize_t readPathPos;
				if (PathContainsEdge(paths[path], pathLength, lastEdge, readPathPos) or
						PathContainsEdge(paths[path], pathLength, firstEdge, readPathPos)) {
					// This guarantees the path should be detached into the
					// new trace, since it overlaps with the resolved edges of the 
					// trace (the first or last).
					readPathIsPaired = 1;
				}
				//
				// Otherwise, it is necessary to look at the mate paths to see if 
				// they land in the first or last edges that are being transformed.
				//
				else if (mateList[path / 2].mateIndex == -1) {
					// This condition is hit if there is no mate pair 
					// stored for this read (mateIndex = -1).
					// The read cannot be placed on a mate-path 
					// because it does not start in the first edge, does 
					// not end in the last edge, and there is no mate
					// linking it to either edge.
					readPathIsPaired = 0;
				}
				else {
					ssize_t forwardRead;
					if (path % 2 == 0) {
						// This is a forward read, look to the end edge to see if
						// the read-path is mapped to this mate-path
						mateReadIndex = mateList[path / 2].mateIndex * 2 + 1;
						forwardRead = 1;
					}
					else {
						//
						// This is a reverse read, look on the forward stand 
						// for the mate-read.
						//
						mateReadIndex = mateList[path / 2].mateIndex * 2;
						forwardRead = 0;
					}
					ssize_t balMateReadIndex;
					if (mateReadIndex % 2 == 0)
						balMateReadIndex = mateReadIndex + 1;
					else
						balMateReadIndex = mateReadIndex - 1;

					ssize_t matePathLength = pathLengths[mateReadIndex];

					//
					//  Check to see if the mate contains one of the resolved
					//  (first or last) edges of the trace.  If so, the path 
					//  should be detached into the trace.
					//
					ssize_t mateFirstEdgePos = -1, mateLastEdgePos  = -1;
						

					// 
					//  If the read is forward, the read must have a mate overlapping the 
					//  last edge (case 1), since the path traversed the mate after traversing
					//  the repeat.  
					//
					//  (case1)                     FORWARDREAD-----------MATEREAD
					//            firstEdge----repeat1--repeat2----------lastEdge
					//  (case2)     MATEREAD---------REVERSEREAD
					//
					//  In (case2), the read is a reverse read, so the mate must be on the path 
					//  before the repeat.
					if ((forwardRead and PathContainsEdge(paths[mateReadIndex], pathLengths[mateReadIndex], 
																								lastEdge, mateLastEdgePos)) or 
							(!forwardRead and PathContainsEdge(paths[mateReadIndex], pathLengths[mateReadIndex], 
																								 firstEdge, mateFirstEdgePos)))  {
						//
						// The mate of this read is mapped to one of the ends of
						// the trace that is resolved, so it should be pulled out
						// into the detached trace.
						//
						if (mateFirstEdgePos != -1 or mateLastEdgePos != -1) {
							readPathIsPaired = 1;
						}
					
						//
						// Now do some error checking on the mated path that links 
						// the read path into the trace. 
						//

						// 
						// It is possible the path of the mate-read is not consistent
						// with the trace-path.  If this is the case, remove the mate-path.
						//
						ssize_t mateForConsistent = 1;
						ssize_t mateRevConsistent = 1;
						ssize_t mateForPathEnd, mateRevPathBegin;
						mateForPathEnd = mateRevPathBegin = -1;
						if (mateFirstEdgePos != -1) {
							if (!traceContainsDuplications) {
								mateForConsistent = 
									ArePathAndTraceForwardConsistent(paths[mateReadIndex], 
																									 mateFirstEdgePos, pathLengths[mateReadIndex], 
																									 trace, 0, mateForPathEnd);
							}
							else {
								mateForConsistent =
									FindLongestPathTraceForwardConsistency(paths[mateReadIndex], 
																												 mateFirstEdgePos, pathLengths[mateReadIndex], 
																												 trace, 0, mateForPathEnd);
							}

							// 
							// The path of the mate begins in this trace, but contains
							// sequence that is not consistent with the trace.  Remove that inconsistent
							// sequence.
							if (!mateForConsistent and removeInconsistentPathIntervals) {
								TrimPathEnds(graph, mateReadIndex, 0, mateForPathEnd);

								ssize_t balMateRevPathBegin = matePathLength - mateForPathEnd;
								TrimPathEnds(graph, balMateReadIndex, balMateRevPathBegin, matePathLength);
							}

						}
						if (mateLastEdgePos != -1) {
							if (!traceContainsDuplications) {
								mateRevConsistent = 
									ArePathAndTraceReverseConsistent(paths[mateReadIndex], 
																									 mateLastEdgePos, pathLengths[mateReadIndex],
																									 trace, trace.edges->size() - 1, mateRevPathBegin);
							}
							else {
								mateRevConsistent =
									FindLongestPathTraceReverseConsistency(paths[mateReadIndex], 
																												 mateLastEdgePos, pathLengths[mateReadIndex],
																												 trace, trace.edges->size() - 1, mateRevPathBegin);
							}

							if (!mateRevConsistent and removeInconsistentPathIntervals) {
								TrimPathEnds(graph, mateReadIndex, mateRevPathBegin, pathLengths[mateReadIndex]);

								//
								// Do the reverse to the balanced path.
								//
								ssize_t balMateForPathEnd = matePathLength - mateRevPathBegin;
								TrimPathEnds(graph, balMateReadIndex, 0, balMateForPathEnd);
							}
						}
					}
				}
			}

			//
			// Look to see if this read path is consistent with the trace.
			// From the previous step, the read path does not overlap
			// with the beginning of the trace, so it must be consistent
			// with the end of the trace.  
			// 

			forConsistent = revConsistent = 0;
			if (!isMatePath or readPathIsPaired) {
				if (!traceContainsDuplications) {
					forConsistent = ArePathAndTraceForwardConsistent(paths[path], 
																													 pathPos, pathLength, 
																													 trace, tracePos, forPathEnd);
					
					revConsistent = ArePathAndTraceReverseConsistent(paths[path], 
																													 pathPos, pathLength, 
																													 trace, tracePos, revPathBegin);
				}
				else {
					ssize_t overlapTracePos = tracePos;
					if (isMatePath) {
						ssize_t longestOverlap;
						ssize_t longestPathBegin, longestPathEnd, longestTraceBegin, longestTraceEnd;
						longestOverlap = FindLongestPathTraceOverlap(paths[path], pathLength, 
																												 longestPathBegin, longestPathEnd,
																												 trace, 
																												 longestTraceBegin, longestTraceEnd, (path % 2)==0);
						
						// If there are duplications in the repeat being resolved (tandem repeat)
						// the path may map to the second duplication. 

						if (longestOverlap > 0) {
							pathPos = longestPathBegin;
							overlapTracePos = longestTraceBegin;
						}
					}

					forConsistent = FindLongestPathTraceForwardConsistency(paths[path], 
																																 pathPos, pathLength, 
																																 trace, overlapTracePos, forPathEnd);
					
					revConsistent = FindLongestPathTraceReverseConsistency(paths[path], 
																																 pathPos, pathLength, 
																																 trace, overlapTracePos, revPathBegin);
				}
			}
			

			//
			// Some mate-paths are mapped to this path, but stray from the trace due 
			// to sequencing erros.  When this is the case, and the flag is set
			// to remove the intervals that are inconsistent, trim a path until it is
			// consistent.
			if (isMatePath and 
					readPathIsPaired and
					removeInconsistentPathIntervals and 
					(!forConsistent || !revConsistent)) {
				
				TrimPathEnds(graph, path, revPathBegin, forPathEnd);

				ssize_t balRevPathBegin, balForPathEnd;
				balRevPathBegin = pathLength - forPathEnd;
				balForPathEnd   = pathLength - revPathBegin;
				TrimPathEnds(graph, balPath, balRevPathBegin, balForPathEnd);
				// 
				// This path has now been trimmed so that it is both
				// forward consistent and reverse consistent with the trace.
				// Set these flags to 1 so that the path may be detached
				// to the new trace edge.
				//
				forConsistent = revConsistent = 1;
			}

			//
			// The path is either consistent up to the beginning of the 
			// trace, or part of it was removed after pathPos.
			//
			if (revPathBegin < pathPos)
				revPathBegin = pathPos;

			ssize_t pathOverlapsFirst = 0, pathOverlapsLast = 0;


			// 
			// Look to see if the path overlaps the first edge on the trace 
			// or the last edge of the trace.  This is important to decide
			// if the read is part of this path, or some other path that
			// contains the repeat that is being resolved here. 


			for (pi = revPathBegin; pi < forPathEnd; pi++) {
				if (paths[path][pi].edge == firstEdge) {
					pathOverlapsFirst = 1;
					break;
				}
				if (paths[path][pi].edge == lastEdge) {
					pathOverlapsLast = 1;
					break;
				}
			}


			/*			if (tracePos > 0 and   // Only move intervals after the 0'th trace edge, since 
					                   // the intervals are moving to that edge.
														 */
			if (forConsistent and  // only move intervals in consistent paths.
					revConsistent and
					revPathBegin < forPathEnd and // A NULL path (r >= f) is not worth detaching
					revPathBegin == pathPos and // Only detach if this path was not trimmed at pathPos
					                        // if it was trimmed, the intervals will be moved later.
					(isMatePath or                 // IF the path is not mate-paired, it must
					 (!isMatePath and              // overlap with either the front or the end
						(pathOverlapsFirst != 0 or // of the trace.
						 pathOverlapsLast != 0)))
					)	{
					

				//
				// Move all intervals to the new edge.
				//
				for (pathPos = revPathBegin; pathPos < forPathEnd; ++pathPos) {

					//
					// Add the interval at path,pathPos to the new edge.
					//
					ssize_t balPathPos = pathLength - pathPos - 1;
					ssize_t pathEdge, pathIntv, balPathEdge, balPathIntv;
					pathEdge = paths[path][pathPos].edge;
					pathIntv = paths[path][pathPos].index;

					balPathEdge = paths[balPath][balPathPos].edge;
					balPathIntv = paths[balPath][balPathPos].index;
					
					
					// This interval is already on the first edge, no need to put it there 
					// again.
					if (pathEdge == firstEdge)
						continue;
					

					//
					// Mark this interval as being detached.  This will be moved
					// to the first edge in a later phase.
					//
					assert(pathIntv < traceIntervalDest[pathEdge].size());
					assert(balPathIntv < traceIntervalDest[balPathEdge].size());
					traceIntervalDest[pathEdge][pathIntv] = firstEdge;
					traceIntervalDest[balPathEdge][balPathIntv] = balLastEdge;
				}
			} // End checking to see if this path should be removed.
		} // Done looking through all intervals on the current edge
	} // Done looking through all edges on the trace.


	ssize_t firstEdgeLength = edges[firstEdge].length;
	ssize_t traceEdgeOffset = firstEdgeLength - vertices[firstDest].vertexSize;

	// 
	// Move all detached intervals to the first edge.
	//
	for (tracePos = 1; tracePos < traceSize; tracePos++ ){
		traceEdge = (*trace.edges)[tracePos];
		ReadIntervalList *traceEdgeIntervals;
		//UNUSED// ReadIntervalList *balTraceEdgeIntervals;
		ssize_t numTraceEdgeIntervals;
		ssize_t numBalTraceEdgeIntervals;
		traceEdgeIntervals    = edges[traceEdge].intervals;
		numTraceEdgeIntervals = edges[traceEdge].intervals->size();

		balTraceEdge = edges[traceEdge].balancedEdge;
		numBalTraceEdgeIntervals = edges[balTraceEdge].intervals->size();
		//		assert(numTraceEdgeIntervals == numBalTraceEdgeIntervals);
		ssize_t i;
		for (i = 0; i < numTraceEdgeIntervals; i++) {

			if (traceIntervalDest[traceEdge][i] == firstEdge and
					(*traceEdgeIntervals)[i].markedForDeletion == 0) {
				//UNUSED// ssize_t pathEdge, pathIndex;
				ssize_t path, pathPos;
				path = (*traceEdgeIntervals)[i].read;
				pathPos = (*traceEdgeIntervals)[i].pathPos;
				
				ssize_t newIntvIndex = edges[firstEdge].intervals->size();
				(*edges[firstEdge].intervals).push_back((*traceEdgeIntervals)[i]);

				//
				// Now fix the fields that should not be the same.
				// The new interval is offset into the new edge according to 
				// the path length before the interval (traceEdgeOffset).
				//
				// Also, the old interval was both marked for deltion and to be
				// detached.  The new interval should not be deleted and should
				// stay put.

				(*edges[firstEdge].intervals)[newIntvIndex].edgePos = traceEdgeOffset + 
					(*traceEdgeIntervals)[i].edgePos;
				(*edges[firstEdge].intervals)[newIntvIndex].markedForDeletion = 0;
				(*edges[firstEdge].intervals)[newIntvIndex].detached = 0;

				// Mark the old interval for removal.
				(*traceEdgeIntervals)[i].markedForDeletion = 1;
				
				// update the path reference to the edge.
				paths[path][pathPos].edge  = firstEdge;
				paths[path][pathPos].index = newIntvIndex;
			}
		}
		traceEdgeOffset += edges[traceEdge].length - vertices[edges[traceEdge].dest].vertexSize;
	}

	// 
	// In the balanced trace, move all detached intervals to the END of
	// the balanced trace (this is the balanced edge of the first edge
	// of the forward trace). 
	//
	//UNUSED// ssize_t balFirstEdgeLength = edges[balFirstEdge].length;
	traceEdgeOffset = 0;

	ssize_t numOrigLastEdgeIntervals = edges[balLastEdge].intervals->size();

	for (tracePos = 0; tracePos < traceSize-1; tracePos++ ){
		balTraceEdge = (*balTrace.edges)[tracePos];
		ReadIntervalList *traceEdgeIntervals;
		ssize_t numTraceEdgeIntervals;
		traceEdgeIntervals    = edges[balTraceEdge].intervals;
		numTraceEdgeIntervals = edges[balTraceEdge].intervals->size();
		ssize_t i;
		for (i = 0; i < numTraceEdgeIntervals; i++) {
			if (traceIntervalDest[balTraceEdge][i] == balLastEdge and 
					(*traceEdgeIntervals)[i].markedForDeletion == 0) {
				//UNUSED// ssize_t pathEdge, pathIndex;
				ssize_t path, pathPos;
				path = (*traceEdgeIntervals)[i].read;
				pathPos = (*traceEdgeIntervals)[i].pathPos;
				
				ssize_t newIntvIndex = edges[balLastEdge].intervals->size();
				(*edges[balLastEdge].intervals).push_back((*traceEdgeIntervals)[i]);

				//
				// Now fix the fields that should not be the same.
				// The new interval is offset into the new edge according to 
				// the path length before the interval (traceEdgeOffset).
				//
				// Also, the old interval was both marked for deltion and to be
				// detached.  The new interval should not be deleted and should
				// stay put.

				(*edges[balLastEdge].intervals)[newIntvIndex].edgePos = traceEdgeOffset + 
					(*traceEdgeIntervals)[i].edgePos;
				(*edges[balLastEdge].intervals)[newIntvIndex].markedForDeletion = 0;
				(*edges[balLastEdge].intervals)[newIntvIndex].detached = 0;

				// Mark the old interval for removal.
				(*traceEdgeIntervals)[i].markedForDeletion = 1;
				
				// update the path reference to the edge.
				paths[path][pathPos].edge  = balLastEdge;
				paths[path][pathPos].index = newIntvIndex;
			}
		}
		traceEdgeOffset += edges[balTraceEdge].length - vertices[edges[balTraceEdge].dest].vertexSize;
	}

	// The last edge has been grown, so modify the edge positions to 
	// reflect this.
	for (i = 0; i < numOrigLastEdgeIntervals; i++ ){ 
		(*edges[balLastEdge].intervals)[i].edgePos += traceEdgeOffset;
	}

	// 1.0  Update the sequence to contain the sequences
	//      of all edges.
	SimpleSequence traceSequences;
	CollectPathSequence(vertices, edges, trace, traceSequences);
	//UNUSED// ssize_t firstEdgeLenght = edges[firstEdge].length;
	delete[] edges[firstEdge].seq.seq;
	edges[firstEdge].seq.seq = traceSequences.seq;
	edges[firstEdge].seq.length = traceSequences.length;
	edges[firstEdge].length = traceSequences.length;
	
	ssize_t removedLastEdge = 0;
	// 1.1 Update the connectivity.
	// link the new dest.
	if (lastEdge != firstEdge) {

		vertices[edges[lastEdge].dest].AddInEdge(firstEdge);
		// Re-point the first edges
		ssize_t firstEdgeOrigDest = edges[firstEdge].dest;
		edges[firstEdge].dest = edges[lastEdge].dest;

		// unlink the last edge
		cout << "first edge len: " << firstEdgeLength << " last: " << edges[lastEdge].length 
				 << " new: " << edges[firstEdge].length << " before disconnecting: "
				 << vertices[edges[lastEdge].src].InDegree() << " " 
				 << vertices[edges[lastEdge].src].OutDegree() << " "
				 << (ssize_t) (vertices[edges[lastEdge].src].InDegree() <= 
						 vertices[edges[lastEdge].src].OutDegree()) << endl;

		
		if (vertices[edges[lastEdge].src].InDegree() <= 
				vertices[edges[lastEdge].src].OutDegree() ){ 
			// unlink the last edge from the last dest
			ssize_t lastSrcOutIndex = vertices[lastSrc].LookupOutIndex(lastEdge);

			
			cout << "removing last src: " << lastSrc 
					 << " [" << vertices[lastSrc].index << "] degree:" 
					 << vertices[lastSrc].OutDegree() << endl;
			//UNUSED// int deg = vertices[lastSrc].OutDegree();
			vertices[lastSrc].out[lastSrcOutIndex] = -1;


			ssize_t destInIndex = vertices[edges[lastEdge].dest].LookupInIndex(lastEdge);
			vertices[edges[lastEdge].dest].in[destInIndex] = -1;
			edges[lastEdge].src  = -1;
			edges[lastEdge].dest = -1;
			removedLastEdge = 1;
		}

		// unlink the first edge from the first src
		ssize_t firstDestInIndex = vertices[firstEdgeOrigDest].LookupInIndex(firstEdge);
		vertices[firstEdgeOrigDest].in[firstDestInIndex] = -1;

	}



	// BALANCE
	// 1.0  Update the sequence to contain the sequences
	//      of all edges.
	
	SimpleSequence balTraceSequences;
	CollectPathSequence(vertices, edges, balTrace, balTraceSequences);
	delete[] edges[balLastEdge].seq.seq;
	edges[balLastEdge].seq.seq    = balTraceSequences.seq;
	edges[balLastEdge].seq.length = balTraceSequences.length;
	edges[balLastEdge].length     = balTraceSequences.length;
	
	// 1.1 Update the connectivity.
	// link the new src.
	if (balLastEdge != balFirstEdge) {

		//
		// Make the first vertex on the trace reference the bal last edge
		// instead of the bal first edge.
		// 
		vertices[edges[balFirstEdge].src].AddOutEdge(balLastEdge);
		// Re-point the first edges
		ssize_t balLastEdgeOrigSrc = edges[balLastEdge].src;
		edges[balLastEdge].src = edges[balFirstEdge].src;

		// unlink the last edge
		/*		if (vertices[edges[balFirstEdge].dest].OutDegree() == 0) {
			cout << "unlinking: " << balFirstEdge << endl;
		*/
		if (removedLastEdge) {
			//
			// The first dest is replaced by the last dest, so remove the link
			// here.
			//
			//			cout << "unlinking first: " << balFirstEdge << " from: " << edges[balFirstEdge].dest << " " << vertices[edges[balFirstEdge].dest].index << endl;
			ssize_t firstDestInIndex = vertices[edges[balFirstEdge].dest].LookupInIndex(balFirstEdge);
			vertices[edges[balFirstEdge].dest].in[firstDestInIndex] = -1;

			// unlink the first edge
			ssize_t srcOutIndex  = vertices[edges[balFirstEdge].src].LookupOutIndex(balFirstEdge);
			
			//UNUSED// ssize_t o;
			cout << "removing bal first src out: " << edges[balFirstEdge].src 
					 << " [" << vertices[edges[balFirstEdge].src].index << "] degree: " 
					 << vertices[edges[balFirstEdge].src].OutDegree() << endl;
			//UNUSED// int deg = vertices[edges[balFirstEdge].src].OutDegree();
			vertices[edges[balFirstEdge].src].out[srcOutIndex] = -1;

			edges[balFirstEdge].src  = -1;
			edges[balFirstEdge].dest = -1;
		}
		// The last source is replaced by the first source, so unlink the
		// out edge here.
		ssize_t lastSrcOutIndex  = vertices[balLastEdgeOrigSrc].LookupOutIndex(balLastEdge);
		
		//UNUSED// ssize_t o;
		ssize_t outDegree = vertices[balLastEdgeOrigSrc].OutDegree();
		cout << "removing bal last edge orig src " << balLastEdge 
				 << " [" << vertices[balLastEdgeOrigSrc].index << "] degree: "
 				 << outDegree << endl;
		//UNUSED// int deg = vertices[balLastEdgeOrigSrc].OutDegree();
		vertices[balLastEdgeOrigSrc].out[lastSrcOutIndex] = -1;
	}


	// The only edge that is guaranteed to be removed is the 
	// last edge, which will have a different balanced edge
	// after the transformation.
	//	edges[firstEdge].balancedEdge = edges[lastEdge].balancedEdge;
	if (removedLastEdge) {
		if (lastEdge != firstEdge)
			edges[lastEdge].balancedEdge = -1;
		
		
		//	edges[balLastEdge].balancedEdge = edges[balFirstEdge].balancedEdge;
		if (balFirstEdge != balLastEdge) 
			edges[balFirstEdge].balancedEdge = -1;
	}

}
