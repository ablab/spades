/***************************************************************************
 * Title:          DAGSimplification.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/29/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
ssize_t IntervalGraph::VertexSetContained(std::set<ssize_t> &dag, ssize_t sourceVertex) {
  // Check the dag for short intruding edges (edges that 
  // 
  ssize_t v;
  std::set<ssize_t>::iterator dagIt, dagEnd;
  for (dagIt = dag.begin(); dagIt != dag.end(); dagIt++) {
    if (*dagIt != sourceVertex) { 
      // Make sure all source vertices of this vertex are in the dag
      ssize_t e, ei;
      for (ei = vertices[*dagIt].FirstIn();
					 ei < vertices[*dagIt].EndIn();
					 ei = vertices[*dagIt].NextIn(ei)) {
				e = vertices[*dagIt].in[ei];
				if (dag.find(edges[e].src) == dag.end()) {
					//	  std::cout << "short intruding edge found " << e << " -> " << edges[e].src << std::endl;
					return 0;
				}
      }
    }
  }
  return 1;
}


ssize_t IntervalGraph::FindShortContainedDAG(ssize_t sourceVertex, ssize_t maxEdgeLength, 
																				 std::set<ssize_t> &dagVertices,
																				 std::set<ssize_t> &dagEdges,
																				 std::vector<ssize_t> &optimalPath) {

  // This code isn't used for now
  std::cout << "this code was only working in small conditions, re-incorporate with care!"
						<< std::endl;
  assert(0);

  // Modify Dijkstra's method a little bit.  

  ssize_t sinkEdge;
  // Initialize the edge leaving this cycle to not found
  sinkEdge = -1;  

  // In Dijkstra, the p-queue contains all edges in the graph
  // We don't want that since we want to explore short paths
  // starting at a source and ending at a unique sink

  // Dag contains a quickly searched list of verties
  // that are part of the dag.
  std::map<ssize_t, ssize_t> optimalVertexPath;
  std::map<ssize_t, ssize_t> optimalEdgePath;

	// Initialize the priority queue to hold the source vertex
  SearchablePriorityQueue<ScoredVertex<ssize_t> > queue;
  ScoredVertex<ssize_t> scoredVertex;
  scoredVertex.vertex = sourceVertex;
  scoredVertex.score = 0;
  queue.push(scoredVertex);
  optimalVertexPath[sourceVertex] = -1;
  optimalEdgePath[sourceVertex] = -1;
  dagVertices.insert(sourceVertex);
  while(queue.size() > 0) {
    scoredVertex = queue.pop();
    //    std::cout << "got vertex " << scoredVertex.vertex << std::endl;
    // First difference from normal Dijkstra: if there are any in edges
    // that are not already part of the dag.
    if (scoredVertex.vertex != sourceVertex) {
      ssize_t inEdgeIndex, inEdge;
      for (inEdgeIndex = vertices[scoredVertex.vertex].FirstIn();
					 inEdgeIndex < vertices[scoredVertex.vertex].EndIn();
					 inEdgeIndex = vertices[scoredVertex.vertex].NextIn(inEdgeIndex)) {
				inEdge = vertices[scoredVertex.vertex].in[inEdgeIndex];
				//	if (DAG.find(edges[inEdge].src) == DAG.end()){
				if (edges[inEdge].length > maxEdgeLength) {
					//	  std::cout << "INTRUDING EDGE FOUND " << inEdge << " -> " << edges[inEdge].src << std::endl;
					dagVertices.clear();
					dagEdges.clear();
					return 0;
				}
      }
      // done checking, no edges in
    }

    // Try and relax any of the edges going out of the current
    // vertex.
    
    ssize_t outEdgeIndex, outEdge;
    for (outEdgeIndex = vertices[scoredVertex.vertex].FirstOut();
				 outEdgeIndex < vertices[scoredVertex.vertex].EndOut();
				 outEdgeIndex = vertices[scoredVertex.vertex].NextOut(outEdgeIndex)) {
      outEdge = vertices[scoredVertex.vertex].out[outEdgeIndex];
      edges[outEdge].traversed = GraphEdge::Marked;
      if (dagEdges.find(outEdge) == dagEdges.end()) {
				if (edges[outEdge].length > maxEdgeLength or
						vertices[edges[outEdge].dest].OutDegree() == 0) {
					if (sinkEdge == -1) 
						sinkEdge = outEdge;
					else if (sinkEdge != outEdge) {
						//	    std::cout << "Found alternative out edge " << outEdge << " " << sinkEdge << std::endl;
						dagVertices.clear();
						dagEdges.clear();	
						return 0;
					}
				}
				else {
					// The out edge isn't too large.  Add the dest vertex to the dag, or relax the score
					ScoredVertex<ssize_t> dest;
					dest.vertex = edges[outEdge].dest;
					dagEdges.insert(outEdge);


					if (queue.find(dest)) {
						// If this is the optimal path, record that
						if (queue[dest].score > scoredVertex.score + edges[outEdge].multiplicity) {
							queue[dest].score = scoredVertex.score + edges[outEdge].multiplicity;
							optimalVertexPath[dest.vertex] = scoredVertex.vertex;
							optimalEdgePath[dest.vertex] = outEdge;
						}
						// If this is tied for the optimal path, one of the two might
						// have preference
						else if (queue[dest].score == scoredVertex.score + edges[outEdge].multiplicity and
										 edges[outEdge].balPreferred == GraphEdge::Marked) {
							std::cout << "finding dag giving preference to "  << outEdge << std::endl;
							optimalVertexPath[dest.vertex] = scoredVertex.vertex;
							optimalEdgePath[dest.vertex] = outEdge;
						}
						/*
							else: the vertex is already in the queue, but is reachable from another
							path.
						*/
					}
					else {
						if (dagVertices.find(dest.vertex) == dagVertices.end()) {
							// This vertex has not been traversed yet.
							// First add it to the dag
							dagVertices.insert(dest.vertex);

							// Then add it to the priority queue.  It may be relaxed
							// in the future if other paths can reach it
							dest.score = scoredVertex.score + edges[outEdge].multiplicity;
							queue.push(dest);
							optimalVertexPath[dest.vertex] = scoredVertex.vertex;
							optimalEdgePath[dest.vertex] = outEdge;
						}
					}
				}
      }
    }
  } 
  
  // Make sure this dag is not intruded by other edges
  if (! VertexSetContained(dagVertices, sourceVertex)) {
    dagVertices.clear();
    dagEdges.clear();
  }

  // Check the dag for cyclicness
  std::list<ssize_t> path;
  std::cout << "checking set of size " << dagVertices.size() << " to see if it is cyclic " << std::endl;

  std::set<ssize_t>::iterator dit;
  for (dit = dagVertices.begin(); dit != dagVertices.end(); ++dit )
    std::cout << *dit << " ";
  std::cout << std::endl;
  if (IsSetCyclic(dagVertices, sourceVertex, path)) {
    std::cout << "the set is cyclic " << std::endl;
    Unflag2(dagEdges);
    dagEdges.clear();
    dagVertices.clear();
    dagEdges.clear();
    return 0;
  }
  // Now, there should only be one dest edge, trace that back to the previous dest edge
  if (optimalVertexPath.find(edges[sinkEdge].src) == optimalVertexPath.end()) {
    // problem, the last edgte in the graph should be in the dag
    std::cout << "ERROR, vertex " << edges[sinkEdge].src << " should be in the dag " << std::endl;
    assert(0);
  }
  
  ssize_t curVertex;
  ssize_t balancedEdge;
  curVertex = edges[sinkEdge].src;
  std::cout << "tracing optimal path from " << curVertex << std::endl;
  ssize_t optEdge;
  while (optimalEdgePath[curVertex] != -1) {
    // Store this edge
    optEdge = optimalEdgePath[curVertex];
    optimalPath.insert(optimalPath.begin(), optEdge);

    // Record that the balanced edge is preferred
    balancedEdge = edges[optEdge].balancedEdge;
    edges[balancedEdge].balPreferred = GraphEdge::Marked;
    
    // Advance to the previous vertex along the optimal path
    curVertex = edges[optEdge].src;
  }
  return 1;
}

ssize_t IntervalGraph::IsSetCyclic(std::set<ssize_t> &vertexSet, 
															 ssize_t curVertex, std::list<ssize_t> &path) {

  std::set<ssize_t>::iterator vertexIt;
  std::list<ssize_t>::iterator pathIt;
  if (std::find(path.begin(), path.end(), curVertex) != path.end()) {
    /*
      std::cout << "cyclic path: to " << curVertex << " ";
      for (pathIt = path.begin(); pathIt != path.end(); pathIt++)
      std::cout << *pathIt << " ";
      std::cout << std::endl;
    */
    return 1;
  }
  else {
    ssize_t outEdge, outEdgeIndex, destVertex;
    path.push_back(curVertex);
    for (outEdgeIndex = vertices[curVertex].FirstOut();
				 outEdgeIndex != vertices[curVertex].EndOut();
				 outEdgeIndex = vertices[curVertex].NextOut(outEdgeIndex)) {
      outEdge = vertices[curVertex].out[outEdgeIndex];
      destVertex = edges[outEdge].dest;
      // Look to see if this vertex is in the set we care about
      vertexIt = vertexSet.find(destVertex);
      if (vertexIt != vertexSet.end()) {
				if (vertices[*vertexIt].flagged2 != GraphVertex::Marked ) {
					if ( IsSetCyclic(vertexSet, destVertex, path)) {
						return 1;
					}
				}
				else {
					//	  std::cout << "not bothering with vertex " << *vertexIt << std::endl;
				}
      }
      else {
				//	std::cout << "did not find " << destVertex << " in dag " << std::endl;
      }
    }

    // Finished checking all paths that exit from this vertex for being cyclic. 
    // that means we don't need to check it again (we might reach this vertex again
    // because there may be undirected cycles in the graph.
    vertexIt = vertexSet.find(curVertex);
    assert(vertexIt != vertexSet.end());
    vertices[*vertexIt].flagged2 = GraphVertex::Marked;
    /*
      std::cout << "path of length: " << path.size() 
      << " to " << curVertex << " is acyclic" << std::endl;
      for (pathIt = path.begin(); pathIt != path.end(); pathIt++)
      std::cout << *pathIt << " ";
      std::cout << std::endl;
    */
    path.pop_back();
    return 0;
  }
}

ssize_t IntervalGraph::PruneDAG(std::set<ssize_t> &dagVertices, 
														std::set<ssize_t> &dagEdges,
														std::vector<ssize_t> &path,
														std::vector<ssize_t> &verticesToDelete,
														std::vector<ssize_t> &edgesToDelete) {
  // remove all edges from a dag unless they are part of 'path'
  
  std::cout << "pruning dag of size " << dagVertices.size() << " " 
						<< dagEdges.size() << " down to " << path.size() + 1 << " vertices " << std::endl;
  std::set<ssize_t>::iterator dagIt, dagEnd;
  std::vector<ssize_t> pathVertices;
  std::vector<ssize_t> sortedPath;
  ssize_t i;
  if (path.size() > 0) {
    pathVertices.push_back(edges[path[0]].src);
  }
  for (i = 0; i < path.size(); i++ ) {
    pathVertices.push_back(edges[path[i]].dest);
    sortedPath.push_back(path[i]);
  }

  // We just need the path to be searchable
  std::sort(pathVertices.begin(), pathVertices.end());
  std::sort(sortedPath.begin(), sortedPath.end());
  
  std::vector<ssize_t> detachedVertices;

  // Remove vertices that are not along the specified path
  for (dagIt = dagVertices.begin(); dagIt != dagVertices.end(); ++dagIt) {
    if (std::binary_search(pathVertices.begin(), pathVertices.end(), *dagIt) == 0) {
      verticesToDelete.push_back(*dagIt);
    }
  }

  // Remove edges that are not along the specified path
  for (dagIt = dagEdges.begin(); dagIt != dagEdges.end(); ++dagIt) {
    if (std::binary_search(sortedPath.begin(), sortedPath.end(), *dagIt) == 0) {
      RemoveEdge(*dagIt, detachedVertices);
      edgesToDelete.push_back(*dagIt);
    }
  }
  
  // we don't really care about the detached vertices, since none of them should
  // remain after removing vertices that are not part of the main path.
}


ssize_t IntervalGraph::VertexOnMaxPath(ssize_t curVertex, std::vector<ssize_t> &optPathTraversals) {
  return (vertices[curVertex].OutDegree() == 1 and
					EdgeOnMaxPath(vertices[curVertex].out[vertices[curVertex].FirstOut()], optPathTraversals));
}

ssize_t IntervalGraph::VertexOnAlternatePath(ssize_t curVertex, std::vector<ssize_t> &optPathTraversals) {
  ssize_t outDegree = vertices[curVertex].OutDegree();
  if (outDegree < 2) 
    return 0;

  ssize_t numMaxPath = 0;
  ssize_t outEdge, outEdgeIndex;
  for (outEdgeIndex = vertices[curVertex].FirstOut();
       outEdgeIndex < vertices[curVertex].EndOut();
       outEdgeIndex = vertices[curVertex].NextOut(outEdgeIndex)) {
    outEdge = vertices[curVertex].out[outEdgeIndex];
    if (optPathTraversals[outEdge] > 0)
      numMaxPath++;
  }
  return (outDegree - numMaxPath == 1);
}

ssize_t IntervalGraph::FindVertexOnMaxPath(ssize_t curVertex, 
																			 std::list<ssize_t> &bestPath,
																			 std::list<ssize_t> &curPath, 
																			 std::vector<ssize_t> &optPathTraversals) {
  
  if (std::find(curPath.begin(), curPath.end(), curVertex) != curPath.end()) {
    std::cout << "found a cycle, not sure what to do in this case " << std::endl;
    std::cout << "work on it if I get here " << std::endl;
    return 0;
  }
  
  curPath.push_back(curVertex);
  ssize_t cyclic = 0;

  if (VertexOnMaxPath(curVertex, optPathTraversals)) {
    // This vertex is along an optimal path, don't continue searching
    if (curPath.size() >= bestPath.size()) {
      // make the current path the best path
      bestPath.clear();
      bestPath.insert(bestPath.end(), curPath.begin(), curPath.end());
      std::cout << "found best path ";
      std::list<ssize_t>::iterator lit;
      for (lit = bestPath.begin(); lit != bestPath.end(); ++lit) {
				std::cout << " " << *lit;
      }
      std::cout << std::endl;
    }
    cyclic = 0;
  }
  else {
    ssize_t outEdge, outEdgeIndex;
    for (outEdgeIndex = vertices[curVertex].FirstOut();
				 outEdgeIndex < vertices[curVertex].EndOut() and cyclic == 0;
				 outEdgeIndex = vertices[curVertex].NextOut(outEdgeIndex)) {
      outEdge = vertices[curVertex].out[outEdgeIndex];
      if (edges[outEdge].traversed == GraphEdge::Marked)
				return 1;

      // Traverse this edge
      edges[outEdge].traversed = GraphEdge::Marked;
      // Use dfs explore.  If any of the routes are cyclic
      // store that.
      if (FindVertexOnMaxPath(edges[outEdge].dest, bestPath, curPath, optPathTraversals) == 0)
				cyclic = 1;
    }
  }
  curPath.pop_back();
  //  std::cout << "done, returning " << !cyclic  << std::endl;
  return (!cyclic);
}


void IntervalGraph::CollectVertices(ssize_t curVertex, ssize_t destVertex, std::set<ssize_t> &dagVertices) {
  std::cout << "collecting from " << curVertex << " to " << destVertex << std::endl;
  if (dagVertices.find(curVertex) == dagVertices.end()) {
    dagVertices.insert(curVertex);
    if (curVertex == destVertex) 
      return;

    ssize_t outEdge, outEdgeIndex;
    for (outEdgeIndex = vertices[curVertex].FirstOut();
				 outEdgeIndex < vertices[curVertex].EndOut();
				 outEdgeIndex = vertices[curVertex].NextOut(outEdgeIndex)) {
      outEdge = vertices[curVertex].out[outEdgeIndex];
      CollectVertices(edges[outEdge].dest, destVertex, dagVertices);
    }
  }
}


ssize_t IntervalGraph::EdgeOnMaxPath(ssize_t edge, std::vector<ssize_t> &optPathTraversals) {
  return optPathTraversals[edge] > 0;
}

ssize_t IntervalGraph::FindRouteToMaxPath(ssize_t curVertex, 
																			std::vector<ssize_t> &optPathTraversals, 
																			std::vector<ssize_t> &vertexPath) {
  std::list<ssize_t> bestPath, curPath;
  if (FindVertexOnMaxPath(curVertex, bestPath, curPath, optPathTraversals)) {
    // the best path stores the longest path to the optimal (anything marked as 
    // passed in optpathtraversals
    vertexPath.clear();
    vertexPath.insert(vertexPath.end(), bestPath.begin(), bestPath.end());
    return 1;
  }
  return 0;
}

void IntervalGraph::CollectEdges(ssize_t lastVertex, std::set<ssize_t> &vertexSet, std::set<ssize_t> &edgeSet) {
  std::set<ssize_t>::iterator vertexIt;
  for (vertexIt = vertexSet.begin(); vertexIt != vertexSet.end(); ++vertexIt) {
    if (*vertexIt == lastVertex)
      continue;
    ssize_t outEdge, outEdgeIndex;
    for (outEdgeIndex = vertices[*vertexIt].FirstOut();
				 outEdgeIndex < vertices[*vertexIt].EndOut();
				 outEdgeIndex = vertices[*vertexIt].NextOut(outEdgeIndex)) {
      outEdge = vertices[*vertexIt].out[outEdgeIndex];
      edgeSet.insert(outEdge);
    }
  }
}

ssize_t IntervalGraph::FindAlternatePaths(ssize_t curVertex,
																			std::vector<ssize_t> &optPathTraversals, 
																			std::vector<ssize_t> &vertexPath,
																			std::vector<ssize_t> &verticesToDelete,
																			std::vector<ssize_t> &edgesToDelete) {

  ssize_t outEdge, outEdgeIndex;
  if (VertexOnAlternatePath(curVertex, optPathTraversals)) {
    // process 
    vertexPath.clear();
    std::cout << "checking for alternate path for " << curVertex << std::endl;
    if (FindRouteToMaxPath(curVertex, optPathTraversals, vertexPath)) {
      // found a route that bypasses the max path, but re-joins it
      // keep searching from ther
      std::cout << "found alternative path ";
      ssize_t v;
      for (v = 0; v < vertexPath.size(); v++) {
				std::cout << " " << vertexPath[v];
      }
      std::cout << std::endl;
      ssize_t last = vertexPath.size()-1;
      if (vertexPath.size() > 1) {

				std::set<ssize_t> bulgeVertices, bulgeEdges;
				std::map<ssize_t, ssize_t> oneToAllMaxPath;
				std::vector<ssize_t> singleMaxPath;
				CollectVertices(vertexPath[0], vertexPath[last], bulgeVertices);

				if (VertexSetContained(bulgeVertices, curVertex)) {

					// Now find a good path through the subset
					std::cout << "The path is contained!!! in a bulge of size " << bulgeVertices.size() << std::endl;

					SingleSourceSubsetMaximumPath(vertices, edges,
																				vertexPath[0], vertexPath[last],
																				bulgeVertices, oneToAllMaxPath);

					// Trace the path from the dest to the source
					TraceOptimalPath(vertices, edges,
													 vertexPath[0], vertexPath[last],
													 oneToAllMaxPath, singleMaxPath);

					std::cout << "retaining path ";
					for (v = 0; v < singleMaxPath.size(); v++ ) {
						std::cout << " " << singleMaxPath[v];
					}
					std::cout << std::endl;
					// Now we need to store the edges that correspond to the dag
					// to check  them if they are suspect (possibly wrong).

					CollectEdges(vertexPath[last], bulgeVertices, bulgeEdges);
					if (IsDAGSuspect(bulgeEdges)) {
						std::cout << "edge set is suspect !!! " << std::endl;
						PruneDAG(bulgeVertices, bulgeEdges, 
										 singleMaxPath, 
										 verticesToDelete, edgesToDelete);
					}
					else {
						ssize_t e;
						std::cout << "edge set is not suspect: ";
						std::set<ssize_t>::iterator beit;
						for (beit = bulgeEdges.begin(); beit != bulgeEdges.end(); ++beit) {
							std::cout << *beit << " ";
						}
						std::cout << std::endl;
					}
				}
				else {
					std::cout << "not contained set of size " << bulgeVertices.size() << ":";
					std::set<ssize_t>::iterator bulgeIt;
					for (bulgeIt= bulgeVertices.begin(); bulgeIt != bulgeVertices.end(); ++bulgeIt) {
						std::cout << " " << *bulgeIt;
					}
					std::cout << std::endl;
				}
				// Continue searching from the end of the 'bulge'
				FindAlternatePaths(vertexPath[last], optPathTraversals, 
													 vertexPath, verticesToDelete, edgesToDelete);
      }
    }
  }
  else {
    // This vertex is not on a single optmial path, search in the graph for more paths.
    // There is no direct path out of this vertex, look for them on other edges
    for (outEdgeIndex = vertices[curVertex].FirstOut();
				 outEdgeIndex < vertices[curVertex].EndOut();
				 outEdgeIndex = vertices[curVertex].NextOut(outEdgeIndex)) {
      outEdge = vertices[curVertex].out[outEdgeIndex];
      
      if (edges[outEdge].traversed == GraphEdge::NotMarked) {
				edges[outEdge].traversed = GraphEdge::Marked;
				//	std::cout << "no alternative path found out of " << curVertex 
				// << " continuing with " << edges[outEdge].dest << std::endl;
				FindAlternatePaths(edges[outEdge].dest, optPathTraversals, vertexPath, verticesToDelete, edgesToDelete);
      }
    }
  }
  // Signal that no optimal path was found from this vertex
  return 1;
}


ssize_t IntervalGraph::FindAlternatePaths() {
  ssize_t v, e;
  // Annotate the sources and sinks in the graph so that they may be 
  // traced back in one-to-all shoretst path traversals.
  
  std::vector<ssize_t> sources, sinks;
  FindSourcesAndSinks(sources, sinks);
  ssize_t source, sink;

  // Count the number of times each highest scoring path
  // passes through an edge
  std::vector<ssize_t> optPathTraversals;
  optPathTraversals.resize(edges.size());
  std::fill(optPathTraversals.begin(), optPathTraversals.end(), 0);

  for (source = 0; source < sources.size(); source++ ) {
    std::cout << "finding paths from source " << sources[source] << std::endl;
    IncrementOptimalPathCount(sources[source], sinks, optPathTraversals);
  }

  Unmark();

  // Starting at each source, do dfs for alternative path.
  std::vector<ssize_t> routeToOptPath;
  std::vector<ssize_t> verticesToDelete, edgesToDelete;
  for (v = 0; v < sources.size(); v++ ) {
    // Find some bulge starting on this vertex
    //    std::cout << "checking out source " << sources[v] << std::endl;
    routeToOptPath.clear();
    FindAlternatePaths(sources[v], optPathTraversals, routeToOptPath, verticesToDelete, edgesToDelete);
  }

  std::cout << "I'm thinking about deleting these vertices: " << std::endl;
  for (v = 0 ; v < verticesToDelete.size(); v++ ) {
    std::cout << verticesToDelete[v]  << std::endl;
  }
  std::cout << "and these edges " << std::endl;
  for (e = 0; e < edgesToDelete.size() ; e++  ) {
    std::cout << edgesToDelete[e] << std::endl;
  }  
  
  Prune(verticesToDelete, edgesToDelete);
  CondenseSimplePaths();
  CheckBalance(edges);
}

ssize_t IntervalGraph::RemoveSuspectBulges() {
  // Precondition: the graph has been processed so that 'suspect' edges are marked
  // These are edges that if there is a bias in sequencing errors, we are likely to mark.

  ssize_t e;  
  ssize_t destVertex, bulgeEnd;
  std::set<ssize_t> dagVertices, dagEdges;
  std::vector<ssize_t> pathThroughDAG;
  std::vector<ssize_t> verticesToRemove, edgesToRemove;
  ssize_t minTrustedLength = 40;
  Untraverse();
  Unmark();
  for (e = 0; e < edges.size(); e++ ) {
    // Don't process the following bulge if this edge is marked
    if (edges[e].traversed == GraphEdge::Marked)
      continue;

    // We will look in bulges to see if they are suspect IF:
    // 1. The vertex at the beginning of the bulge is of in
    //    detree 1 (it *might* not be the beginning of a repeat)
    // 2. The edge preceeding the bulge is long.  We trust long 
    //    edges.  Maybe I should change this to the edge not being suspect. 
    
    destVertex = edges[e].dest;
    // case 1, in degree is 1, and this is a branching vertex (out>1)
    if (vertices[destVertex].InDegree() == 1 and
				vertices[destVertex].OutDegree() > 1) {
      // case 2, we trust the in edge
      if (edges[e].marked != GraphEdge::Marked and
					edges[e].length > minTrustedLength) {
				//	std::cout << "checking edge " << e << " to see if it leads to a bulge " << std::endl;
				dagVertices.clear();
				dagEdges.clear();
				pathThroughDAG.clear();
				// For now use 100 for the edge length.  later look for trusted edges
				if (FindShortContainedDAG(destVertex, minTrustedLength, 
																	dagVertices, dagEdges,
																	pathThroughDAG)) {

					// Some dag (bulge network) exists.  Mark the entire dag as traversed so we
					// dont try and process it again
					TraverseEdgeSet(dagEdges);

					// For now, simply replace that network by it's highest multiplicity path, if the 
					// dag is suspect
					std::set<ssize_t>::iterator dagit;
					std::cout << "dag contains: " << std::endl;
					for (dagit = dagEdges.begin(); dagit != dagEdges.end(); ++dagit) {
						std::cout << *dagit << " ";
					}
					std::cout << std::endl;
					std::cout << "the optimal path contains: " << std::endl;
					ssize_t i;
					for (i = 0; i < pathThroughDAG.size(); i++) { 
						std::cout << pathThroughDAG[i] << " ";
					}
					std::cout << std::endl;
					if (IsDAGSuspect(dagEdges)) {
						// ok, the dag is suspect.  Remove all the edges that are not on the 
						// optimal path
	    
						// sanity check... is there an optimal path?
						assert(pathThroughDAG.size() > 0);

						PruneDAG(dagVertices, dagEdges, pathThroughDAG, verticesToRemove, edgesToRemove);
					}
				}
      }
    }
  }
  std::cout << "Removing suspect bulges found " << verticesToRemove.size() << " vertices "
						<< " and " << edgesToRemove.size() << " edges" << std::endl;
  Prune(verticesToRemove, edgesToRemove);
  CondenseSimplePaths();
}

ssize_t IntervalGraph::IsDAGSuspect(std::set<ssize_t> &dagEdges) {
  std::set<ssize_t>::iterator dagIt, dagEnd;
  ssize_t suspect = 0;
  for (dagIt = dagEdges.begin(); dagIt != dagEdges.end(); dagIt++ ) {
    if (edges[*dagIt].suspect == GraphEdge::Marked) {
      return 1;
    }
  }

  // There is no suspect edge in this dag, say that
  // it is not suspect
  return 0;
}
