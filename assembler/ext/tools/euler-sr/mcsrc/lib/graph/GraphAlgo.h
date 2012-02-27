/***************************************************************************
 * Title:          GraphAlgo.h 
 * Author:         Mark Chaisson, Glenn Tesler
 * Created:        2007
 * Last modified:  01/19/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef GRAPH_ALGO_H_
#define GRAPH_ALGO_H_

#include <vector>
#include <iostream>
#include <set>
#include <map>
#include "rbtree.h"
#include "utils.h"

class GraphVertex {
public:
	// HAS BITFIELD
  enum Mark_e { Marked, NotMarked};
  unsigned char marked: 1;

  void Unmark() { marked = NotMarked;}
  virtual ssize_t FirstOut() = 0;
  virtual ssize_t FirstIn() = 0;
  virtual ssize_t EndOut() = 0;
  virtual ssize_t EndIn() = 0;
  virtual ssize_t NextOut(ssize_t curEdge) = 0;
  virtual ssize_t NextIn(ssize_t curEdge) = 0;
  virtual ssize_t InDegree() = 0;
  virtual ssize_t OutDegree() = 0;
};

class GraphEdge {
public:
	// HAS BITFIELD
  enum Mark_e { Marked, NotMarked};
  unsigned char marked: 1;
  unsigned char mst: 1;

  void Unmark() { marked = NotMarked;}
};


template<typename V, typename E, typename Funct>
void TraverseRandomAccess(std::vector<V> &vertices, std::vector<E> &edges,
			  Funct &funct) {
  ssize_t v;
  for (v = 0; v < vertices.size(); v++) {
    funct(v);
  }
}

template<typename V, typename E, typename Funct> 
void TraverseDFS(std::vector<V> &vertices, std::vector<E> &edges, ssize_t vertexIndex, Funct &funct) {
  ssize_t edgeIndex;
  ssize_t e;
  if (vertices[vertexIndex].marked == GraphVertex::Marked) 
    return;
  
  // the vertex is new to the dfs, run 'f' on it.
  funct(vertexIndex);

  // Don't want to visit this vetex again, so mark it
  vertices[vertexIndex].marked = GraphVertex::Marked;

  // search back
  ssize_t srcIndex;
  for (e = vertices[vertexIndex].FirstIn(); 
       e < vertices[vertexIndex].EndIn(); 
       e = vertices[vertexIndex].NextIn(e)) {
    edgeIndex = vertices[vertexIndex].in[e];
    srcIndex = edges[edgeIndex].src;
    if (vertices[srcIndex].marked == GraphVertex::NotMarked) {
      TraverseDFS(vertices, edges, srcIndex, funct);
    }
  }

  // search forward
  ssize_t destIndex;
  for (e = vertices[vertexIndex].FirstOut(); 
       e < vertices[vertexIndex].EndOut(); 
       e = vertices[vertexIndex].NextOut(e)) {
    edgeIndex = vertices[vertexIndex].out[e];
    destIndex = edges[edgeIndex].dest;
    if ( vertices[destIndex].marked == GraphVertex::NotMarked ) {
      TraverseDFS(vertices, edges, destIndex, funct);
    }
  }
}

template<typename V, typename E, typename Funct> 
void TraverseGraphDFS(std::vector<V> &vertices, std::vector<E> &edges, Funct &funct) {
  ssize_t v;
  for (v = 0; v < vertices.size(); v++) {
    if ( vertices[v].marked == GraphVertex::NotMarked ) {
      TraverseDFS(vertices, edges, v, funct);
    }
  }
}

template<typename T>
void Unmark(std::vector<T> &elements) {
  ssize_t i;
  for (i = 0; i < elements.size(); i++) {
    elements[i].Unmark();
  }
}


template<typename T, typename S>
class IndexedScorePriorityQueue {
public:
  S* scores;
  std::set<T> indices;
};

template<typename T> 
class SearchablePriorityQueue {
public:
  T min;
  std::set<T> heap;
  ssize_t size() {
    return heap.size();
  }
  void push(const T &value) {
    heap.insert(value);
  }
  T pop() {
    // adjust the bounds if T is the min value
    typename std::set<T>::iterator last, nextToLast;

    last = heap.begin();
    T value = *last;
    heap.erase(last);
    return value;
  }

  ssize_t find(const T &key) {
    typename std::set<T>::iterator heapIt;
    for (heapIt = heap.begin(); heapIt != heap.end(); ++heapIt) {
      if (*heapIt == key)
				return 1;
    }
    return 0;
  }
  T&operator[](const T &key) {
    ssize_t i;
    typename std::set<T>::iterator it, end;
    end = heap.end();
    for (it = heap.begin(); it != heap.end(); ++it) {
      if (*it == key) {
				return (T&) *it;
      }
    }
    assert(it != heap.end());
  }
};

template<typename V>
class ScoredVertex {
public:
  // This class adds a score component to a vertex
  // It's convenient when only a small component of a graph will be
  // scored, so that the entire graph does not need to have an extra field 
  // added.
  
  V vertex; // a way of indexing a vertex
  ssize_t score; // score of the vertex
  ScoredVertex() {
    vertex = -1;
    score = 0;
  }

  ScoredVertex(const ScoredVertex<V> &v) {
    *this = v;
  }
  ScoredVertex(ssize_t v, ssize_t s) {
    vertex = v;
    score  = s;
  }
  ScoredVertex<V> &operator=(const ScoredVertex<V> &v) {
		if (this != &v) {
			vertex = v.vertex;
			score  = v.score;
		}
    return *this;
  };
  bool operator>(const ScoredVertex<V> &v) const {
    return score > v.score;
  }
  bool operator<(const ScoredVertex<V> &v) const {
    return score < v.score;
  }
  bool operator==(const ScoredVertex<V> &v) const {
    return vertex == v.vertex;
  }
};
template<typename T>
std::ostream & operator<<(std::ostream& out, const ScoredVertex<T> &sv) {
  out << sv.vertex << " " << sv.score << " ";
  return out;
}

/*
  template<typename T>
  class CompareScoredVertex {
  public:
  bool operator<(const ScoredVertex<T> &A, const ScoredVertex<T> &B) const {
    return A.score < B.score;
  }
  };
*/

typedef RBTreeNode<ScoredVertex<ssize_t> > HeapScoredVertex;


template<typename V, typename E>
void BellmanFordMaximumPath(std::vector<V> &vertices,
														std::vector<E> &edges,
														ssize_t source, std::vector<ssize_t> &predecessor) {

  std::cout << "THIS IS NOT FINISHED NOR TESTED" << std::endl;
  assert(0);
  ssize_t v, ei, e;

  if (edges.size() < 1) return;
  // Preprocessing step, find the maximum cost edge in the graph
  ssize_t maxCost = edges[0].Cost();
  for (e = 0; e < edges.size(); e++ ) {
    if (maxCost < edges[e].Cost())
      maxCost = edges[e].Cost();
  }
  
  ssize_t infinity = 999999999;
  std::vector<ssize_t> distances;
  distances.resize(vertices.size());
  predecessor.resize(vertices.size());
  std::fill(distances.begin(), distances.end(), infinity);
  std::fill(predecessor.begin(), predecessor.end(), -1);

  distances[source] = 0;
  ssize_t dest, src;
  for (v = 0; v < vertices.size(); v++ ) {
    for (e = 0; e < edges.size(); e++ ) {
      dest = edges[e].dest;
      src  = edges[e].src;

      if (distances[dest] > distances[src] + (maxCost - edges[e].Cost())) {
				distances[dest] = distances[src] + (maxCost - edges[e].Cost());
				predecessor[dest] = ei;
      }
    }
  }
  
  // Check for negative weight cycles
  for (e = 0; e < edges.size(); e++ ) {
    src = edges[e].src;
    dest = edges[e].dest;
    if (distances[dest] > distances[src] + (maxCost - edges[e].Cost())) {
      std::cout << "Error, negative weight cycle on edge " << edges[e] << std::endl;
      std::cout << distances[dest] << " " << distances[src] << " " << maxCost - edges[e].Cost() << std::endl;
      assert(0);
    }
  }
}


template<typename V, typename E>
void SingleSourceSubsetMaximumPath(std::vector<V> &vertices, std::vector<E> &edges, 
																	 ssize_t sourceVertex, ssize_t destVertex,
																	 std::set<ssize_t> vertexSubset, std::map<ssize_t,ssize_t> &shortestPathEdges) {
  ssize_t maxCost = 0;
  ssize_t v, e;
  std::set<ssize_t>::iterator vertexIt;
  ssize_t outEdge, outEdgeIndex;
  
  // Find the maximum cost vertex, but over only a subset of vertices
  for (vertexIt = vertexSubset.begin();
       vertexIt != vertexSubset.end();
       ++vertexIt) {
    if (*vertexIt != destVertex) {
      for (outEdgeIndex = vertices[*vertexIt].FirstOut();
					 outEdgeIndex != vertices[*vertexIt].EndOut();
					 outEdgeIndex = vertices[*vertexIt].NextOut(outEdgeIndex)) {
				outEdge = vertices[*vertexIt].out[outEdgeIndex];
				
				if (edges[outEdge].Cost() >maxCost) 
					maxCost = edges[outEdge].Cost();
      }
    }
  }
	
  // For now, the code is more-or-less copied from the 
  // singtle source max path. It just searches throug the vertex subset alone

  // prepare references to scores
  std::map<ssize_t, HeapScoredVertex*> references;

  // prepare the priority queue
  RBTree<ScoredVertex<ssize_t> > queue;
  HeapScoredVertex *minVertex;
  
  references[sourceVertex] = queue.Insert(ScoredVertex<ssize_t>(sourceVertex, 0));
  
  // Find One-to-all shortest path
  
  ssize_t dest;
  ssize_t first = 1;
  ScoredVertex<ssize_t> closestVertex;
  ssize_t vertex, closestScore;
  ssize_t inEdgeIndex;
  while (queue.size() > 0) {
    // fetch the highest 'closest' vertex

    if (queue.Pop(closestVertex) == 0) {
      std::cout << "INTERNAL ERROR: the size of the queue was " << queue.size() 
								<< " but nothing could be popped" << std::endl;
      assert(0);
    }
    if (closestVertex.vertex  == destVertex)
      return;

    // Relax vertices out of 'vertex'
    for (outEdgeIndex = vertices[closestVertex.vertex].FirstOut();
				 outEdgeIndex != vertices[closestVertex.vertex].EndOut();
				 outEdgeIndex = vertices[closestVertex.vertex].NextOut(outEdgeIndex)) {
      outEdge = vertices[closestVertex.vertex].out[outEdgeIndex];
      dest = edges[outEdge].dest;
      // If this vertex is already in the priority queue, 
      // We need to store a new score if the dest hasn't had a score
      // computed for it, or the score is worse than the new score
      if (references.find(dest) == references.end() or
					(references[dest]->data.score > closestVertex.score + (maxCost - edges[outEdge].Cost()))) {
				// Need to update the score on the dest vertex
				// First remove the old score from the priority queue

				if (references.find(dest) != references.end()) {
					queue.Delete(references[dest]);
				}

				// Update the optimal score in the queue
				ssize_t destScore;
				destScore        = closestVertex.score + (maxCost - edges[outEdge].Cost());
				ssize_t prevSize     = queue.size();
				references[dest] = queue.Insert(ScoredVertex<ssize_t>(dest, destScore));

				// Store the path information for the optimal path
				inEdgeIndex = vertices[dest].LookupInIndex(outEdge);
				assert(inEdgeIndex >= 0);
				shortestPathEdges[dest] = inEdgeIndex;
      }
    }
  }
}



template<typename V, typename E>
void SingleSourceMaximumPath(std::vector<V> &vertices, std::vector<E> &edges, 
														 ssize_t source, std::vector<ssize_t> &shortestPathEdges, ssize_t destVertex=-1) {

  ssize_t v, e;

  if (edges.size() < 1) return;
  // Preprocessing step, find the maximum cost edge in the graph
  ssize_t maxCost = edges[0].Cost();
  for (e = 0; e < edges.size(); e++ ) {
    if (maxCost < edges[e].Cost())
      maxCost = edges[e].Cost();
  }

  // Set all edges to not traversed
  shortestPathEdges.resize(vertices.size());
  for (v =0 ; v < vertices.size(); v++ ) shortestPathEdges[v] = -1;
  
  // prepare references to scores
  std::vector<HeapScoredVertex*> references;

#ifdef _DEBUG_
  std::vector<char> visited;
  visited.resize(vertices.size());
  std::fill(visited.begin(), visited.end(), 0);
#endif

  references.resize(vertices.size());
  for (v = 0; v < vertices.size(); v++) references[v] = NULL;
  
  // prepare the priority queue
  RBTree<ScoredVertex<ssize_t> > queue;
  //UNUSED// HeapScoredVertex *minVertex;
  
  references[source] = queue.Insert(ScoredVertex<ssize_t>(source, 0));
  
  // Find One-to-all shortest path
  
  ssize_t outEdge, outEdgeIndex, dest;
  //UNUSED// int first = 1;
  ScoredVertex<ssize_t> closestVertex;
  //UNUSED// int vertex, closestScore;
  ssize_t inEdgeIndex;
  while (queue.size() > 0) {
    // fetch the highest 'closest' vertex
    /*
      std::cout << "cur queue " << queue.size() << std::endl;
      std::cout << queue << std::endl;
    */
    if (queue.Pop(closestVertex) == 0) {
      std::cout << "INTERNAL ERROR: the size of the queue was " << queue.size() 
								<< " but nothing could be popped" << std::endl;
      assert(0);
    }

		//      std::cout << "popped vertex: " << closestVertex.vertex << std::endl;
		//      std::cout << queue << std::endl;

#ifdef _DEBUG_
    if (visited[closestVertex.vertex] == 1) {
      std::cout << "ERROR, re-visiting " << closestVertex.vertex << std::endl;
      assert(0);
    }
    visited[closestVertex.vertex] = 1;
#endif
    // Relax vertices out of 'vertex'
    for (outEdgeIndex = vertices[closestVertex.vertex].FirstOut();
				 outEdgeIndex != vertices[closestVertex.vertex].EndOut();
				 outEdgeIndex = vertices[closestVertex.vertex].NextOut(outEdgeIndex)) {
      outEdge = vertices[closestVertex.vertex].out[outEdgeIndex];
      dest = edges[outEdge].dest;
      // If this vertex is already in the priority queue, 
      // We need to store a new score if the dest hasn't had a score
      // computed for it, or the score is worse than the new score
      if (references[dest] == NULL or 
					(references[dest] != NULL and
					 references[dest]->data.score > closestVertex.score + (maxCost - edges[outEdge].Cost()))) {
				// Need to update the score on the dest vertex
				// First remove the old score from the priority queue
	
				if (references[dest] != NULL) {
					//	  std::cout << "Deleting " << dest << " from queue " << std::endl;
					queue.Delete(references[dest]);
				}

				// Update the optimal score in the queue
				ssize_t destScore;
				destScore = closestVertex.score + (maxCost - edges[outEdge].Cost());
				//	std::cout << "Added  " << dest << " to the queue " << std::endl;
				references[dest] = queue.Insert(ScoredVertex<ssize_t>(dest, destScore));
				// Store the path information for the optimal path
				inEdgeIndex = vertices[dest].LookupInIndex(outEdge);
				assert(inEdgeIndex >= 0);
				shortestPathEdges[dest] = inEdgeIndex;
      }
    }
  }
}

template<typename V, typename E, typename List>
void TraceOptimalPath(std::vector<V> &vertices,
											std::vector<E> &edges,
											ssize_t sourceVertex, ssize_t destVertex,
											List pathList,
											std::vector<ssize_t> &edgePath) {

  assert(destVertex >= 0);
  ssize_t curVertex = destVertex;
  ssize_t edge;
  //  assert(std::find(pathList.begin(), pathList.end(), curVertex) != pathList.end());
  while(curVertex != sourceVertex) {
    edge = vertices[curVertex].in[pathList[curVertex]];
    //    std::cout << "edge: " << edge << " " << startToAll[curVertex] << std::endl;
    edgePath.insert(edgePath.begin(), edge);
    curVertex = edges[edge].src;
    
    // Some sanity checks. This will bomb if we are in an ininite loop
    assert(edgePath.size() < edges.size());
    //    assert(std::find(pathList.begin(), pathList.end(), curVertex) != pathList.end());
  }
}


template<typename V, typename E>
ssize_t StoreFinishingTimes(std::vector<V> &vertices, std::vector<E> &edges,
												std::vector<ssize_t> &finishingTimes, std::vector<ssize_t> &finishingOrder) {
  finishingTimes.clear();
  finishingTimes.resize(vertices.size());
  finishingOrder.clear();
  finishingOrder.resize(vertices.size());
  Pause();
  Unmark(vertices);
  ssize_t v;
  ssize_t time = 0;
  ssize_t lastVertex = -1;
  for (v = 0; v < vertices.size(); v++ ) {
    if (vertices[v].marked == GraphVertex::NotMarked) {
      lastVertex = StoreFinishingTimes(vertices, edges, 
																			 finishingTimes, finishingOrder,
																			 v, time, lastVertex);
    }
  }
  return lastVertex;
}

template<typename V, typename E>
ssize_t StoreFinishingTimes(std::vector<V> &vertices, std::vector<E> &edges,
												std::vector<ssize_t> &f, std::vector<ssize_t> &finishingOrder, 
												ssize_t curVertex, ssize_t &curTime, ssize_t &previous) {
  // OK, this is now written to be particular to the method for fining
  // strongly connected components, but it should be ok.
  
  ssize_t outEdge, outEdgeIndex;
  //UNUSED// int lastVertex = -1;
  ssize_t destVertex = -1;
  assert(vertices[curVertex].marked != GraphVertex::Marked);
  vertices[curVertex].marked = GraphVertex::Marked;
  for (outEdgeIndex = vertices[curVertex].FirstOut();
       outEdgeIndex < vertices[curVertex].EndOut();
       outEdgeIndex = vertices[curVertex].NextOut(outEdgeIndex)) {
    outEdge    = vertices[curVertex].out[outEdgeIndex];
    destVertex = edges[outEdge].dest;
    // Only traverse edges that are part of the higest scoring tree
    // These are marked.
    // Only visit vertices that are not yet visited
    if (edges[outEdge].marked == GraphEdge::Marked and
				vertices[destVertex].marked != GraphVertex::Marked) {
      StoreFinishingTimes(vertices, edges, f, finishingOrder, destVertex, curTime, previous);
    }
  }
  ++curTime;
  finishingOrder[curVertex] = previous;
	//UNUSED?//  vertices[curVertex].length = curTime;
  f[curVertex]    = curTime;
  previous        = curVertex;
  return curVertex;
}


template<typename V, typename E>
ssize_t DFSBackwards(std::vector<V> &vertices, std::vector<E> &edges,
								 ssize_t curVertex, std::vector<ssize_t> &traced, std::vector<ssize_t> &finishingTimes) {
  
  ssize_t inEdge, inEdgeIndex;
  ssize_t srcVertex;
  ssize_t last = -1;
  assert(vertices[curVertex].marked != GraphVertex::Marked);
  if (vertices[curVertex].marked != GraphVertex::Marked) {
    traced.push_back(curVertex);
    assert(curVertex >= 0 and curVertex < vertices.size());
    vertices[curVertex].marked = GraphVertex::Marked;
    /*
			std::cout << "dfsb " << curVertex << " " << finishingTimes[curVertex] << std::endl;
    */
    for (inEdgeIndex = vertices[curVertex].FirstIn();
				 inEdgeIndex < vertices[curVertex].EndIn();
				 inEdgeIndex = vertices[curVertex].NextIn(inEdgeIndex)) {
      inEdge = vertices[curVertex].in[inEdgeIndex];
      srcVertex = edges[inEdge].src;
      assert(srcVertex >= 0 and srcVertex < vertices.size());
      /*
				std::cout << inEdge << " " << (ssize_t) edges[inEdge].marked <<  " " 
				<< (ssize_t) vertices[srcVertex].marked << " " << (ssize_t) GraphVertex::Marked <<  std::endl;
      */
      if (vertices[srcVertex].marked == GraphVertex::NotMarked and
					edges[inEdge].marked  == GraphEdge::Marked) {
				last = DFSBackwards(vertices, edges, srcVertex, traced, finishingTimes);
      }
    }
  }
  if (last == -1)
    return curVertex;
  else
    return last;
}

template<typename V, typename E>
ssize_t FindStronglyConnectedComponents(std::vector<V> &vertices, std::vector<E> &edges,
																		std::vector<std::vector<ssize_t > > &sccs) {
  std::vector<ssize_t> finishingTimes;
  std::vector<ssize_t> finishingOrder;

  ssize_t curVertex, lastVertex;
  
  lastVertex = StoreFinishingTimes(vertices, edges,
																	 finishingTimes, finishingOrder);

  Unmark(vertices);
  ssize_t nProcessed= 0;
  std::vector<ssize_t> scc;
  
  curVertex = lastVertex;
  /*
		ssize_t v;
		std::cout << "graph traversal result:" << std::endl;
		for (v = 0 ; v < vertices.size(); v++ ) {
		std::cout << v << " " << finishingTimes[v] << " " << finishingOrder[v] << std::endl;
		}
  */
  ssize_t prevOrder = finishingTimes[curVertex] + 1;
  while (finishingOrder[curVertex] != -1) {
    /*
      std::cout << "backwards dfs for " << curVertex << " f: " << finishingTimes[curVertex] 
      << " next: " << finishingOrder[curVertex] << std::endl;
    */
    assert(prevOrder > finishingTimes[curVertex]);
    prevOrder = finishingTimes[curVertex];
    // All sorts of cases to do sanity checks onl
    // First, we should always be starting the dfs on a vertex that has not 
    // been touched
    assert(curVertex < vertices.size() and curVertex >= 0);
    assert(vertices[curVertex].marked == GraphVertex::NotMarked);


    // Next, the number of vertices processed can't be more than the number 
    // of vertices in the graph.  If it is, this is probably in some infinite loop.
    assert(nProcessed < vertices.size());

    scc.clear();
    lastVertex = DFSBackwards(vertices, edges, curVertex, scc, finishingTimes);
    

    // Finally, the dfs should have processed at least one vertex
    // and, to make sure the last assignment is workin properly, 
    // the last index returned should have been the last one added to the list as well.
    assert(scc.size() > 0);
    assert(lastVertex == scc[scc.size()-1]);

    // Move curVertex on to the next unprocessed vertex.
    // This should be the vertex preceeding the earliest finish
    /*
      ssize_t c;
      ssize_t soonestFinish = -1;
      ssize_t soonestFinishIndex = -1;
    */
    
    while (vertices[curVertex].marked == GraphVertex::Marked and
					 finishingTimes[curVertex] >= 0 and 
					 finishingOrder[curVertex] >= 0) {
      curVertex = finishingOrder[curVertex];
    }
    //    curVertex = finishingOrder[scc[soonestFinishIndex]];
    /*
      std::cout << "dfs started with " << prevOrder << " ended " << finishingTimes[curVertex] 
      << " " << scc.size() <<  std::endl;
    */
    if (scc.size() > 1) {
      sccs.push_back(scc);
    }
  }
  return sccs.size();
}


// Return index of balancing vertex, or -1 if can't determine it.
// Optionally return e, an edge with the vertex
// 
// Since there isn't a field for it, the way to find it is
// to look at the dual of an edge on the vertex.
// This won't work for isolated vertices, however.
// Return -1 if can't determine the balancing vertex.
template<typename V, typename E>
	ssize_t LookupBalancedVertex(std::vector<V> &vertices,
															 std::vector<E> &edges,
															 ssize_t vertex) {

	// Check outgoing edges
	for (ssize_t outEdgeIndex = vertices[vertex].FirstOut();
			 outEdgeIndex != vertices[vertex].EndOut();
			 outEdgeIndex = vertices[vertex].NextOut(outEdgeIndex)) {
		ssize_t outEdge = vertices[vertex].out[outEdgeIndex];

		// not an edge
		if (edges[outEdge].dest == -1)
			continue;

		// found an edge.
		// vertex is the src on outEdge,
		// and balVertex is the dest on balEdge
		ssize_t balEdge = edges[outEdge].balancedEdge;
		ssize_t balVertex = edges[balEdge].dest;
		return balVertex;
	}

	// Check incoming edges
	for (ssize_t inEdgeIndex = vertices[vertex].FirstIn();
			 inEdgeIndex != vertices[vertex].EndIn();
			 inEdgeIndex = vertices[vertex].NextIn(inEdgeIndex)) {
		ssize_t inEdge = vertices[vertex].in[inEdgeIndex];

		// not an edge
		if (edges[inEdge].src == -1)
			continue;

		// found an edge.
		// vertex is the dest on outEdge,
		// and balVertex is the src on balEdge
		ssize_t balEdge = edges[inEdge].balancedEdge;
		ssize_t balVertex = edges[balEdge].src;
		return balVertex;
	}

	return -1;
}


#endif
