/***************************************************************************
 * Title:          MSTAlgo.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/09/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MST_ALGO_H_
#define MST_ALGO_H_

#include "DisjointSet.h"
#include "GraphAlgo.h"
#include <vector>
#include <queue>
#include <set>
namespace GraphAlgo {
  enum MSTEdge { MSTIn, MSTOut };

  template<typename Edge_t> 
		class CompareEdgePtrs {
		public:
    ssize_t operator()(Edge_t* ptr1, Edge_t* ptr2) {
      return (ptr1->Cost() < ptr2->Cost());
    }
  };
  
  
  template<typename Vertex_t, typename Edge_t> 
		void CalcMST(std::vector<Vertex_t> &vertices, std::vector<Edge_t> &edges) {
    Unmark(vertices);
    Unmark(edges);
    
    std::vector<Edge_t*> edgePtrs;
    edgePtrs.resize(edges.size());
    ssize_t e;
    for (e = 0; e < edges.size(); e++ )
      edgePtrs[e] = &(edges[e]);

    CompareEdgePtrs<Edge_t> compareEdgePtrs;
    std::sort(edgePtrs.begin(), edgePtrs.end(), compareEdgePtrs);

    std::vector<DisjointSet::DJVertex<Vertex_t > > djForest;
    DisjointSet::CreateForest(vertices, djForest);
  
    ssize_t srcVertex, destVertex;
    Edge_t *edge;
    for (e = 0; e < edgePtrs.size(); e++ ) {
      edge = edgePtrs[e];
      srcVertex = edge->src;
      destVertex = edge->dest;
      if (DisjointSet::Find(&(djForest[srcVertex])) != 
					DisjointSet::Find(&(djForest[destVertex]))) {
				DisjointSet::Union(&(djForest[srcVertex]), 
													 &(djForest[destVertex]));
				edge->mst = MSTIn;
      }
      else {
				
      }
    }
  }

  template<typename Vertex_t, typename Edge_t>
		void MarkLowestCostInEdges(std::vector<Vertex_t> &vertices, std::vector<Edge_t> &edges) {
    ssize_t v;
    ssize_t inEdge, inEdgeIndex, lowestInEdgeIndex;
    ssize_t minInEdgeCost;

    for (v = 0; v < vertices.size(); v++) {
      minInEdgeCost = 99999999;
      lowestInEdgeIndex = -1;
      for (inEdgeIndex = vertices[v].FirstIn();
					 inEdgeIndex < vertices[v].EndIn();
					 inEdgeIndex = vertices[v].NextIn(inEdgeIndex)) {
				inEdge = vertices[v].in[inEdgeIndex];
				if (edges[inEdge].Cost() < minInEdgeCost) {
					minInEdgeCost = edges[inEdge].Cost();
					lowestInEdgeIndex = inEdgeIndex;
				}
      }
      assert(vertices[v].InDegree() == 0 or lowestInEdgeIndex >= 0);
      if (lowestInEdgeIndex >= 0) {
				//	std::cout << "vertex: " << v << "low edge: " << vertices[v].in[lowestInEdgeIndex] << std::endl;
				edges[vertices[v].in[lowestInEdgeIndex]].marked = GraphEdge::Marked;
				// make this red when we print 
				edges[vertices[v].in[lowestInEdgeIndex]].flagged = GraphEdge::Marked;
      }
    }
  }

  template<typename E>
		class CollapsedVertex {
		public:
    ssize_t proxyVertex;
		// vertexToEdge is a convenient way of representing the cycle.  
		// Maybe I should just rename it 'cycle', but oh well. It's a map from
		// vertices in the cycle to the in edge that is part of the cycle into
		// each cycle vertex.
    std::map<ssize_t, ssize_t> vertexToEdge;
		std::set<ssize_t> edges;
		std::set<ssize_t> vertices;
    std::vector<std::pair<ssize_t, ssize_t> > inEdges, outEdges;
		//    std::vector<E*> deletedEdges;
    std::vector<ssize_t> srcIndices;
		std::vector<ssize_t> srcVertex, destVertex;
    std::vector<ssize_t> destIndices;
    std::vector<ssize_t> edgeIndices;
		std::vector<unsigned char> optimal;

    void RecordRemovedEdge(ssize_t edgeIndex, ssize_t src, ssize_t srcIndex, ssize_t dest, ssize_t destIndex, unsigned char marked) {
			srcVertex.push_back(src);
			destVertex.push_back(dest);
			//      deletedEdges.push_back(edge);
      edgeIndices.push_back(edgeIndex);
      srcIndices.push_back(srcIndex);
      destIndices.push_back(destIndex);
			optimal.push_back(marked);
    }
		ssize_t EdgeStartsInCycle(E& edge) {
			ssize_t srcVertex = edge.src;
			return (vertexToEdge.find(srcVertex) != vertexToEdge.end());
		}
		ssize_t EdgeEndsInCycle(E& edge) {
			ssize_t destVertex = edge.dest;
			return (vertexToEdge.find(destVertex) != vertexToEdge.end());
		}
		ssize_t VertexInCycle(ssize_t v) {
			return (vertices.find(v) != vertices.end());
		}
		ssize_t EdgeInCycle(ssize_t e) {
			return (edges.find(e) != edges.end());
		}
  };


  template<typename Vertex_t, typename Edge_t> 
		void ContractCycle(std::vector<Vertex_t> &vertices, std::vector<Edge_t> &edges, 
											 std::vector<ssize_t> &cycle,
											 CollapsedVertex<Edge_t> &collapsedVertex) {

    ssize_t c, v, e;
    ssize_t edge;
    ssize_t cycleCost = 0;
		// ************************************************************
    // First step.
		// Record the edges that are in the cycle, and
    // calculate the cost of the cycle (the sum of the costs of all edges
		// in the cycle.
		// ************************************************************
		
		ssize_t cycleVertex;
    for (v = 0; v < cycle.size(); v++ ) {
			collapsedVertex.vertices.insert(cycle[v]);
      for (e = 0; e < vertices[cycle[v]].EndIn(); e++ ) {
				cycleVertex = cycle[v];
				edge = vertices[cycleVertex].in[e];
				if (edge >= 0 and edges[edge].marked == GraphEdge::Marked) {
					//					std::cout << "v to e: " << cycle[v] << " " << edge << std::endl;
					collapsedVertex.vertexToEdge[cycleVertex] = edge;
					collapsedVertex.edges.insert(edge);
					//					std::cout << "cycle: " << cycle[v] << " " << edge << std::endl;
					cycleCost += edges[edge].Cost();
				}
      }
    }
    /*
      std::cout << "The cost of the cycle is: " << cycleCost << std::endl;
    */
    // Update the cost of the out edges
    ssize_t inEdge, inEdgeIndex;
    //UNUSED// int outEdge;
		ssize_t outEdgeIndex;
    ssize_t destVertex;
    ssize_t sourceVertex;
    //UNUSED// Edge_t *deletedEdge;
    ssize_t removedSrcIndex, removedDestIndex;
    // The proxy vertex will be the vertex that will be pointed to by 
    // edges that point to the consolidated vertex.
    // This way no extra space is needed.
    // I'm not sure what will happen if a consolidated vertex gets used twice.
    // That's to think about when the first part of this works.

		// ************************************************************
		// Step 2.
		// Update the costs of all edges that enter the cycle.
		// The cost is reflected to be the maximum cost of traversing the 
		// entire cycle from that edge.
		//

    //UNUSED// int newEdgeNeeded;
    // Recompute the costs of the edges that reach the cycle
    for (c = 0; c < cycle.size(); c++ ) {
      cycleVertex = cycle[c];
      // Update the cost and destination of all the edges in the graph to 
      // point to the proxy vertex
      ssize_t cycleInEdge;
      for (inEdgeIndex = vertices[cycleVertex].FirstIn();
					 inEdgeIndex < vertices[cycleVertex].EndIn();
					 inEdgeIndex = vertices[cycleVertex].NextIn(inEdgeIndex)) {
				cycleInEdge = vertices[cycleVertex].in[inEdgeIndex];
				sourceVertex = edges[cycleInEdge].src;
				// Make sure the edge isn't one of the cycle edges
				if (!collapsedVertex.EdgeStartsInCycle(edges[cycleInEdge])) {
					//				if (collapsedVertex.vertexToEdge.find(sourceVertex) == collapsedVertex.vertexToEdge.end()) {
					// The in edge is not part of the cycle
					//UNUSED// int prevCost = edges[cycleInEdge].Cost();
					if (cycleCost - edges[collapsedVertex.vertexToEdge[cycleVertex]].Cost() >= 0) {
						/*
						std::cout << "warning, this may not be good " << std::endl;
						std::cout << cycleCost - edges[collapsedVertex.vertexToEdge[cycleVertex]].Cost() << std::endl;
						*/
					}
					ssize_t curCost = edges[cycleInEdge].Cost();
					edges[cycleInEdge].Cost(curCost + cycleCost - edges[collapsedVertex.vertexToEdge[cycleVertex]].Cost());
					/*
						std::cout << "updated cost of " << cycleInEdge << " " << prevCost 
						<< " -> " << edges[cycleInEdge].Cost() << std::endl;
					*/
				}
      }
    }


		// ************************************************************
		// Step 3.
		// Route edges from vertices outside the cycle into the collapsed
		// vertex.
		// Let v be a vertex that has from 1 to K edges going into the cycle.
		// Pick the highest scoring edge  (v,u*), and move it from the cycle 
		// vertex to the proxy vertex.  
		// Delete all other lower scoring edges (but record them in the collapsed
		// vertex.

    ssize_t proxyVertex;
		proxyVertex = cycle[0];
		//		std::cout << "collapsing a cycle of size " << cycle.size() << " to " << proxyVertex << std::endl;
		
		std::vector<ssize_t> newProxyInEdges;

    collapsedVertex.proxyVertex = proxyVertex;

    ssize_t endCycleVertexIn;
    // Now, given a vertex outside the cycle, v:
    //   Find the edge (v, u*) that has the lowest cost for entering the cycle.
    //   Delete all edges (v,u) that enter the cycle
    //   Add  (v,p) to the graph, so that the proxy is pointed to by v.
    for (c = 0; c < cycle.size(); c++ ) {
      cycleVertex = cycle[c];

      // First find a vertex that has the cycle as the destination of
      // one of it's out edges.  This vertex will be removed from pointing to the cycle

      // Update the cost and destination of all the edges in the graph to 
      // point to the proxy vertex
      ssize_t cycleInEdge;

      endCycleVertexIn = vertices[cycleVertex].EndIn(); 
      for (inEdgeIndex = vertices[cycleVertex].FirstIn(); 
					 inEdgeIndex < endCycleVertexIn;
					 inEdgeIndex = vertices[cycleVertex].NextIn(inEdgeIndex)) {
				cycleInEdge = vertices[cycleVertex].in[inEdgeIndex];
	
				// sourceVertex --> cycleInEdge --> cycle
				sourceVertex = edges[cycleInEdge].src;

				if (collapsedVertex.VertexInCycle(sourceVertex)) {
					// Delete this edge, but don't bother routing it to the proxy
					//					deletedEdge = new Edge_t;
					//					*deletedEdge = edges[cycleInEdge];
					/*
						std::cout << "removing cycleedge " << cycleInEdge << " src: "
						<< edges[cycleInEdge].src << " "
						<< edges[cycleInEdge].dest << std::endl;
					*/
					// detach the edge from the graph, but remember where it came from
					ssize_t srcVertex, destVertex;
					srcVertex = edges[cycleInEdge].src;
					destVertex = edges[cycleInEdge].dest;
					RemoveEdge(vertices, edges, cycleInEdge, removedSrcIndex, removedDestIndex);
					// store information for the edge to reconstruct the graph later
					collapsedVertex.RecordRemovedEdge(cycleInEdge, 
																						srcVertex, removedSrcIndex, 
																						destVertex, removedDestIndex, edges[cycleInEdge].marked);
					continue;
				}
				// Otherwise, 'sourceVertex' is out of the cycle
				// look at all the edges that to out of source 
				// to see if they are in the cycle.  If any do, remove
				// the highest cost edges, and reroute the lowest
				// to go to the proxy vertex

				// Process all out edges into the cycle
				ssize_t minCycleInEdgeIndex = -1;
				ssize_t minCycleInEdgeScore = 0;
				ssize_t minCycleInEdge = -1;
				//UNUSED// int edgesThatMatch = 0;
				ssize_t minCycleDest = -1;

				// Now unlink all edges (sourceVertex, cycle), where cycle
				// is a vertex in the cycle.
				// Remove (sourceVertex, cycle) if it is not the optimal edge.
				// Replace (sourceVertex, cycle) with (sourceVertex, proxyVertex)
				// if the original edge was the optimal edge reaching the cycle.

				ssize_t sourceToCycleEdge;
				ssize_t deletedAnEdge = 0;
				for (outEdgeIndex = vertices[sourceVertex].FirstOut(); 
						 outEdgeIndex < vertices[sourceVertex].EndOut(); 
						 outEdgeIndex = vertices[sourceVertex].NextOut(outEdgeIndex)) {

					sourceToCycleEdge = vertices[sourceVertex].out[outEdgeIndex];

					// If this edge points in the cycle:
					//   We'll first look and see if it's the best edge from 'sourceVertex'
					//   into the cycle.
					//   Then unlink the edge from the graph (save unlinked edges???)
					if (collapsedVertex.EdgeEndsInCycle(edges[sourceToCycleEdge]) and
							!collapsedVertex.EdgeInCycle(sourceToCycleEdge)) {
						//					if (collapsedVertex.vertexToEdge.find(edges[sourceToCycleEdge].dest) !=
						//							collapsedVertex.vertexToEdge.end()) {
						// We have an edge:
						// sourceVertex -> sourceToCycleEdge -> (cycle)
						if ( minCycleInEdgeIndex == -1 or 
								 minCycleInEdgeScore > edges[sourceToCycleEdge].Cost()) {
							minCycleInEdgeIndex = outEdgeIndex;
							minCycleInEdgeScore = edges[sourceToCycleEdge].Cost();
							minCycleInEdge      = sourceToCycleEdge;
							minCycleDest        = edges[sourceToCycleEdge].dest;
						}

						// Unlink the previous destination vertex

						//						deletedEdge = new Edge_t;
						//						*deletedEdge = edges[sourceToCycleEdge];
						// detach the edge from the graph, but remember where it came from
						/*
							std::cout << "removing outside ->in edge " << sourceToCycleEdge << " "
							<< edges[sourceToCycleEdge].src << " "
							<< edges[sourceToCycleEdge].dest << std::endl;
						*/
						ssize_t srcVertex, destVertex;
						srcVertex  = edges[sourceToCycleEdge].src;
						destVertex = edges[sourceToCycleEdge].dest;
						
						RemoveEdge(vertices, edges, sourceToCycleEdge, removedSrcIndex, removedDestIndex);
						// store information for the edge to reconstruct the graph later
						collapsedVertex.RecordRemovedEdge(sourceToCycleEdge, 
																							srcVertex, removedSrcIndex, 
																							destVertex, removedDestIndex, edges[sourceToCycleEdge].marked);

						// this is just for a debug assertion later on
						deletedAnEdge = 1;
					}
				}
			
				// Make sure that at least one edge was found

				if (minCycleInEdgeIndex >= 0) {
					// Record where this edge came from

					// This is the lowest-scoring edge.  We need to save which 
					// vertex it is referencing so that when the collapsed vertex
					// is expanded, the min-cost edge may be marked.

					collapsedVertex.inEdges.push_back(std::pair<ssize_t, ssize_t>(minCycleInEdge, minCycleDest));
	  
					// Now add the edge to the proxy vertex.

					// Re-route the edge
					edges[minCycleInEdge].src  = sourceVertex;
					edges[minCycleInEdge].dest = proxyVertex;

					// Store the edge vertex adjacency lists
					newProxyInEdges.push_back(minCycleInEdge);
					//					vertices[proxyVertex].in.push_back(minCycleInEdge);
					vertices[sourceVertex].out[minCycleInEdgeIndex] = minCycleInEdge;
					/* print some debugging info, remove soon*/
					/*
						std::cout << "adding edge " << minCycleInEdge << " " << sourceVertex << " " << proxyVertex 
						<< " " << minCycleInEdgeIndex 
						<< " " << vertices[proxyVertex].in.size()-1 << std::endl;
					*/

					// Store the new minimum cost
					edges[minCycleInEdge].Cost(minCycleInEdgeScore);
				} // Done processing the vertex that has an edge e->(cycle)
				else {
					// If one edge was deleted, an optimal one must have been found.
					assert(deletedAnEdge == 0);
				}
	      
      }

    }
		// Done processing all vertices that have edges that go into the cycle
		// Santiy check.  All edges into the proxy vertex should have been removed.
		ssize_t proxyIn;
		ssize_t numIn = 0;
		//		std::cout << "old side: " << vertices[proxyVertex].in.size() << " ";
		for (proxyIn = 0; proxyIn < vertices[proxyVertex].in.size(); proxyIn++) {
			//			assert(vertices[proxyVertex].in[proxyIn] == -1 || 
			//						 collapsedVertex.EdgeInCycle(vertices[proxyVertex].in[proxyIn]));
			assert(vertices[proxyVertex].in[proxyIn] == -1);
			if (vertices[proxyVertex].in[proxyIn] != -1) ++numIn;
		}

		if (vertices[proxyVertex].in.size() < newProxyInEdges.size())
			vertices[proxyVertex].in.resize(newProxyInEdges.size());

		for (proxyIn = 0; proxyIn < newProxyInEdges.size(); proxyIn++) {
			vertices[proxyVertex].in[proxyIn] = newProxyInEdges[proxyIn];
		}
		//		std::cout << " new: " << vertices[proxyVertex].in.size() << std::endl;


		// ************************************************************
		// Step 4.
		// For all edges (c, v), where c is a cycle vertex, and v is outside
		// the cycle, unlink (c,v), and add (p,v), where p is the proxy
		// vertex to the collapsed cycle.
    // Now process all edges that go out of the cycle.

		
		std::vector<ssize_t> newProxyOutEdges;
    // Now, for each vertex that has an in-edge from the cycle, pick the top-scoring edge
    for (c = 0; c < cycle.size(); c++ ) {
      cycleVertex = cycle[c];
      // We don't need to re-route any edges that are already in the proxy

      // Update the cost and destination of all the edges in the graph to 
      // point to the proxy vertex
      ssize_t cycleOutEdge;
      ssize_t cycleOutEnd = vertices[cycleVertex].EndOut();
      for (outEdgeIndex = vertices[cycleVertex].FirstOut(); 
					 outEdgeIndex < cycleOutEnd; 
					 outEdgeIndex = vertices[cycleVertex].NextOut(outEdgeIndex)) {
				cycleOutEdge    = vertices[cycleVertex].out[outEdgeIndex];
				destVertex = edges[cycleOutEdge].dest;
				// This edge is part of the cycle, it's ok.
				if (collapsedVertex.vertexToEdge.find(destVertex) != 
						collapsedVertex.vertexToEdge.end())
					continue;

				// re-route edges from cycle to destVertex.
				// Process all edges going into this vertex, we'll keep the minimum cost one.
				ssize_t minCycleToDestInIndex, minCycleToDestInCost = -1;
				ssize_t minCycleToDestInEdge = -1;
				minCycleToDestInIndex = -1;
				for (inEdgeIndex = vertices[destVertex].FirstIn();
						 inEdgeIndex < vertices[destVertex].EndIn();
						 inEdgeIndex = vertices[destVertex].NextIn(inEdgeIndex)) {
					inEdge = vertices[destVertex].in[inEdgeIndex];
					// Check to see if this vertex is the destination from
					// the collapsed vertex.
					// If so, store minimal cost.
	  
					if (collapsedVertex.vertexToEdge.find(edges[inEdge].src) != 
							collapsedVertex.vertexToEdge.end() and 
							(minCycleToDestInIndex == -1 or 
							 minCycleToDestInCost > edges[inEdge].Cost())) {

						minCycleToDestInIndex = inEdgeIndex;
						minCycleToDestInCost  = edges[inEdge].Cost();
						minCycleToDestInEdge  = inEdge;
					}
				}

				if (minCycleToDestInEdge >= 0) {

					// Now unlink all edges going into the dest that 
					// come from the cycle, unless 
					// it happens to be the edge from the proxy, and is 
					// the minimal edge on the proxy
	  
	  
					// Record the optimal cycle source vertex for reconstruction 
					collapsedVertex.outEdges.push_back(std::pair<ssize_t, ssize_t>(edges[minCycleToDestInEdge].src, 
																																 minCycleToDestInEdge));
					ssize_t cycleToDestEdge;
					for (inEdgeIndex = vertices[destVertex].FirstIn();
							 inEdgeIndex < vertices[destVertex].EndIn();
							 inEdgeIndex = vertices[destVertex].NextIn(inEdgeIndex)) {

						cycleToDestEdge = vertices[destVertex].in[inEdgeIndex];
						// Check to see if this vertex is the destiation from
						// the collapsed vertex.
						// If so, if if this is not originating from the proxy vertex
						// unlink it.  We will need to add an edge to the proxy
						// If the best edge is coming from the proxy, just update
						// the score
						if (collapsedVertex.vertexToEdge.find(edges[cycleToDestEdge].src) !=
								collapsedVertex.vertexToEdge.end()) {
							// The edge comes from the cycle
							// Remove it if it is not from the proxy
							/*
								std::cout << "unlinking " << edges[cycleToDestEdge].src << " - " 
								<< cycleToDestEdge << " " << destVertex << std::endl;
							*/
							// And unlink the edge
							//							deletedEdge = new Edge_t;
							//							*deletedEdge = edges[cycleToDestEdge];

							ssize_t srcVertex, destVertex;
							srcVertex = edges[cycleToDestEdge].src;
							destVertex = edges[cycleToDestEdge].dest;

							RemoveEdge(vertices, edges, cycleToDestEdge, removedSrcIndex, removedDestIndex);
							collapsedVertex.RecordRemovedEdge(cycleToDestEdge, 
																								srcVertex, removedSrcIndex, 
																								destVertex, removedDestIndex, edges[cycleToDestEdge].marked);

						}
					}
					newProxyOutEdges.push_back(minCycleToDestInEdge);
					//					vertices[proxyVertex].out.push_back(minCycleToDestInEdge);
					// Link this edge into the graph
					edges[minCycleToDestInEdge].src = proxyVertex;

					// Althoug this shouldn't change, it was unlinked just above
					// so restore the edge links here.
					edges[minCycleToDestInEdge].dest = destVertex;

					vertices[destVertex].in[minCycleToDestInIndex] = minCycleToDestInEdge;
	  
					// Update the cost of this edge
					edges[minCycleToDestInEdge].Cost(minCycleToDestInCost);
				} // done processing a vertex that was reached by the cycle
      } // Done processing all out edges from the cycle
    }

		// sanity check.  all proxy out edges should ave been replaced.
		ssize_t proxyOut;
		for (proxyOut = 0; proxyOut < vertices[proxyVertex].out.size(); proxyOut++)
			assert(vertices[proxyVertex].out[proxyOut] == -1 || 
						 collapsedVertex.EdgeInCycle(vertices[proxyVertex].out[proxyOut]));
		
		if (vertices[proxyVertex].out.size() < newProxyOutEdges.size())
			vertices[proxyVertex].out.resize(newProxyOutEdges.size());

		for (proxyOut = 0; proxyOut < newProxyOutEdges.size(); proxyOut++)
			vertices[proxyVertex].out[proxyOut] = newProxyOutEdges[proxyOut];

    // Now, remove the edges that correspond to the cycle from the graph.
		// Getting rid of the non-proxy cycle edges.
		/*
    ssize_t removeEdge;
    for (c = 0; c < cycle.size(); c++ ) {
      cycleVertex = cycle[c];
      for (inEdgeIndex = vertices[cycleVertex].FirstIn();
					 inEdgeIndex < vertices[cycleVertex].EndIn();
					 inEdgeIndex = vertices[cycleVertex].NextIn(inEdgeIndex)) {
				inEdge = vertices[cycleVertex].in[inEdgeIndex];
				removeEdge = 0;
				if (cycleVertex == proxyVertex) {
					if (collapsedVertex.vertexToEdge.find(edges[inEdge].src) !=
							collapsedVertex.vertexToEdge.end())
						removeEdge = 1;
				}
				else {
					// This edge is not on the proxy.  Remove it
					// First, sanity check. 
					// All edges into or out of the cycle should have been processed, so this edge
					// should be back in the cycle.
					assert(collapsedVertex.VertexInCycle(edges[inEdge].src));
					removeEdge = 1;
				}
				if (removeEdge) {
					deletedEdge  =  new Edge_t;
					*deletedEdge = edges[inEdge];
					RemoveEdge(vertices, edges, inEdge, removedSrcIndex, removedDestIndex);
					collapsedVertex.RecordRemovedEdge(deletedEdge, inEdge, removedSrcIndex, removedDestIndex);
				}
      }
    }
		*/
  }
  
  template<typename Vertex_t, typename Edge_t>
		void ExpandCollapsedVertices(std::vector<Vertex_t> &vertices, std::vector<Edge_t> &edges,
																 std::vector<CollapsedVertex<Edge_t>*> &collapsedVertices) {
    
    ssize_t c;
    //UNUSED// int v, e;
    CollapsedVertex<Edge_t> *cv;
    ssize_t proxyVertex;
    ssize_t inEdge,inEdgeIndex, outEdge, outEdgeIndex;
    ssize_t removedInIndex, removedOutIndex;
    ssize_t optInEdge, optOutEdge;
    optOutEdge = -1;
    for (c = collapsedVertices.size()-1; c >= 0;  c-- ) {
      cv = collapsedVertices[c];
      // First get rid of all the edges that point to the proxy vertex
      proxyVertex = cv->proxyVertex;

      optInEdge = -1;
      for (inEdgeIndex = vertices[proxyVertex].FirstIn();
					 inEdgeIndex < vertices[proxyVertex].EndIn();
					 inEdgeIndex = vertices[proxyVertex].NextIn(inEdgeIndex)) {
				inEdge = vertices[proxyVertex].in[inEdgeIndex];
				if (edges[inEdge].marked == GraphEdge::Marked) {
					// Sanity check: we shouldn't have more than one edge
					// marked as optimal
					assert(optInEdge == -1);
					optInEdge = inEdge;
				}
				RemoveEdge(vertices, edges, inEdge, removedInIndex, removedOutIndex);
      }
      for (outEdgeIndex = vertices[proxyVertex].FirstOut();
					 outEdgeIndex < vertices[proxyVertex].EndOut();
					 outEdgeIndex = vertices[proxyVertex].NextOut(outEdgeIndex)) {
				outEdge = vertices[proxyVertex].out[outEdgeIndex];
				ssize_t r;
				for (r = 0; r < cv->edgeIndices.size(); r++ ) {
					if (cv->edgeIndices[r] == outEdge) {
						cv->optimal[r] = edges[outEdge].marked;
						//(cv->deletedEdges[r])->flagged = edges[outEdge].flagged;
					}
				}
						
				if (edges[outEdge].marked == GraphEdge::Marked) {
					optOutEdge = outEdge;
				}
				RemoveEdge(vertices, edges, outEdge, removedInIndex, removedOutIndex);
      }

      // Now add back all the edges that were deleted
      ssize_t r;
      ssize_t srcVertex, destVertex;
      ssize_t restoredEdge;
      for (r = 0; r < cv->edgeIndices.size(); r++ ) {
				restoredEdge = cv->edgeIndices[r];
				edges[restoredEdge].src = cv->srcVertex[r];
				edges[restoredEdge].dest = cv->destVertex[r];
				edges[restoredEdge].marked = cv->optimal[r]; //cv->deletedEdges[r]->marked; // cv->optimal[r];
				
				srcVertex = edges[restoredEdge].src;
				destVertex = edges[restoredEdge].dest;
				vertices[srcVertex].out[cv->srcIndices[r]]  = restoredEdge;
				vertices[destVertex].in[cv->destIndices[r]] = restoredEdge;
				/*
					std::cout << "restored edge " << restoredEdge << " " << srcVertex << " " << destVertex << " " 
					<< cv->srcIndices[r] << " " << cv->destIndices[r] << std::endl;
				*/
      }

      // Now record which edges are part of the MST
      //UNUSED// int i, o;
      ssize_t vertex;
      ssize_t in, out;
      ssize_t optInInEdgeIndex, optInInEdge;
      for (in = 0; in < cv->inEdges.size(); in++ ) {
				inEdge = cv->inEdges[in].first;
				vertex = cv->inEdges[in].second;
				/*
					std::cout << "setting in edge " << vertex << " " << inEdge << std::endl;
				*/
				inEdgeIndex = vertices[vertex].LookupInIndex(inEdge);
				// Record this as part of the MST:
				if (vertices[vertex].in[inEdgeIndex] == optInEdge) {

					// Now, construct the minimum arborescence of the tree by unlinking the
					// optimal in vertex of the cycle with it's predecessor
					for (optInInEdgeIndex = vertices[vertex].FirstIn();
							 optInInEdgeIndex < vertices[vertex].EndIn();
							 optInInEdgeIndex = vertices[vertex].NextIn(optInInEdgeIndex)) {
						optInInEdge = vertices[vertex].in[optInInEdgeIndex];
						edges[optInInEdge].marked = GraphEdge::NotMarked;
						edges[optInInEdge].flagged = GraphEdge::NotMarked;
					}

					// Now add the optimal in edge to the MST (maybe I should be using the mst flag, 
					// since that's what it's for...  I think I'm using marked for plotting purposes)
					edges[vertices[vertex].in[inEdgeIndex]].marked  = GraphEdge::Marked;
					edges[vertices[vertex].in[inEdgeIndex]].flagged = GraphEdge::Marked;
				}
      }

      // Now link the out-edges of the cycle
      for (out = 0; out < cv->outEdges.size(); out++ ) {
				vertex = cv->outEdges[out].first;
				outEdge = cv->outEdges[out].second;
				outEdgeIndex = vertices[vertex].LookupOutIndex(outEdge);
				// Record this as part of the MST:
				/*
					edges[vertices[vertex].out[outEdgeIndex]].marked  = GraphEdge::Marked;
					edges[vertices[vertex].out[outEdgeIndex]].flagged = GraphEdge::Marked;
				*/
      }
			delete collapsedVertices[c];
    }
  }
		    
  
  template<typename Vertex_t, typename Edge_t> 
		void CalcDirectedMST(std::vector<Vertex_t> &vertices, std::vector<Edge_t> &edges) {
    /*  
				This is an eimplementation of the Chu-Liu-Edwards Directed Minimum Spanning Tree.(DMST)
				The method is a greedy method that iteratively:
				1. Assigns to the DMST the minimum-cost in-edge into every vertex.  
				2. Collapses eacy cycle in the resulting DMST into a super-vertex
				Until there are no cylces in the DMST.  Then all super vertices are replaced 
				by the original collapsed vertices connected by the edges that formed the cycle - the 
				highest scoring edge in the cycle.

				The rules for collapsing a cycle are the following:
				Let C be the set of vertices in the cycle.
				Let u and w be vertices outside the cycle, and vc be a cycle vertex.
				u is any vertex where there is an edge (u, vc), and 
				w be any vertex where there is an edge (vc, w).

				Let cost(C)
				Redefine the cost of each edge (u, vc) to be cost(u,vc)  + cost(C).

				Remove all edges (u, vc1), ... (u, vc2) that are of higher cost than some edge (u, vc*). 
	
				... more documentation to come.
				A little bit about how I did the collapsing:
				There are two options for creating the super-vertices: adding a new vertex to the graph, or 
				transforming an exising vertex into a super-vertex.  Right now I'm doing the second so that
				no new vertices need to be added.  This might change since in the end it may use more memory
				than necessary.  Given the cycle-set of vertices, C, I store a data structure 'CollapsedVertex'
				that has the following properties:
				-proxyVertex: a vertex that is transformed into the super-vertex.  
	
    */
    std::vector<char> MSTIn;
    std::vector<ssize_t>  edgeCost;

    edgeCost.resize(edges.size());
    MSTIn.resize(edges.size());

    ssize_t dMSTFound = 0;
    std::vector<std::vector<ssize_t> > sccs;
    std::vector<CollapsedVertex<Edge_t>*> collapsedVertices;
    CollapsedVertex<Edge_t>* collapsedVertex;
    /*
			std::cout << "finding dmst " << std::endl;
    */
    ssize_t iteration = 0;
    ssize_t cycleIt = 0;
    while (dMSTFound == 0) {
      Unmark(edges);
      Unmark(vertices);
      ssize_t v,e;
      for (v= 0; v < vertices.size();v++) vertices[v].flagged = GraphVertex::NotMarked;
      for (e = 0; e < edges.size(); e++) edges[e].flagged = GraphEdge::NotMarked;
      sccs.clear();
      MarkLowestCostInEdges(vertices, edges);
			//			std::cout << "finding sccs" << std::endl;
      FindStronglyConnectedComponents(vertices, edges, sccs);
			//			std:::cout << std::endl;
      ssize_t cycle;
      ssize_t cycleVertex;
			for (cycle = 0 ; cycle < sccs.size(); cycle++ ) {
				for (cycleVertex = 0; cycleVertex < sccs[cycle].size(); cycleVertex++ ) {
					vertices[sccs[cycle][cycleVertex]].flagged = GraphVertex::Marked;
				}
      }
      /*
				std::stringstream itername;
				itername << "iteration." << iteration << ".dot";
				std::string gname = itername.str();
				GVZPrintBGraph(vertices, edges, gname);
				std::cout << "before contracting, printed graph to " << gname <<std::endl;
      */
      if (sccs.size() == 0) {
				/*
					std::cout << "NO CYCLES EXIST!!! " << std::endl;
				*/
				dMSTFound = 1;
      }
      else {
				/*
					std::cout << "the following cycles exist: " << std::endl;
				*/
				for (cycle = 0 ; cycle < sccs.size(); cycle++ ) {
					/*
						std::cout << cycle << " : ";
						for (cycleVertex = 0; cycleVertex < sccs[cycle].size(); cycleVertex++ ) {
						std::cout << " " << sccs[cycle][cycleVertex];
						}
						std::cout << std::endl;
					*/
					collapsedVertex = new CollapsedVertex<Edge_t>;
					ContractCycle(vertices, edges, sccs[cycle], *collapsedVertex );
					collapsedVertices.push_back(collapsedVertex);
					/*
						std::ofstream contractedGraph;
						std::stringstream name;
						name << "contracted." << cycleIt << ".dot";
						std::string graphName = name.str();
						std::cout << "contracted cycle of size " << sccs[cycle].size() 
						<< " and writing to " << graphName << std::endl;
						GVZPrintBGraph(vertices, edges, graphName);
					*/
					++cycleIt;
				}
      }
      ++iteration;
      
    } // Done contrating cycles
    ExpandCollapsedVertices(vertices, edges, collapsedVertices);
    ValidateDMST(vertices, edges);
  }

  template<typename V, typename E> 
		ssize_t ValidateDMST(std::vector<V> &vertices, std::vector<E> &edges) {
    ssize_t v;
    ssize_t nMarked;
    ssize_t inEdge, inEdgeIndex;
    for (v = 0; v < vertices.size(); v++ ) {
      nMarked = 0;
      for (inEdgeIndex = vertices[v].FirstIn();
					 inEdgeIndex < vertices[v].EndIn();
					 inEdgeIndex = vertices[v].NextIn(inEdgeIndex)) {
				inEdge = vertices[v].in[inEdgeIndex];
				/*	std::cout << v << " " << (int) edges[inEdge].marked << " " 
						<< (ssize_t) GraphEdge::Marked 
						<< " " << (ssize_t) edges[inEdge].flagged << std::endl;
				*/
				if (edges[inEdge].marked == GraphEdge::Marked)
					nMarked++;
      }
      if (nMarked != 1 and vertices[v].InDegree() > 0) {
				std::cout << "error: vertex: " << v 
									<< " does not have a marked in-edge " << std::endl;
				//				return 0;
      }
    }
    return 1;
  }
};

#endif
