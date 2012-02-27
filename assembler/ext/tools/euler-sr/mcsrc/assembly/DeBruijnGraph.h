/***************************************************************************
 * Title:          DeBruijnGraph.h 
 * Author:         Mark Chaisson, Boyko Kakaradov, Glenn Tesler
 * Created:        2007
 * Last modified:  12/14/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef DE_BRUIJN_GRAPH_H_
#define DE_BRUIJN_GRAPH_H_

#include "compatibility.h"
#include "DNASequence.h"
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "ReadPos.h"
#include "SeqUtils.h"
#include "utils.h"
#include "graph/GraphAlgo.h"
#include "graph/MSTAlgo.h"
#include "SortedTupleList.h"


template<typename ReadPosType>
std::ostream &operator<<(std::ostream &out, const std::vector<ReadPosType> &rp) {
  out << rp.size() << std::endl;
  ssize_t i;
  for (i = 0; i < rp.size(); i++ ) {
    out << rp[i] << std::endl;
  }
	return out;
}


template<typename ReadPosType>
std::istream &operator>>(std::istream &in, std::vector<ReadPosType> &rp) {
  ssize_t size;
  in >> size;
  ssize_t i;
  rp.resize(size);
  for (i = 0; i < size; i++ ) {
    in >> rp[i];
    if (! in ) {
      std::cout << "Error reading position " << i << " of " << size << std::endl;
      exit(1);
    }
  }
  return in;
}

template<typename V>
std::istream &ReadVertexList(std::istream &in, std::vector<V> &vertices) {
  ssize_t size;
  std::string word;
  in >> word >> size;
  ssize_t index;
  vertices.resize(size);
  ssize_t i;
  for (i = 0; i < size; i++ ) {
    in >> index;
    ReadVertex(in, vertices[i]);
		vertices[i].index = index;
  }
  return in;
}

template<typename V>
std::ostream &WriteVertexList(std::ostream &out, std::vector<V> &vertices) {
  ssize_t i;
  ssize_t numVertices = 0;
  for (i = 0; i < vertices.size(); i++ ) {
    if  (vertices[i].OutDegree() > 0 or
				 vertices[i].InDegree() > 0) {
      numVertices++;
    }
  }
 
  out << "Number_of_vertices " << numVertices << std::endl;

  for (i = 0; i < vertices.size(); i++ ) {
    if (vertices[i].OutDegree() > 0 or
				vertices[i].InDegree() > 0) {
      out << i << " ";
      WriteVertex(out, vertices[i]);
      out << std::endl;
    }
  }
  return out;
}

template<typename V>
std::ostream &WriteVertex(std::ostream &outs, V &vertex) {
	vertex.Write(outs);
  return outs;
}

template<typename V>
void CondenseEdgeLists(std::vector<V> &vertexList) {
	ssize_t v;
	for (v = 0; v < vertexList.size(); v++) {
		vertexList[v].CondenseEdgeLists();
	}
}

template<typename V>
std::istream &ReadVertex(std::istream &ins, V &vertex) {
  ins >> vertex.index;
	std::string token;
	ins >> token;
	ssize_t outIndex, inIndex;
	ssize_t i;
	if (token == "out") {
		ssize_t nOut;
		ins >> nOut;
		for (i = 0; i < nOut; i++ ){
			ins >> outIndex;
			vertex.out.push_back(outIndex);
		}
		ins >> token;
		if (token != "in") {
			std::cout << "ERROR reading vertex, expected 'in'."<< std::endl;
			exit(1);
		}
		ssize_t nIn;
		ins >> nIn;
		for(i = 0; i < nIn; i++ ) {
			ins >> inIndex;
			vertex.in.push_back(inIndex);
		}
	}
	else {
		vertex.out[0] = atosz(token.c_str());
		ins >> vertex.out[1] >> vertex.out[2] >> vertex.out[3]
				>> vertex.in[0] >> vertex.in[1] >> vertex.in[2] >> vertex.in[3];
	}
  return ins;
}

template<typename E>
std::ostream &WriteEdgeList(std::ostream &out, std::vector<E> &edges) {
  out << "Number_of_edges " << edges.size() << std::endl;
  ssize_t i;
  for (i = 0; i < edges.size(); i++) { 
    if (! edges[i].IsNullified() )
      out << i << " " << edges[i] << std::endl;
  }
  return out;
}
template<typename E>
std::istream &ReadEdgeList(std::istream &in, std::vector<E> &edges ) {
  ssize_t size, index;
  std::string word;
  in >> word >> size;
  edges.resize(size);
  ssize_t i;
  for (i = 0; i < size; i++) {
		edges[i].Init();
    in >> index >> edges[i];
		edges[i].index = index;
		//    edges[i].index = index;
  }
  return in;
}


template<typename V>
ssize_t StoreZeroDegreeVertexIndices(std::vector<V> &vertices, std::vector<ssize_t> &indices ) {
  ssize_t v;
  for (v = 0; v < vertices.size(); v++ ) {
    if (vertices[v].InDegree() == 0 and
				vertices[v].OutDegree() == 0) 
      indices.push_back(v);
  }
  return indices.size();
}

template<typename V, typename E>
	ssize_t RemoveZeroDegreeVertices(std::vector<V> &vertices, std::vector<E> &edges) {
  ssize_t nZero = 0;
  std::vector<ssize_t> removedIndices;
  nZero = StoreZeroDegreeVertexIndices(vertices, removedIndices);

  // If there is nothing to fix, end now
  if (nZero == 0) {
    return 0;
  }
  DeleteElementsFromList(vertices, removedIndices);
  UpdateRemovedVertexIndices(edges, removedIndices);
  
  return nZero;
}

template<typename V>
ssize_t DeleteElementsFromList(std::vector<V> &list, std::vector<ssize_t> &removedIndices) {
  // Now copy the removed vertices from the list
  if (removedIndices.size() == 0) 
    return 0;
  ssize_t curRemoved = 0;
  ssize_t v;
  for (v = 0; v < list.size(); v++ ) {
    //    stopfn(v);
    if (curRemoved < removedIndices.size() and
				v == removedIndices[curRemoved]) {
      curRemoved++;
    }
    else {
      //      std::cout << "prev: " << v << " new: " << v - curRemoved << std::endl;
      list[v - curRemoved] = list[v];
    }
  }
  list.resize(list.size() - removedIndices.size());
  return removedIndices.size();
}

template<typename E>
ssize_t UpdateRemovedEdgeBalancedIndices(std::vector<E> &edges,
																		 std::vector<ssize_t> &removedEdges) {
  // Once edges have been removed from the list, 
  // the indices of the balanced edges need to be updated
  // to the new ordering of edges.
  // In the case of updating vertex in/out edge indices and edge src/dest
  // vertex indices, we don't expect to find an index in the removedIndices 
  // list because only unremoved indices remain.  However, since we're simply
  // re-packing the edges here (indices 1 2 X ALPHABET_SIZE X X 5 6 X 7 -> 
  // 1 2 3 4 5 6), the balanced edges have nothing to do with whether or not
  // they are from a removed edge.

  ssize_t e;
  ssize_t balancedEdge;
  if (removedEdges.size() == 0) 
    return 0;

  std::vector<ssize_t>::iterator edgeIt;
  ssize_t numLowerRemovedEdges;
  for (e = 0; e < edges.size(); e++ ) {
    balancedEdge = edges[e].balancedEdge;
    assert(balancedEdge != -1);
    if (std::binary_search(removedEdges.begin(),
													 removedEdges.end(), balancedEdge) != 0) {
      std::cout << "edge: " << balancedEdge 
								<< " should not be in the edge list " << std::endl;
      assert(0);
    }
    if (edges[e].balancedEdge >= removedEdges[0]) {
      if (balancedEdge > removedEdges[removedEdges.size()-1]) {
				numLowerRemovedEdges = removedEdges.size();
      }
      else {
				edgeIt = std::lower_bound(removedEdges.begin(),
																	removedEdges.end(),
																	balancedEdge);
				numLowerRemovedEdges = edgeIt - removedEdges.begin();
      }
      edges[e].balancedEdge -= numLowerRemovedEdges;
    }
  }
  return 0;
}

template<typename V>
void UpdateRemovedEdgeIndices(std::vector<V> &vertices,
															std::vector<ssize_t> &removedIndices) {
  ssize_t v;
  ssize_t i, o;
  ssize_t edgeIndex;
  ssize_t numberRemovedEdgesBelow;
  std::vector<ssize_t>::iterator nextEdgeIndex;
  for (v = 0; v < vertices.size(); v++ ) {
    //squeeze edge indices back down after removing them
    for (i = 0; i < vertices[v].EndIn(); i++ ) {
      edgeIndex = vertices[v].in[i];
      if (std::binary_search(removedIndices.begin(), 
														 removedIndices.end(), edgeIndex) != 0) {
				std::cout << "edge index: " << edgeIndex 
									<< " should not be part of vertex "
									<< v << std::endl;
				assert(0);
      }
      if (edgeIndex >= removedIndices[0]) {
				nextEdgeIndex = std::lower_bound(removedIndices.begin(),
																				 removedIndices.end(),
																				 edgeIndex);
				if(nextEdgeIndex == removedIndices.end())
					numberRemovedEdgesBelow = removedIndices.size();
				else
					numberRemovedEdgesBelow = nextEdgeIndex - removedIndices.begin();
				vertices[v].in[i] -= numberRemovedEdgesBelow;
      }
    }
    // repeated code, except using 'out'
    for (o = 0; o < vertices[v].EndOut(); o++ ) {
      edgeIndex = vertices[v].out[o];
      if (std::binary_search(removedIndices.begin(), 
														 removedIndices.end(), edgeIndex) != 0) {
				std::cout << "out edge: " << edgeIndex 
									<< " should not be a part of " 
									<< v << std::endl;
				assert(0);
      }
      if (edgeIndex >= removedIndices[0]) {
				nextEdgeIndex = std::lower_bound(removedIndices.begin(),
																				 removedIndices.end(),
																				 edgeIndex);
				if(nextEdgeIndex == removedIndices.end())
					numberRemovedEdgesBelow = removedIndices.size();
				else
					numberRemovedEdgesBelow = nextEdgeIndex - removedIndices.begin();
				vertices[v].out[o] -= numberRemovedEdgesBelow;
      }
    }
  }
}

template<typename E>
void UpdateRemovedVertexIndices(std::vector<E> &edges, std::vector<ssize_t> &removedIndices) {
  // Now offset the vertex indices stored in the edge list
  ssize_t e;
  ssize_t vertex;
  ssize_t numberRemovedVerticesBelow;
  std::vector<ssize_t>::iterator nextIndex;
  
  for (e = 0; e < edges.size(); e++) {
    vertex = edges[e].src;

    // No edge should reference a removed vertex.
    if (std::binary_search(removedIndices.begin(), removedIndices.end(), vertex) != 0) {
      std::cout << "vertex " << vertex << " should not be src of  " << e << std::endl;
      assert(0);
    }

    // If the vertex has vertices with indices lower than it that 
    // have been removed, we need to fix the index so that only
    // present indices remain.
    if (vertex >= removedIndices[0]) {
      nextIndex = lower_bound(removedIndices.begin(), removedIndices.end(), vertex);
      if (nextIndex == removedIndices.end())
				numberRemovedVerticesBelow = removedIndices.size();
      else 
				numberRemovedVerticesBelow = nextIndex - removedIndices.begin();
      edges[e].src -= numberRemovedVerticesBelow;
    }

    // Now the same thing with the destination vertex
    vertex = edges[e].dest;
    if (binary_search(removedIndices.begin(), removedIndices.end(), vertex) != 0) {
      std::cout << "vertex " << vertex  << " should not be dest of " << e << std::endl;
      assert(0);
    }
    if (vertex >= removedIndices[0]) {
      nextIndex = lower_bound(removedIndices.begin(), removedIndices.end(), vertex);
      if (nextIndex == removedIndices.end())
				numberRemovedVerticesBelow = removedIndices.size();
      else 
				numberRemovedVerticesBelow = nextIndex - removedIndices.begin();
      edges[e].dest -= numberRemovedVerticesBelow;
    }
  }
}



template<typename T, typename E>
	void PrintBGraph(std::vector<T> &vertices, std::vector<E> &edges, int vertexSize, std::ostream &out) {
  WriteVertexList(out, vertices);
  WriteEdgeList(out, edges);
	out << vertexSize << std::endl;
}

template<typename T, typename E>
	void PrintBGraph(std::vector<T> &vertices, std::vector<E> &edges, int vertexSize, std::string &outName,
									 std::ostream &report = std::cout
									 ) {
  std::ofstream out;
  openck(outName, out, std::ios::out, report);
  PrintBGraph(vertices, edges, vertexSize, out);
  out.close();
}


template<typename V, typename E>
	ssize_t ReadBGraph(std::string &inName, std::vector<V> &vertices, std::vector<E> &edges,
								 std::ostream &report = std::cout
								 ) {
  std::ifstream in;
  openck(inName, in, std::ios::in, report);
  ReadVertexList(in, vertices);
  ReadEdgeList(in, edges);
	int vertexSize;
	if (!(in >> vertexSize)) {
		vertexSize = 1;
	}
  in.close();
	return vertexSize;
}

template<typename V, typename E>
	void CheckVertices(std::vector<V> &vertices, std::vector<E> &edges) {
  ssize_t v;
  ssize_t nEdges = edges.size();
  for (v = 0; v < vertices.size(); v++ ) {
    assert(vertices[v].in[0] < nEdges);
    assert(vertices[v].in[1] < nEdges);
    assert(vertices[v].in[2] < nEdges);
    assert(vertices[v].in[3] < nEdges);
    assert(vertices[v].out[0] < nEdges);
    assert(vertices[v].out[1] < nEdges);
    assert(vertices[v].out[2] < nEdges);
    assert(vertices[v].out[3] < nEdges);
  }
}

template<typename V, typename E>
	ssize_t CheckEdges(std::vector<V> &vertices, std::vector<E> &edges) {
  ssize_t e;
	ssize_t src, dest, srcEdgeIndex, destEdgeIndex;
  for (e = 0; e < edges.size(); e++) {
		if (edges[e].src >= 0 and edges[e].dest >= 0) {
			src = edges[e].src;
			srcEdgeIndex = vertices[src].LookupOutIndex(e);
			assert(srcEdgeIndex >= 0);
			dest = edges[e].dest;
			destEdgeIndex = vertices[dest].LookupInIndex(e);
			assert(destEdgeIndex >= 0);
		}
  }
  return 1;
}


template<typename E>
ssize_t CheckBalancedEdgeList(std::vector<E> &edges, std::vector<ssize_t> &indices) {
  ssize_t i;
  for (i = 0; i < indices.size(); i++ ) {
    if (!std::binary_search(indices.begin(), indices.end(), edges[indices[i]].balancedEdge)) {
      std::cout << "edge " << indices[i] << " is in the list but not " << edges[indices[i]].balancedEdge << std::endl;

      return 0;
    }
  }
  return 1;
}


template<typename E>
void CheckBalance(std::vector<E> &edges) {
  ssize_t e;
  ssize_t balancedEdge;
  for (e = 0; e < edges.size(); e++ ) {
    balancedEdge = edges[e].balancedEdge;
    assert(balancedEdge >= 0);
    assert (edges[balancedEdge].balancedEdge == e);
  }
}

template<typename Vertex_t, typename Edge_t>
	class PrintComponentFunctor {
 public:
  ssize_t size;
  ssize_t degree;
  ssize_t component;
  ssize_t edgeLength;
  std::vector<Vertex_t> *vertices;
  std::vector<Edge_t> *edges;

  PrintComponentFunctor() {
    component = 0;
  }
  void Reset() {
    size = 0; 
    degree = 0;
    edgeLength = 0;
  }
  void operator()(ssize_t vertexIndex) {

    ssize_t i;
    for (i = 0; i < (*vertices)[vertexIndex].EndOut(); i++ ) {
      if ((*vertices)[vertexIndex].out[i] >= 0) {
				edgeLength += (*edges)[(*vertices)[vertexIndex].out[i]].length;
				/*
					std::cout << " e: " << (*vertices)[vertexIndex].out[i] << " (";
					std::cout << (*edges)[(*vertices)[vertexIndex].out[i]].length << " b: ";
					std::cout << (*edges)[(*vertices)[vertexIndex].out[i]].balancedEdge << ") ";
				*/
      }
    }
    size ++;
    degree += (*vertices)[vertexIndex].OutDegree();
  }
};

template<typename Vertex_t, typename Edge_t>
	void PrintComponents(std::vector<Vertex_t> &vertices,
											 std::vector<Edge_t> &edges) {
  PrintComponentFunctor<Vertex_t, Edge_t> printComponent;

  printComponent.vertices = &vertices;
  printComponent.edges = &edges;
  ssize_t v;
  printComponent.component = 0;
  for (v = 0; v < vertices.size(); v++ ) {
    if (vertices[v].marked == GraphVertex::NotMarked) {
      printComponent.Reset();
      TraverseDFS(vertices,edges,v,printComponent);
      std::cout << printComponent.component << "\t" 
								<< printComponent.size << "\t" 
								<< printComponent.degree << "\t"
								<< printComponent.edgeLength << std::endl;
      printComponent.component++;
    }
  }
}

template<typename Edge_t>
void ReadSequences(std::string edgeFileName, std::vector<Edge_t> &edges,
									 std::ostream &report = std::cout
									 ) {
  SimpleSequenceList sequences;
  ReadSimpleSequences(edgeFileName, sequences, report);
  assert(sequences.size() == edges.size());
  ssize_t e;
  for (e = 0; e < sequences.size(); e++ ) {
    edges[e].seq.seq = sequences[e].seq;
    edges[e].seq.length = sequences[e].length;
    edges[e].length  = sequences[e].length;
  }
}

template<typename Vertex_t, typename Edge_t>
	void PrintEdges(std::vector<Vertex_t> &vertices, std::vector<Edge_t> &edges, std::string &edgeFileName,
									std::ostream &report = std::cout
									) {
  //UNUSED+// ssize_t  v;
  ssize_t e;
  DNASequence tmpSeq;
  std::string blank;
  std::ofstream seqOut;
  openck(edgeFileName, seqOut, std::ios::out, report);
  //UNUSED// ssize_t edgeIndex;
	for (e = 0; e < edges.size(); e++) {
		
		//  for (v = 0; v < vertices.size(); v++ ) {
		//    for (e = 0; e < vertices[v].EndOut(); e++ ) {
		//      if (vertices[v].out[e] >= 0) {
		//				edgeIndex = vertices[v].out[e];
		//		tmpSeq.seq = edges[edgeIndex].seq.seq;
		tmpSeq.seq = edges[e].seq.seq;
		tmpSeq.length = edges[e].seq.length;
		tmpSeq.ClearName();
		*tmpSeq.titlestream << edges[e].src << " -> " << edges[e].dest 
												<< " (" << e << ")";
		tmpSeq.PrintSeq(seqOut);
		seqOut << std::endl;
		//	}
		//    }
  }  
}

template<typename Vertex_t, typename Edge_t>
	void PrintGraph(std::vector<Vertex_t> &vertices, std::vector<Edge_t> &edges, std::string graphOutName,
									std::ostream &report = std::cout
									) {
  std::ofstream graphOut;
  openck(graphOutName, graphOut, std::ios::out, report);
  ssize_t v, e;
  ssize_t numVertices = 0;
  for (v = 0; v < vertices.size(); v++ ){
    if (vertices[v].OutDegree() != 0 or
				vertices[v].InDegree() != 0 )
      numVertices++;
  }
  graphOut << "Number_of_Vertex " << numVertices << std::endl;

  ssize_t balEdge;
  ssize_t curVertex = 0;
  for (v = 0; v < vertices.size(); v++ ) {
    if (vertices[v].InDegree() != 0 or 
				vertices[v].OutDegree() != 0) {
      graphOut << "Vertex " << curVertex << " "  << vertices[v].OutDegree() << " "
							 << vertices[v].InDegree() << std::endl;
      graphOut << "Last_edge";
      for (e = 0; e < vertices[v].EndIn(); e++) { 
				if (vertices[v].in[e] >= 0) {
					graphOut << " " << vertices[v].in[e];
					balEdge = edges[vertices[v].in[e]].balancedEdge;
					graphOut << " " << balEdge;
				}
      }
      graphOut << std::endl;
      graphOut << "Next_edge";
      for (e = 0; e < vertices[v].EndOut(); e++) { 
				if (vertices[v].out[e] >= 0) {
					graphOut << " " << vertices[v].out[e];
					balEdge = edges[vertices[v].out[e]].balancedEdge;
					graphOut << " " << balEdge;
				}
      }
      graphOut << std::endl;
      curVertex++;
    }
  }
  graphOut.close();
}

/* BEGIN INSERT BOYKO */
template<typename Vertex_t, typename Edge_t>
	class GVZPrintComponentFunctor {
 public:
	ssize_t size;
	ssize_t degree;
	ssize_t component;
	ssize_t edgeLength;
	std::vector<Vertex_t> *vertices;
	std::vector<Edge_t> *edges;
	std::ofstream gvzOut;
	std::ostream *report;

	GVZPrintComponentFunctor(ssize_t component,
													 std::vector<Vertex_t> *vertices,
													 std::vector<Edge_t> *edges,
													 std::string &bgraphFileName,
													 std::ostream *report
													 ) {
		this->component = component;
		this->vertices = vertices;
		this->edges = edges;
		this->report = report;

		std::string gvzOutFileName = bgraphFileName;
		gvzOutFileName += "_files/";
		if (component > 0) {
			gvzOutFileName += NumToStr(component);
		} else if (component < 0) {
			gvzOutFileName += NumToStr(-component);
			gvzOutFileName += "d";
		} else {
			// TODO: special case for component = 0
			gvzOutFileName += NumToStr(component);
		}
		gvzOutFileName += ".dot";
		openck(gvzOutFileName.c_str(), gvzOut, std::ios::out, *report);
		gvzOut << "digraph G {"<< std::endl
					 << "\tsize=\"8,8\";"<< std::endl
					 << "\toverlap=\"false\";" << std::endl;
	}

	~GVZPrintComponentFunctor() {
		gvzOut << "}" << std::endl;
		gvzOut.close();
	}

	void Reset() {
		size = 0;
		degree = 0;
		edgeLength = 0;
	}

	void operator()(ssize_t vertexIndex) {
		gvzOut << "\t" << vertexIndex << " [";
		// TODO: to get the sequence I must use BVertex instead of BranchingVertex
		// gvzOut << "label=\"" << (*vertices)[vertexIndex].seq.seq << "\" ";
		if ((*vertices)[vertexIndex].flagged == GraphVertex::Marked)
			gvzOut << " , color=\"gray\",style=filled ";
		gvzOut << "];" << std::endl;

		for (ssize_t e = 0; e < (*vertices)[vertexIndex].EndOut(); e++) {
			ssize_t edgeIndex = (*vertices)[vertexIndex].out[e];
			if (edgeIndex >= 0) {
				edgeLength += (*edges)[edgeIndex].length;
				gvzOut << "\t" << vertexIndex << " -> " << (*edges)[edgeIndex].dest
				       << " [label=\"" << (*edges)[edgeIndex].index
							 << " l=" << (*edges)[edgeIndex].length
							 << " m=" << (*edges)[edgeIndex].multiplicity << "\"";
				if ((*edges)[edgeIndex].flagged == GraphEdge::Marked)
					gvzOut << " ,color=red ";
				gvzOut << "] weight=\"0.5\" fontsize=\"10\" len=\"3\"";
				// BEGIN BOYKO: add tooltips like so:
				// gvzOut << " [tooltip=\"" << edges[edgeIndex].seq.seq << "\"]" ;
				// END BOYKO
				gvzOut << ";" << std::endl;
				/*
					gvzOut << " e: " << (*vertices)[vertexIndex].out[e] << " (";
					gvzOut << (*edges)[(*vertices)[vertexIndex].out[e]].length << " b: ";
					gvzOut << (*edges)[(*vertices)[vertexIndex].out[e]].balancedEdge << ") ";
				*/
			}
		}
		size++;
		degree += (*vertices)[vertexIndex].OutDegree();
	}
};

template<typename V, typename E>
	void GVZPrintBGraph(std::vector<V> &vertices, std::vector<E> &edges, std::string &graphOutName,
											std::ostream &report = std::cout
											) {
  ssize_t v;
  std::string command = "mkdir -p " + graphOutName + "_files/";
  std::system(command.c_str()); // TODO: error check & log; add routine to utils.h
	Unmark(vertices); // Clear all marks.  Loop will mark each vertex as used.


	// special file for singletons
	GVZPrintComponentFunctor<V, E>
		gvzPrintComponentSingle(0, &vertices, &edges, graphOutName, &report);

  ssize_t component = 1;
  for (v = 0; v < vertices.size(); v++) {
		//		if (vertices[v].flagged == GraphVertex::Marked) continue; // next vertex
		if (vertices[v].marked == GraphVertex::Marked) continue; // next vertex

		// Collect all isolated vertices and edges into one file.
		ssize_t d_in, d_out;
		bool is_isolated = false;
		d_in = vertices[v].InDegree();
		d_out = vertices[v].OutDegree();

		if (d_in + d_out <= 1) {
			// cases:
			// 1. isolated vertex with no edges
			// 2. two distinct vertices connected by an edge
			// 3. isolated vertex with one or more loops --- TODO

			if (d_in == 0 && d_out == 0) {
				is_isolated = true;
			} else if (d_in + d_out == 1) {
				ssize_t v2 = -1;
				if (d_out == 1) {
					ssize_t outEdgeIndex = vertices[v].FirstOut();
					assert(outEdgeIndex < vertices[v].EndOut());
					ssize_t edge = vertices[v].out[outEdgeIndex];
					v2 = edges[edge].dest;
				} else { // d_out == 0 && d_in == 1
					ssize_t inEdgeIndex = vertices[v].FirstIn();
					assert(inEdgeIndex < vertices[v].EndIn());
					ssize_t edge = vertices[v].in[inEdgeIndex];
					v2 = edges[edge].src;
				}
				if (v2 != -1
						&& vertices[v2].InDegree() + vertices[v2].OutDegree() <= 1) {
					is_isolated = true;
				}
			}

			if (is_isolated) {
				// isolated vertex or edge.  Print it (and dual if known) to
				// singleton file.

				TraverseDFS(vertices,edges,v,gvzPrintComponentSingle);

				ssize_t balVertex = LookupBalancedVertex(vertices, edges, v);
				if (balVertex != -1
						&& vertices[balVertex].marked != GraphVertex::Marked) {

					TraverseDFS(vertices,edges,balVertex,gvzPrintComponentSingle);
				}

				continue;
			}
		}

		// Not isolated.  Print component and dual component to separate files.
		// First print component.

		GVZPrintComponentFunctor<V, E>
			gvzPrintComponent(component, &vertices, &edges, graphOutName, &report);

		TraverseDFS(vertices,edges,v,gvzPrintComponent);

		// Then print dual component, if different

		ssize_t balVertex = LookupBalancedVertex(vertices, edges, v);
		if (balVertex != -1
				&& vertices[balVertex].marked != GraphVertex::Marked) {


			GVZPrintComponentFunctor<V, E>
				gvzPrintComponentDual(-component, &vertices, &edges, graphOutName,
															&report);

			TraverseDFS(vertices,edges,balVertex,gvzPrintComponentDual);
		}

		component++;
	}
}
/* END INSERT BOYKO */

#if 0 // DISABLE OLD GVZPrintBGraph - REPLACED BY BOYKO'S
template<typename V, typename E>
	void GVZPrintBGraph(std::vector<V> &vertices, std::vector<E> &edges, std::string &graphOutName,
											std::ostream &report = std::cout
											) {
  std::ofstream graphOut;
  openck(graphOutName, graphOut, std::ios::out, report);
  
  ssize_t v, e;
  graphOut << "digraph G {"<<std::endl
					 << "\tsize=\"8,8\";"<<std::endl
					 << "\toverlap=\"false\";" << std::endl;
	

  ssize_t edgeIndex;
  for (v = 0; v < vertices.size(); v++ ) {
    for (e = 0; e < vertices[v].EndOut(); e++ ) {
      if (vertices[v].out[e] != -1) {
				edgeIndex = vertices[v].out[e];
				graphOut << "\t" << v << " -> " << edges[edgeIndex].dest;
				graphOut << " [label=\"" << edgeIndex <<" {" << edges[edgeIndex].index << "}";
				graphOut << " (" << edges[edgeIndex].seq.length << ")\""; 
				if (edges[edgeIndex].flagged == GraphEdge::Marked)
					graphOut << " ,color=red ";
	    
				graphOut <<  "] weight=\"0.5\" fontsize=\"10\" len=\"3\"" << std::endl;

				// BEGIN BOYKO: add tooltips like so:
				// graphOut << " [tooltip=\"" << edges[edgeIndex].seq.seq << "\"]" ;
				// END BOYKO

				graphOut << ";" << std::endl;

      }
    }
    graphOut << "\t" << v << " [label=\"";
    graphOut << " " << v << " : " << vertices[v].index << "\" ";
    if (vertices[v].flagged == GraphVertex::Marked)
      graphOut << " " << ", color=\"gray\",style=filled  ";
    graphOut << "]" << std::endl;

		// BEGIN BOYKO: add tooltips like so:
		// graphOut << " [tooltip=\"" << edges[edgeIndex].seq.seq << "\"]" ;
		// END BOYKO

    graphOut << ";" << std::endl;

  }
  graphOut << "}" << std::endl;
}
#endif // DISABLE OLD GVZPrintBGraph - REPLACED BY BOYKO'S


template<typename Vertex_t, typename Edge_t> 
	void CalcBalancedMST(std::vector<Vertex_t> &vertices, std::vector<Edge_t> &edges) {
  Unmark(vertices);
  Unmark(edges);
    
  std::vector<Edge_t*> edgePtrs;
  edgePtrs.resize(edges.size());
  ssize_t e;
  for (e = 0; e < edges.size(); e++ )
    edgePtrs[e] = &(edges[e]);

  GraphAlgo::CompareEdgePtrs<Edge_t> compareEdgePtrs;
  std::sort(edgePtrs.begin(), edgePtrs.end(), compareEdgePtrs);

  std::vector<DisjointSet::DJVertex<Vertex_t > > djForest;
  DisjointSet::CreateForest(vertices, djForest);
  
  ssize_t srcVertex, destVertex;
  Edge_t *edge;
  ssize_t balEdge, balSrc, balDest;
  ssize_t balancedEdge;
  for (e = 0; e < edgePtrs.size(); e++ ) {
    // add e
    edge = edgePtrs[e];
    srcVertex = edge->src;
    destVertex = edge->dest;
    if (DisjointSet::Find(&(djForest[srcVertex])) != 
				DisjointSet::Find(&(djForest[destVertex]))) {
      DisjointSet::Union(&(djForest[srcVertex]), 
												 &(djForest[destVertex]));
      std::cout << "adding " << edge->index << " to mst " << std::endl;
      edge->mst = GraphAlgo::MSTIn;

      balancedEdge = edge->balancedEdge;
      balSrc = edges[balancedEdge].src;
      balDest = edges[balancedEdge].dest;
      if (balancedEdge != edge->index) {
				// the balanced edge should not have been added yet to the graph
				if (DisjointSet::Find(&(djForest[balSrc])) != 
						DisjointSet::Find(&(djForest[balDest]))) {
					DisjointSet::Union(&(djForest[balSrc]), 
														 &(djForest[balDest]));
					std::cout << "adding balanced edge: " << balancedEdge << " to mst " << std::endl;
					edges[balancedEdge].mst = GraphAlgo::MSTIn;
				}
				else {
					if (edges[balancedEdge].balancedEdge != edge->index) {
						std::cout << "Is there something strange about this? added a balanced edge," << std::endl;
						std::cout << " but not the normal edge " << std::endl;
					}
				}
      }
    }
  }
}

template<typename V, typename E>
	void RemoveEdge(std::vector<V> &vertices, std::vector<E> &edges, ssize_t edge,
									ssize_t &srcIndex, ssize_t &destIndex) {
  ssize_t src, dest;
  src = edges[edge].src;
  dest = edges[edge].dest;
  srcIndex = vertices[src].LookupOutIndex(edge);
  assert(srcIndex >=0);
	assert(dest >= 0);
	assert(dest <= vertices.size());
  destIndex = vertices[dest].LookupInIndex(edge);
  assert(destIndex >=0);
  vertices[src].out[srcIndex] = -1;
  vertices[dest].in[destIndex] = -1;
  edges[edge].src = -1;
  edges[edge].dest = -1;
}
#endif
