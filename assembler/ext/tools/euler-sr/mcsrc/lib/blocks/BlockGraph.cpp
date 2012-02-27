/***************************************************************************
 * Title:          BlockGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "BlockGraph.h"
#include "Block.h"

#include <iostream>
#include <queue>

ssize_t BlockGraph::AddVertex(ssize_t start, ssize_t end) {
  BlockMap *map = new BlockMap;
  vertexMaps.push_back(map);
  Block *vertex = new Block(start, end, vertexMaps.size()-1);
  map->AddReference(vertex, vertexMaps.size()-1);
  return vertexMaps.size()-1;
}

void BlockGraph::GlueSequences(ssize_t refSeq, ssize_t startRef, ssize_t endRef,
			       ssize_t qrySeq, ssize_t startQry, ssize_t endQry) {
  // Glue the query strand into the reference. 
  assert(vertexMaps.size() > refSeq);
  assert(vertexMaps.size() > qrySeq);
   
  // Find out what blocks these go to.
  //UNUSED// Block *refSplitIn, *refSplitIsect, *refSplitOut;
  //UNUSED// Block *qrySplitIn, *qrySplitIsect, *qrySplitOut;
  BlockPath refPath, qryPath;
  
  FindPath(refSeq, startRef, endRef, refPath);

  Block *refStart, *refEnd;
  refStart = refPath.vertices[0];
  refEnd   = refPath.vertices[refPath.size()-1];
  Block *in, *out;
  in  = NULL;
  out = NULL;

  if (refStart != refEnd) {
    // The sequence to be glued maps to several edges
    // map to the start edge first.
    SplitVertex(refStart, refSeq, startRef-1, in, out); 
    SplitVertex(refEnd, refSeq, endRef, in, out);
  }
  else
    SplitVertex(refStart, refSeq, startRef-1, endRef);


  std::cout << "tree before splitting query vertex: " << std::endl;
  vertexMaps[qrySeq]->Print(std::cout);

  in  = NULL;
  out = NULL;

  FindPath(qrySeq, startQry, endQry, qryPath);
  Block *qryStart, *qryEnd;
  qryStart = qryPath.vertices[0];
  qryEnd   = qryPath.vertices[qryPath.size()-1];
  if (qryStart != qryEnd) {
    // The sequence to be glued maps to several edges
    // map to the start edge first.
    SplitVertex(qryStart, qrySeq, startQry-1, in, out); 
    SplitVertex(qryEnd, qrySeq, endQry, in, out); 
  }
  else 
    SplitVertex(qryStart, qrySeq, startQry-1, endQry);


  // Reset the paths to include new vertices
  refPath.vertices.clear();
  refPath.indices.clear();
  FindPath(refSeq, startRef, endRef, refPath);

  qryPath.vertices.clear();
  qryPath.indices.clear();
  FindPath(qrySeq, startQry, endQry, qryPath);
  std::cout << "reference path " << std::endl;
  std::cout << refPath << std::endl;
  std::cout << "query path " << std::endl;
  std::cout << qryPath << std::endl;
  
  //UNUSED// ssize_t m;
  // Glue the two paths together
  MergePaths(refPath, qryPath);
}

void BlockGraph::MergePaths(BlockPath &refPath, BlockPath &qryPath) {
  // Merge the query path into the reference path.
  
  /******
   * Sanity check.  Make sure the lengths of the paths are the same.
   */
  ssize_t startRef, endRef, startQry, endQry;
  ssize_t refStartIndex, refEndIndex;
  ssize_t qryStartIndex, qryEndIndex;
  assert(refPath.indices[refPath.size()-1] != -1 );

  // Find the coordinates of the ref seq
  //  refStartIndex = refPath.vertices[0]->GetSequenceIndex(refPath.indices[0]);
  refStartIndex = refPath.indices[0];
  startRef = refPath.vertices[0]->start[refStartIndex];
  
  //  refEndIndex = refPath.vertices.back()->GetSequenceIndex(refPath.indices.back());
  refEndIndex = refPath.indices.back();
  endRef   = refPath.vertices[refPath.size()-1]-> end[refEndIndex];

  // Double check the coordinates of the qry seq.
  //  qryStartIndex= qryPath.vertices[0]->GetSequenceIndex(qryPath.indices[0]);
  qryStartIndex= qryPath.indices[0];
  startQry = qryPath.vertices[0]->start[qryStartIndex];

  //  qryEndIndex = qryPath.vertices.back()->GetSequenceIndex(qryPath.indices.back());
  qryEndIndex = qryPath.indices.back();
  endQry   = qryPath.vertices.back()->end[qryEndIndex];

  std::cout << endRef << " " << startRef << " " << endQry << " " << startQry << std::endl;
  assert( (endRef - startRef) == (endQry - startQry));

  /**************
   * done with sanity check.
   */

  //UNUSED// Block *refVertex, *qryVertex, *endRefVertex, *endQryVertex;
  ssize_t   startRefPos, startQryPos, endRefPos, endQryPos;
  ssize_t   refSegmentLength, qrySegmentLength, refVertexLength, qryVertexLength;
  ssize_t   refIndex, qryIndex;
  refIndex = 0; 
  qryIndex = 0;
  PathIterator refIt, qryIt, nextIt;
  IndexIterator refIndIt, qryIndIt, nextIndIt;
  Block *in, *out;
  refIt = refPath.vertices.begin();
  qryIt = qryPath.vertices.begin();
  
  refIndIt = refPath.indices.begin();
  qryIndIt = qryPath.indices.begin();

  ssize_t step = 0;
  ssize_t refStep = 0;
  ssize_t qryStep = 0;
  while (refIt != refPath.vertices.end() || 
	 qryIt != qryPath.vertices.end()) {
    std::cout << "merging " << step << " th segment " << std::endl;
    step++;
    refIndex = (*refIndIt);
    qryIndex = (*qryIndIt);
    startRefPos = (*refIt)->start[refIndex];
    startQryPos = (*qryIt)->start[qryIndex];
    endRefPos = (*refIt)->end[refIndex];
    endQryPos = (*qryIt)->end[qryIndex];
    refSegmentLength = endRefPos - startRefPos;
    qrySegmentLength = endQryPos - startQryPos;
    std::cout << "trying to merge segments: "
	      << refSegmentLength 
	      << " " << qrySegmentLength << std::endl;
    if (refSegmentLength == qrySegmentLength) {
      if ( *refIt != *qryIt) {
	// If the two path segments are transitively equivalent 
	// they may have already been merged, so do not try and merge 
	// them again.
	std::cout << "merging two equal length segments " << std::endl;
	(*refIt)->Merge(*qryIt);
	// The merged edge is going to be removed, remove it's reference 
	// from the map.
	ssize_t s;
	ssize_t seqIdx;
	for (s = 0; s < (*qryIt)->sequence.size(); s++) {
	  seqIdx = (*qryIt)->sequence[s];
	  vertexMaps[seqIdx]->RemoveReference(*qryIt, (*qryIt)->sequence[s]);
	  vertexMaps[seqIdx]->AddReference(*refIt, (*qryIt)->sequence[s]);
	}

	// this node does not exist any more
	delete *qryIt;
	*qryIt = NULL;
      }
      else {
	std::cout << " not merging since they are the same! " << std::endl;
      }
      // advance
      ++refIt;
      ++qryIt;
      ++refIndIt;
      ++qryIndIt;
      ++refStep;
      ++qryStep;
    }
    else {
      // One of the vertices is shorter than the other. 
      // Split the longer vertex
      if (refSegmentLength < qrySegmentLength ) {
	refVertexLength = (*refIt)->end[0] - (*refIt)->start[0] + 1;
	std::cout << "ref segment ends before query " << std::endl
		  << " splitting the query vertex into lengths " 
		  << refVertexLength << " at: " 
		  << (*qryIt)->start[0] + refVertexLength - 1 << std::endl;
	SplitVertex(*qryIt,
		    (*qryIt)->sequence[*qryIndIt],
		    (*qryIt)->start[0] + refVertexLength - 1, in, out);

	// update the path, insert new vertex after the current vertex 
	// in a roundabout way by inserting it before the vertex after the 
	// current vertex.
	*qryIt = in;
	ssize_t qid;
	qid = *qryIndIt;
	++qryIt;
	++qryIndIt;
	qryIt = qryPath.vertices.insert(qryIt, out);
	qryIndIt = qryPath.indices.insert(qryIndIt, qid);
	--qryIt;
	--qryIndIt;
      }
      else {
	// split the reference vertex into two.
	qryVertexLength = (*qryIt)->end[0] - (*qryIt)->start[0] + 1;
	std::cout << "qry segment ends before ref " << std::endl
		  << " splitting the query vertex into lengths " 
		  << qryVertexLength << " at: " 
		  << (*refIt)->start[0] + qryVertexLength - 1 << std::endl;
	SplitVertex(*refIt, 
		    (*refIt)->sequence[*refIndIt], 
		    (*refIt)->start[0] + qryVertexLength - 1, in, out);

	// update the path
	*refIt = in;

	// Insert the new vertex after the current vertex
	++refIt;
	ssize_t rit;
	rit = *refIndIt;
	++refIndIt;
	refIt = refPath.vertices.insert(refIt, out);
	refIndIt = refPath.indices.insert(refIndIt, rit);
	--refIt;
	--refIndIt;
      }
    }
  }
}

Block* BlockGraph::FindNext(Block *block, ssize_t seq, ssize_t pos, Block *&next, ssize_t &s) {
  // use graph traversal to find the next block
  ssize_t n,c;
  //UNUSED// ssize_t foundNext;
  next = NULL;
  s = -1;
  // look through all the out nodes
  c = block->GetSequenceIndexPos(seq, pos);
  
  for (n = 0; n < block->out.size(); n++) {
    next = block->out[n];
    s = 0;
    while (s > -1) {
      s = next->GetSequenceIndexPos(seq, block->end[c]+1, s);
      if (s != -1) {
	assert(next->start[s] == block->end[c]+1);
	return next;
      }
    }
  }
  next = NULL;
  s = -1;
  return NULL;
}
 
void BlockGraph::FindPath(ssize_t seq, ssize_t start, ssize_t end, BlockPath &path) {
  Block *vertex, *nextVertex, *endVertex;
  ssize_t pos;
  ssize_t seqIndex, nextSeqIndex;
  std::cout << "finding path in: " << seq  << " from " << start << " to " << end << std::endl;
  assert(seq < vertexMaps.size());

  vertex = vertexMaps[seq]->FindEnclosingBlock(start);
  endVertex   = vertexMaps[seq]->FindEnclosingBlock(end);
  assert(vertex != NULL);
  assert(endVertex != NULL);
  
  pos = start;
  seqIndex = vertex->GetSequenceIndex(seq);

  while (vertex != endVertex && vertex != NULL) {
    path.vertices.push_back(vertex);
    path.indices.push_back(seqIndex);
    FindNext(vertex, seq, pos, nextVertex, nextSeqIndex);
    vertex = nextVertex;
    seqIndex = nextSeqIndex;
    if (vertex != NULL && vertex != endVertex)
      pos = vertex->start[seqIndex];
  }
  if (vertex != NULL) {
    path.vertices.push_back(vertex);
    path.indices.push_back(seqIndex);
  }
}

void BlockGraph::SplitVertex(Block *&vertex, ssize_t sequence, ssize_t pos, 
			     Block *&in, Block *&out) {
  // Create (up to) two vertices from previous one
  // pos is with respect to 'sequence'
  vertex->Split(pos, sequence, in, out);

  // No split happens if pos is not contained in the 
  // vertex boundaries.
  // This may be a good place for a sanity check
  if (in == NULL && out == NULL) 
    return;

  // Since a split happened, two new vertices must 
  // have been created.
  assert(in != NULL);
  assert(out != NULL);

  // Update all maps tha point to this sequence.

  // The reference is stored in map # sequence[0]
  //UNUSED+// ssize_t seqIndex;
  ssize_t s ;
  assert(vertex->sequence.size() == vertex->ref.size());
  for (s = 0; s < vertex->sequence.size(); s++) {
    std::cout << "removing vertex " << vertex << " from map " << vertex->sequence[s] << std::endl;
    BlockMap *map = vertexMaps[vertex->sequence[s]];
    // This vertex does not exist any more
    bool removeResult;
    removeResult = map->RemoveReference(vertex, vertex->sequence[s]);
    assert(removeResult == true);
    if (in != NULL) 
      map->AddReference(in, vertex->sequence[s]);
    
    if (out != NULL)
      map->AddReference(out, vertex->sequence[s]);
  }
  // vertex does not exist any more
  delete vertex;
  vertex = out;
}

void BlockGraph::SplitVertex(ssize_t sequence, ssize_t start, ssize_t end) {
  Block *block = NULL;
  block = vertexMaps[sequence]->FindEnclosingBlock(start);
  SplitVertex(block, sequence, start, end);
}

void BlockGraph::SplitVertex(Block *&vertex, ssize_t seq, ssize_t start, ssize_t end) {
  // Split a vertex into (possibly) three components, and update the 
  // map to reference these.

  // Split the block according to the start and end coordinates
  // The vertex represents an alignment from pos start to pos end, 
  // so the first vertex (first split should go from the beginning 
  // of the vertex to start-1, and the second split goes from
  // start all the way until end.

  Block *in, *out; // discard in and out for now.  not hungry
  ssize_t seqIndex;
  in = NULL;
  out = NULL;
  seqIndex= vertex->GetSequenceIndex(seq);
  assert (seqIndex >= 0 && seqIndex < vertex->size());
  if (vertex->start[seqIndex] <= start)
    SplitVertex(vertex, seq, start, in, out);
  if (vertex->end[seqIndex] > end)
    SplitVertex(vertex, seq, end, in, out);
}

void BlockGraph::PrintPath(ssize_t map, ssize_t startPos, ssize_t endPos, std::ostream &out) {
  Block *vertex;
  Block *end;
  vertex = vertexMaps[map]->FindEnclosingBlock(startPos);
  end = vertexMaps[map]->FindEnclosingBlock(endPos);

  ssize_t done;
  if (vertex == NULL)
    done = 1;
  else 
    done = 0;

  ssize_t i;
  //  std::cout << "got end: " << end << std::endl;
  while (! done && vertex != NULL ) {
    for (i = 0; i < vertex->size(); i++) 
      out << vertex->sequence[i] << "\t" << vertex->start[i] << "\t" << vertex->end[i] << "\t";
    out << std::endl;
    if (vertex == end)
      done = 1;
    else {
      vertex = vertexMaps[map]->FindEnclosingBlock(vertex->end[0] + 1);
    }
  }
}

void BlockGraph::ResetVisited() {
 BlockMapIterator mapIt;
  Block *vertex;
  ssize_t m;

  // Reset the visited status of every vertex
  for (m = 0; m < vertexMaps.size(); m++) {
    vertexMaps[m]->BeginIterator(mapIt);
    while (mapIt.NotNil()) {
      vertex = mapIt.Data();
      /*      std::cout << "checking vertex: (" << vertex->sequence[0] << ") " << vertex->start[0] 
		<< " " << vertex->end[0] << std::endl;
      */
      assert(vertex != NULL);
      vertex->marked = 0;
      vertex->number = 0;
      mapIt.Next();
    }
  }
}

void BlockGraph::PrintGraph(std::ostream &out) {
  
  // Print a gvz compatible output of the graph.

  ssize_t m;
  ssize_t n;
  ssize_t o, i;
  BlockMapIterator mapIt;
  Block *vertex;
  // Mark each node as not visited.
  ResetVisited();
  n = 1;
  out << "digraph G {" << std::endl;
  out << "\tsize=\"8,8\";" << std::endl;
  out << "label=\"test graph\"" << std::endl;
  for (m = 0; m < vertexMaps.size(); m++) {
    vertexMaps[m]->BeginIterator(mapIt);
    while (mapIt.NotNil()) {
      vertex = mapIt.Data();
      if (! vertex->marked ) {
	if ( vertex->number == 0 ) {
	  vertex->number = n;
	  ++n;
	}
	for (o = 0; o < vertex->out.size(); o++) {
	  if ( vertex->out[o]->number == 0) {
	    vertex->out[o]->number = n;
	    ++n;
	  }
	  out << "\t" << vertex->number << " -> " //<<  vertex->out[o] << " " 
	      <<  vertex->out[o]->number << "[label=\"";
	  for (i = 0; i < vertex->size(); i++) 
	    out << " " << vertex->sequence[i] << " " << vertex->start[i] << " " << vertex->end[i] << ", ";
	  out << " -> " << vertex->out[o]->sequence[0] << " " 
	      <<  vertex->out[o]->start[0] << " " << vertex->out[o]->end[0] ;
	  out << "\"];" << std::endl;
	}
	vertex->marked = 1;
      }
      mapIt.Next();
    }
  }

  out << "}" << std::endl;
}

