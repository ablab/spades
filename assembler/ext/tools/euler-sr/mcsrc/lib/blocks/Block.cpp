/***************************************************************************
 * Title:          Block.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "Block.h"

std::ostream &operator<<(std::ostream &out, const Block &b)
{ 
  ssize_t i;
  for (i =0; i < b.size(); i++) 
    out << "[ " << b.sequence[i] << ": " << b.start[i] << " " << b.end[i] << " " <<  "] ";
  return out;
}

void Block::AddCoords(ssize_t qryStrand, ssize_t qryStart, ssize_t qryEnd) {
  std::cout << "adding coords " << qryStrand << " into block: " << this << std::endl;
  
  sequence.push_back(qryStrand);
  start.push_back(qryStart);
  end.push_back(qryEnd);
  ref[qryStrand] = NULL; // start this block off w/o a reference
}

void Block::Merge(Block *vertex) {
  ssize_t i, o;
  // add the paths in vertex to this block 
  for (i = 0; i < vertex->sequence.size(); i++) {
    std::cout << "merging sequence " << vertex->sequence[i] << " into block: " << this << std::endl;
    sequence.push_back(vertex->sequence[i]);
    start.push_back(vertex->start[i]);
    end.push_back(vertex->end[i]);
    ref[vertex->sequence[i]] = NULL; // wait for this to be assigned.
  }

  // fix the edges of the vertices that point to vertex
  Block *prev, *next;
  for (i = 0; i <vertex->in.size(); i++) {
    // fix the in edges
    prev = vertex->in[i];
    in.push_back(prev);

    // fix the out edges of the vertices that are connected to 
    // vertex.
    for (o = 0; o < prev->out.size(); o++) {
      if ( prev->out[o] == vertex ) {
	// prev points to the vertex that needs to be moved.
	prev->out[o] = this;
      }
    }
  } // done fixing in edges

  for (o = 0; o < vertex->out.size(); o++) {
    next = vertex->out[o];
    // fix this nodes out edges
    out.push_back(next);
    for (i = 0; i < next->in.size(); i++) { 
      if (next->in[i] == vertex) {
	next->in[i] = this;
      }
    }
  }
}

ssize_t  Block::GetRefPos(ssize_t pos, ssize_t seqIndex) {
  // ******* sanity ********
  assert(seqIndex >= 0 && seqIndex < size());
  // either the position is within the range specified by the start sequence, or is
  // immediately before
  assert((pos >= (start[seqIndex]-1) && pos <= end[seqIndex]));
   // ******* sanity ********

  // map the position in non-reference sequence to the 
  // position in the reference seq.
  pos = start[0] + (pos - start[seqIndex]);
  
  // ******* sanity ********
  //assert(pos <= end[0]);
  assert(pos >= (start[0]-1));
  // ********end sanity ****

  return pos;
}
void Block::Split(ssize_t pos, ssize_t seq, Block *&inVertex, Block *&outVertex) {
  // Take a vertex and split it in two at position out within the
  // sequence, seq.  
  // 
  ssize_t seqIndex;
  seqIndex  = GetSequenceIndex(seq);
  // if the pos is not in the reference sequence, map it to the
  // reference (0th sequence) pos for this block
  if (seqIndex != 0) {
    pos = GetRefPos(pos, seqIndex);
  }

  // Make sure that the pos is within the coordinates of the referene
  // sequence of this block.
  ssize_t inLength, outLength;
  ssize_t i,o;
  inLength = 0;
  outLength = 0;
  
  inVertex = NULL;
  outVertex = NULL;
  
  // Make sure that we're not splitting a single nucleotide
  assert(start[0] < end[0] or pos < start[0]);
  
  if (pos < start[0]) {
    return;
  }
  if (pos >= end[0]) {
    return;
  }

  inVertex = new Block(start[0], pos, sequence[0]);
  inLength = pos - start[0] + 1;
  // add aligned coords
  for (i = 1; i < start.size(); i++)
    inVertex->AddCoords(sequence[i], start[i], start[i] + inLength - 1);
  // link in
  for (i = 0; i < in.size(); i++) 
    inVertex->in.push_back(in[i]);

  outVertex = new Block(pos+1, end[0], sequence[0]);
  // second vertex starts after pos and goes until the end.
  // add aligned coords
  outLength = end[0] - (pos+1) + 1; 
  for (i = 1; i < start.size(); i++)
    outVertex->AddCoords(sequence[i],  // map to this sequence
			 start[i] + inLength,  // start here
			 start[i] + inLength + outLength-1 // end here.
			 );
  // link out
  for (i = 0; i < out.size(); i++) 
    outVertex->out.push_back(out[i]);
  
  outVertex->in.push_back(inVertex);
  inVertex->out.push_back(outVertex);

  // Fix external edges
  for (i = 0; i < in.size(); i++) {
    for (o = 0; o < in[i]->out.size(); o++)
      if (in[i]->out[o] == this) {
	in[i]->out[o] = inVertex; 
      }
  }
  for (o = 0; o < out.size(); o++) {
    for (i = 0; i < out[o]->in.size(); i++)
      if ( out[o]->in[i] == this ) {
	out[o]->in[i] = outVertex;
      }
  }
}

void Block::Split(ssize_t startPos, ssize_t endPos,  // The start and end indices in the sequence.
		                       // They are exact indices into the sequence, rather
		                       // than offsets into the vertex
		                       // The three possible resulting blocks
		  ssize_t seq,
		  Block *&in, Block *&intersect, Block *&out) {
  //  vertex  
  //  s0                     e0
  //  (----------------------)
  //  is turned into:
  //  s0    s-1,s...e,e+1... e0                          
  //  (-----)->(----)->(------)
  //    in      isect    out
  //  Alghough if s-1 < s0, no in is created
  //   and if e+1 > e0, no out is created.
  
  // located these coordinates within the vertex, split this vertex into three.
  in = out = NULL;  // Might not create these vertexs if the bounaries are an 
  // exact overlap.
  ssize_t i;
  ssize_t seqIndex;
  seqIndex= GetSequenceIndex(seq);
  if (seqIndex != 0) {
    startPos = GetRefPos(startPos, seqIndex);
    endPos   = GetRefPos(endPos, seqIndex);
  }

  if (startPos-1 > start[0]) {
    in = new Block(start[0], startPos-1, sequence[0]);

    // Update the coordinates of the other paths that map through here
    for (i = 1; i < sequence.size(); i++)
      in->AddCoords(sequence[i], 
		    start[i], 
		    start[i] + (startPos-1 -start[0]+1));
  }
  // Middle vertex is always added
  intersect = new Block(startPos, endPos, sequence[0]);

  for (i = 1; i < sequence.size(); i++) 
    intersect->AddCoords(sequence[i], 
			 start[i] + (startPos - start[0]+1),
			 start[i] + (endPos - start[0]+1));
  // Create edges
  if (in != NULL) {
    in->out.push_back(intersect);
    intersect->in.push_back(in);
  }
    
  // Create out vertex
  if (endPos+1 < end[0]) {
    out = new Block(endPos+1, end[0], sequence[0]);
    // Copy overlapping coordinates
    for (i = 1; i < sequence.size(); i++) {
      out->AddCoords(sequence[i], 
		     start[i] + ((endPos+1) - start[0] + 1),
		     end[i]);
    }
  }
  // Create edges
  if (out != NULL) {
    out->in.push_back(intersect);
    intersect->out.push_back(out);
  }
}

ssize_t Block::GetSequenceIndexPos(ssize_t seq, ssize_t pos, ssize_t startIndex) {
  ssize_t i;
  for (i = startIndex; i < start.size(); i++) {
    if (seq == sequence[i] &&
	pos >= start[i] &&
	pos <= end[i])
      return i;
  }
  return -1;
}

ssize_t Block::GetSequenceIndex(ssize_t seq, ssize_t startIndex) {
  ssize_t i;
  for (i = startIndex; i < sequence.size(); i++) {
    if (sequence[i] == seq)
      return i;
  }
  return -1;
}
