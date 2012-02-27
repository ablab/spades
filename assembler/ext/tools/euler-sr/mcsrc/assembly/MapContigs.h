/***************************************************************************
 * Title:          MapContigs.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MAP_CONTIGS_H_
#define MAP_CONTIGS_H_
#include "SortedTupleList.h"
#include "ReadPos.h"
#include "DNASequence.h"
#include <list>

// class that is an ungapped alignment.
class Block {
 public:
	ssize_t refPos;
	ssize_t qryPos;
	ssize_t length;
	ssize_t strand;
};

typedef std::list<Block> BlockList;

void MapContig(DNASequence &contig,
							 SimpleSequenceList &sequences,
							 std::vector<CountedReadPos> &refPositions,
							 int tupleSize, ssize_t minBlockLength,
							 BlockList &blocks);

#endif
