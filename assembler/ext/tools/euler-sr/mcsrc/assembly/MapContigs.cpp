/***************************************************************************
 * Title:          MapContigs.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SortedTupleList.h"
#include "DNASequence.h"
#include "utils.h"
#include "SeqReader.h"
#include "ReadPos.h"
#include "MapContigs.h"


void MapContig(DNASequence &contig,
							 SimpleSequenceList &sequences,
							 std::vector<CountedReadPos> &refPositions,
							 int tupleSize, ssize_t minBlockLength,
							 BlockList &blocks) {
	ssize_t p;
	ssize_t index;
	ssize_t p0;
	p0 = 0;
	ssize_t curBlockStart = -1;
	ssize_t curBlockEnd   = 0;
	ssize_t queryStart    = 0;
	// Store the starting and

	ssize_t curRefPos, curRefStrand, prevRefPos, prevRefStrand;
	curRefPos = -1;
	curRefStrand = -1;
	Block curBlock;
	//UNUSED// ssize_t forIndex, revIndex;
	DNASequence contigRC;

	for (p = 0; p < contig.length - tupleSize + 1; p++ ) {
		index = LocateTuple(sequences, refPositions, tupleSize, 
													 (char*) &(contig.seq[p]));

		if (index != -1) {
			prevRefPos    = curRefPos;
			prevRefStrand = curRefStrand;
			curRefPos     = refPositions[index].pos;
			curRefStrand  = refPositions[index].read;
			//			std::cout << prevRefPos << " " << curRefPos << " " << std::endl;
		}
		else {
			continue;
		}

		// Found adjacent entries
		// the second condition will almost never matter
		if (curRefPos != -1 and  // found a match
				curRefPos == prevRefPos + 1 and  // a contiguous match
				curRefStrand == prevRefStrand) {
			// Not in strip, store the start of it.
			if (curBlockStart == -1) {
				curBlockStart = prevRefPos;
				queryStart    = p-1;
			}
			// update the end of the strip
			curBlockEnd = curRefPos;
		}
		else {
			// Two non-adjacent entries.

			if (curBlockStart != -1) {
				// Previous pos was in a strip, store that.
				curBlock.refPos = curBlockStart;
				curBlock.qryPos = queryStart;
				curBlock.length = curBlockEnd - curBlockStart + 1;
				curBlock.strand = curRefStrand;
				/*
				std::cout << "bloc " <<  blocks.size()  << " " << curBlockStart 
									<< " " << queryStart << " " << curBlock.length << std::endl;
				*/
				blocks.push_back(curBlock);
			}
			curBlockStart = -1;
		}	
	}
	// process the last block.
	if (curBlockStart != -1) {
		curBlock.refPos = curBlockStart;
		curBlock.qryPos = queryStart;
		curBlock.length = curBlockEnd - curBlockStart + 1;
		curBlock.strand = curRefStrand;
		/*
		std::cout << "block " <<  blocks.size()  << " " << curBlockStart 
							<< " " << queryStart << " " << curBlock.length << std::endl;
		*/
		blocks.push_back(curBlock);
	}

	// Discard blocks below a certain length;
	ssize_t changeMade;
	do {
		
		// Assume all blocks are large enough to not be removed.
		changeMade = 0;
		BlockList::iterator blockIt, nextBlockIt;
		blockIt = blocks.begin();
		while (blockIt != blocks.end()) {
			if ((*blockIt).length < minBlockLength) {
				blockIt = blocks.erase(blockIt);
				changeMade = 1;
			}
			else {
				/*
				std::cout << blocks.size() << " " << (*blockIt).refPos << " " << (*blockIt).qryPos
									<< " " << (*blockIt).length << std::endl;
				*/
				++blockIt;
			}					
		}
		//		std::cout << "merging adjacent blocks."  << std::endl;
		// Merge newly adjacent blocks
		blockIt = blocks.begin();
		nextBlockIt = blockIt; nextBlockIt++;

		// brief check to see if all blocks have been
		// removed (all are short)
		if (blockIt == blocks.end())
			break;

		while (nextBlockIt != blocks.end()) {
			// Found two newly adjacent blocks.
			/*
			std::cout << (*nextBlockIt).refPos << " " << (*blockIt).refPos << " "
								<< (*nextBlockIt).qryPos << " " << (*blockIt).qryPos << " " 
								<< (*nextBlockIt).refPos - (*blockIt).refPos << " " 
								<< (*nextBlockIt).qryPos - (*blockIt).qryPos << std::endl;
			*/
			if ((*nextBlockIt).refPos - (*blockIt).refPos ==
					(*nextBlockIt).qryPos - (*blockIt).qryPos) {

				(*blockIt).length = (*nextBlockIt).length + 
					(*nextBlockIt).refPos - (*blockIt).refPos;

				blockIt++;
				// erase the next element
				blockIt = blocks.erase(blockIt);
				if (blockIt != blocks.end()) {
					nextBlockIt = blockIt;
					nextBlockIt++;
				}
				else {
					nextBlockIt = blocks.end();
				}
				changeMade = 1;
			}
			else {
				++blockIt;
				++nextBlockIt;
			}
		}
		/*
		std::cout << "remaining blocks: " << std::endl;
		for(blockIt = blocks.begin(); blockIt != blocks.end(); ++blockIt) {
			std::cout << (*blockIt).refPos << " " << (*blockIt).qryPos << " " 
								<< (*blockIt).length << std::endl;

		}
		std::cout << "end" << std::endl;
		*/
	}
	while (changeMade);
	BlockList::iterator blockIt;
	for(blockIt = blocks.begin(); blockIt != blocks.end(); ++blockIt) {
		(*blockIt).length += tupleSize;
	}
}

