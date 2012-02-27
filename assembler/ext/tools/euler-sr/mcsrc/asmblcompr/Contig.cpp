/***************************************************************************
 * Title:          Contig.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "Contig.h"

#include "MUMAlignBlock.h"

ssize_t Contig::SumNumAlignBlocks() {
  ssize_t c;
  ssize_t size = 0;
  for (c = 0; c < alignedClusters.size(); c++ ) {
    size += alignedClusters[c]->totalSize();
  }
  return size;
}


void Contig::FindIntervals() {
  ssize_t c, a, b, p;
  MUMCluster *cluster;
  MUMAlignBlock *block;
  //UNUSED// int startIntervalFound, endIntervalFound;
  ssize_t start, end;
	//UNUSED// int between;
  ssize_t startInside, endInside;
  ssize_t startBetween, endBetween;

  for (c = 0; c < alignedClusters.size(); c++ ) {
    cluster = alignedClusters[c];
    for (b = 0; b < cluster->size(); b++) {
      block = cluster->blocks[b];
      for (a = 0; a < block->size(); a++) {
				startInside = -1;
				endInside = -1;
				startBetween = -1;
				endBetween   = -1;
				start = block->refPos[a];
				end   = start + block->length[a];

				// Create the first interval
				if (intervals.size() == 0) {
					intervals.startPositions.push_back(start);
					intervals.endPositions.push_back(end);
					continue;
				}
	
				// Try to find the start of this interval
				for (p = 0; p < intervals.size(); p++) {
					if (intervals.startPositions[p] > start) {
						startBetween = p;
						break;
					}
					else if (intervals.startPositions[p] <= start and
									 intervals.endPositions[p] >= start) {
						startInside = p;
						break;
					}
				}
	
				if (startInside == -1 or p == intervals.size()) {
					// Did not find an interval for the starting position, create one
					if (p == intervals.size())
						// Didn't find an interval
						startInside = intervals.size();
					else
						startInside = startBetween;
					intervals.startPositions.insert(intervals.startPositions.begin() + startInside, start);
					intervals.endPositions.insert(intervals.endPositions.begin() + startInside, start);
					// Inserted a new interval, don't need to check it to see if the end interval is here
				}
	
				// Now try to find a location for the end interval
				ssize_t p2;
				p2 = intervals.size() - 1;
				for (; p2 >= p ; p2--) {
					if (intervals.endPositions[p2] < end) {
						endBetween = p2 ; // end after this
						break;
					}
					else if (intervals.endPositions[p2] <= end and
									 intervals.startPositions[p2] >= end) {
						endInside = p2;
						break;
					}
				}
	
				if (startBetween == -1 and 
						endBetween == -1 and
						endInside == startInside) {
					// This interval is eclipsed by a preexisting one.
					continue;
				}

				if (endBetween != -1 or p2 < p) {
					// endInside references before where you want to insert
					if (p2 < p) 
						endInside = p + 1 ;
					else
						endInside = endBetween + 1 ;
					if (endInside == intervals.startPositions.size()) {
						intervals.startPositions.insert(intervals.startPositions.end(), end);
						intervals.endPositions.insert(intervals.endPositions.end(), end);
					}
					else {
						intervals.startPositions.insert(intervals.startPositions.begin() + endInside, end);
						intervals.endPositions.insert(intervals.endPositions.begin() + endInside, end);
					}
					// The new index is one greater than where it was between.
				}

				// Now merge everything between the two intervals.
				assert(startInside != endInside);

				intervals.endPositions[startInside] = intervals.endPositions[endInside];
	    
				intervals.startPositions.erase(intervals.startPositions.begin() + startInside + 1,
																			 intervals.startPositions.begin() + endInside + 1);
				intervals.endPositions.erase(intervals.endPositions.begin() + startInside + 1,
																		 intervals.endPositions.begin() + endInside + 1);
      }
    }
  }
}

ssize_t Contig::CalculateCoverage() {
  // Simple statistic for determining the coverage 
  covered = 0;
  ssize_t i;
  for (i = 0; i < intervals.size(); i++) 
    covered += intervals.endPositions[i] - intervals.startPositions[i] + 1;

  return covered;
}
