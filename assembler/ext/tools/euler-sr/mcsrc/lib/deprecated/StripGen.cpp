/***************************************************************************
 * Title:          StripGen.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "StripGen.h"

#include <stdlib.h>
#include <cmath>

long totalCalls, totalInserts;


void EnumerateStrips(T_Strips &strips) {
  T_StripQQryPos queue;
  BuildStripPqueue(strips, queue);
  AssignEnumerations(queue);
}



ssize_t DiscardAndMerge(T_Strips &strips, ssize_t distThreshold) {
  ssize_t startCount, numDiscardedStrips;
  startCount = strips.size();
  numDiscardedStrips = DiscardSmallStrips(strips, distThreshold);
  EnumerateStrips(strips);
  Merge(strips);
  return numDiscardedStrips;
}

void Merge(T_Strips &strips) {
  T_Strips::iterator cur, prev;
  // Merge adjacent strips
  prev = cur = strips.begin();
  ++cur;
  for (; cur != strips.end();) {
    /*    std::cout << "merging strips: " << strips.size() << " "
	  << prev->startQryEnum << " " << cur->startQryEnum;*/
    if (MergeStrips(prev, cur, 0) == 1) {
      cur = strips.erase(cur);
      //      std::cout << " (merged) " << std::endl;
    }
    else {
      prev = cur;
      ++cur;
      //      std::cout << std::endl;
    }
  }
}

void PrintEnumerations(T_Strips &strips, std::ostream &out) {
  T_StripQRef pqueue;
  long curEnum = 1;
  out << strips.size() << std::endl;
  T_Strips::iterator it, end = strips.end();
  for (it = strips.begin(); it != end; ++it)
    pqueue.push(it);

  T_Strips::iterator stripIt;
  while (pqueue.size() > 0) {
    stripIt = pqueue.top();
    pqueue.pop();
    out << stripIt->startQryEnum << "\t" << stripIt->startRefPos 
	<< "\t" << stripIt->endRefPos << "\t" << stripIt->startQryPos << "\t" 
	<< stripIt->endQryPos << "\t" << stripIt->size << "\t" 
	<< stripIt->numMerged << std::endl;
  }
}


ssize_t MergeStrips(T_Strips::iterator prev, T_Strips::iterator cur, ssize_t forceMerge= 0) {
  // Merge current strip into previous strip.  The current strip may be 
  // removed after this since it is now redundant. 

  // Condition 1: the removed strip is a reversed strip in the middle
  // of two other strips 
  ssize_t doMerge = 0;
  if (prev->sign == cur->sign) {
    if (prev->sign == -1 && (prev->endQryEnum - 1 ==
			     cur->startQryEnum))
      doMerge = 1;
    else if (prev->endQryEnum + 1 == cur->startQryEnum) 
      doMerge = 1;
  }
  if (doMerge || forceMerge) {
    prev->size += (cur->size - 
		   std::max(prev->endRefPos + prev->hashLength - cur->startRefPos, (long)0));
    prev->endRefEnum = cur->endRefEnum;
    prev->endQryEnum = cur->endQryEnum;
    prev->endRefPos = cur->endRefPos;
    prev->endQryPos = cur->endQryPos;
    prev->numMerged += 1 + cur->numMerged;
    return 1;
  }
  else
    return 0;
}

ssize_t DiscardSmallStrips(T_Strips& strips, ssize_t distThresh) {
  //  std::cout << "discarding small stips " << std::endl;
  T_Strips::iterator cur, prev, end;
  Strip removedStrip;
  prev = strips.end();
  ssize_t numSkipped = 0;
  //  prev = strips.end();
  ssize_t size;
  //  for (cur = strips.begin(); cur != strips.end(); ) {
  ssize_t startRefPos;
  cur = strips.begin();
  startRefPos = cur->startRefPos;
  startRefPos;
  while (cur != strips.end()) {
    if (cur->size < distThresh) {
      // Remove the current stip since it is too small.
      numSkipped++;
      // Will need to keep track of the removed strip to make sure
      // that the signs of the adjacent strips are consistent.
      cur = strips.erase(cur);  // advance cur to the next strip
      if(prev != strips.end() && cur != strips.end() && MergeStrips(prev, cur)) {
	cur = strips.erase(cur);
      }
    }
    else {
      prev = cur;
      ++cur;
    }
    strips.size();
  }
  return numSkipped;
}

void ReadStrips(std::string inName, 
		ssize_t *&enumerations, ssize_t *&startRefLocations, ssize_t *&endRefLocations,
		ssize_t *&startQryLocation, ssize_t *&endQryLocation, ssize_t *&stripSize, 
		ssize_t *&stripsMerged, ssize_t &length) {
  std::ifstream in;
  in.open(inName.c_str());
  if (!in.good()) {
    std::cout << "Could not open " << inName << std::endl;
    exit(0);
  }
  ReadStrips(in, enumerations, startRefLocations, endRefLocations,
		startQryLocation, endQryLocation, stripSize, 
	     stripsMerged, length);
  in.close();
  in.clear();
}


void ReadStrips(std::ifstream &inFile, 
		ssize_t *&enumerations, ssize_t *&startRefLocations, ssize_t *&endRefLocations,
		ssize_t *&startQryLocation, ssize_t *&endQryLocation, ssize_t *&stripSize, 
		ssize_t *&stripsMerged, ssize_t &length) {
  inFile >> length;
  enumerations = new ssize_t[length];
  startRefLocations = new ssize_t[length];
  endRefLocations = new ssize_t[length];
  startQryLocation = new ssize_t[length];
  endQryLocation = new ssize_t[length];
  stripSize = new ssize_t[length];
  stripsMerged = new ssize_t[length];
  ssize_t i;  
  for (i = 0; i < length; i++) {
    inFile >> enumerations[i] >> startRefLocations[i] 
	   >> endRefLocations[i] >> startQryLocation[i] 
	   >> endQryLocation[i] >> stripSize[i] >> stripsMerged[i] ;
  }
}

ssize_t GetStrips(std::ifstream &inFile, 
	      T_Strips &strips, ssize_t &numRead, ssize_t mergeAdjacent) {

  ssize_t *enumerations, *startRefLocations,  *endRefLocations;
  ssize_t *startQryLocation,  *endQryLocation, *stripSize;
  ssize_t *stripsMerged, length;

  ReadStrips(inFile, enumerations, startRefLocations, endRefLocations,
		       startQryLocation, endQryLocation, stripSize,
		       stripsMerged, length);

  numRead = length;
  ssize_t numStrips;
  numStrips = GetStrips(enumerations, startRefLocations, endRefLocations,
			startQryLocation, endQryLocation, stripSize,
			stripsMerged, length, strips, mergeAdjacent);

  // This information is all contained in the array strips, get rid of it.
  delete []enumerations;
  delete []startRefLocations;
  delete []endRefLocations;
  delete []startQryLocation;
  delete []endQryLocation;
  delete []stripSize;
  delete []stripsMerged;
    
  return numStrips;
}

ssize_t GetStripSize(ssize_t *strip, ssize_t pos, ssize_t hashLength) {
  if (strip == NULL)
    return hashLength;
  else
    return strip[pos];
}

ssize_t GetStripsMerged(ssize_t *merged, ssize_t pos) {
  if (merged == NULL) 
    return 0;
  else
    return merged[pos];
}

ssize_t GetMUMs(ssize_t *enumerations, ssize_t *refLocations, ssize_t *qryLocations,
	    ssize_t length, T_Strips &strips, ssize_t hashLength) {
  ssize_t numStrips = 0;
  ssize_t numMerged = 0;
  ssize_t pos       = 0;
  ssize_t stripSize = 0;
  ssize_t stripStart= 0;

  // Temporary strip
  Strip strip;

  // Process remainder of strips.
  for(pos = 1; pos <= length; pos++) {
    stripSize++;
   if (pos == length || 
       enumerations[pos] != (enumerations[pos-1] + 1) ||
       refLocations[pos] != refLocations[pos-1] + 1 ||
       qryLocations[pos] != qryLocations[pos-1] + 1) {
     // Finish storing the current strip
     // Create a strip.
     CreateStrip(stripStart, pos, stripStart+ 1, pos, 
		 enumerations, refLocations, refLocations,
		 qryLocations, qryLocations, stripSize, numMerged, strip, hashLength);
     strips.push_back(strip);
     numStrips++;

     // Reset everthing
     stripSize = 0;
     numMerged = 0;
     
     // The next strip starts at the current location.
     stripStart = pos;
   }
  }
  return numStrips;
}


ssize_t GetStrips(ssize_t *enumerations, ssize_t *startRefLocations, ssize_t *endRefLocations,
	       ssize_t *startQryLocations, ssize_t *endQryLocations, ssize_t *stripSizes, 
	      ssize_t *stripsMerged, ssize_t length, T_Strips &strips, ssize_t hashLength, ssize_t mergeAdjacent) {
  ssize_t numStrips = 0;
  ssize_t numMerged = 0;
  ssize_t pos       = 0;
  ssize_t stripSize = 0;
  ssize_t stripStart= 0;

  // Temporary strip
  Strip strip;

  // Process remainder of strips.
  for(pos = 1; pos <= length; pos++) {
    stripSize += GetStripSize(stripSizes,pos-1, hashLength);
    numMerged += GetStripsMerged(stripsMerged,pos-1);
    // Don't adjusst the length of the last tuple.
    if (pos < (length-1) && enumerations[pos] == enumerations[pos-1]+1) {
      // Adjust for overlap of k-mers
      if (pos == length -1) {
	std::cout << "adjusting strip size from " << stripSize << " by " 
		  << std::max(0, endRefLocations[pos-1] + hashLength - startRefLocations[pos])
		  << std::endl;
      }
      stripSize -= std::max(0, endRefLocations[pos-1] + hashLength - startRefLocations[pos]);
    }
    else if (pos == length || enumerations[pos] != (enumerations[pos-1] + 1) || !mergeAdjacent) {
      // Finish storing the current strip
      // Create a strip.
      // std::cout << "creating strip " << stripStart << " " << pos << std::endl;
      CreateStrip(stripStart, pos, stripStart+ 1, pos, 
		  enumerations, startRefLocations, endRefLocations, 
		  startQryLocations, endQryLocations, stripSize, 
		  numMerged, strip, hashLength);
      strips.push_back(strip); 
      numStrips++;
      
      // Reset everthing
      stripSize = 0;
      numMerged = 0;
      
      // The next strip starts at the current location.
      stripStart = pos;
    }
  }
  return numStrips;
}


void CreateStrip(ssize_t stripStart, ssize_t pos, ssize_t enumRefStart, ssize_t enumRefEnd,
		 ssize_t *enumerations, ssize_t *startRefLocations, ssize_t *endRefLocations, 
		 ssize_t *startQryLocations, ssize_t *endQryLocations, ssize_t stripSize, ssize_t numMerged, 
		 Strip &strip, ssize_t hashLength) 
{
  strip.startRefEnum = enumRefStart;
  strip.endRefEnum   = enumRefEnd;   // should be pos -1 + 1, -1 for
                                     // the index of the pr evious position
                                     // and +1 since the enumeration = position + 1

  strip.startRefPos  = startRefLocations[stripStart];
  strip.endRefPos    = endRefLocations[pos-1] + (hashLength-1);
  
  strip.startQryEnum = enumerations[stripStart];
  strip.endQryEnum   = enumerations[pos-1];
  
  strip.startQryPos  = startQryLocations[stripStart];
  strip.endQryPos    = endQryLocations[pos-1] + (hashLength-1) * Sign(strip.startQryEnum);
  
  strip.size         = stripSize;
  strip.numMerged    = numMerged;
  strip.hashLength   = hashLength;
  strip.sign = Sign(enumerations[pos-1]);
}

