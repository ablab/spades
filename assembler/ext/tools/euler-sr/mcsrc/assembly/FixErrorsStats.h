/***************************************************************************
 * Title:          FixErrorsStats.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/05/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef FIX_ERRORS_STATS_H_
#define FIX_ERRORS_STATS_H_

class Stats {
public:
  ssize_t numEdge;
  ssize_t numIns;
  ssize_t numDel;
  ssize_t numMut;
  ssize_t numNotSolid;
	ssize_t numTies;
  ssize_t numNoPathFound;
  ssize_t numMultiplePaths;
	ssize_t numErrorAtEnd;
  Stats() {
    Reset();
  }
  void Reset() {
    numEdge = numIns = numDel = numMut = 0;
    numNotSolid = 0;
    numNoPathFound = 0;
    numMultiplePaths = 0;
		numErrorAtEnd = 0;
		numTies = 0;
  }
	Stats &Append(Stats &s) {
    numEdge += s.numEdge;
    numIns += s.numIns;
    numDel += s.numDel;
    numMut += s.numMut;
		numTies += s.numTies;
    numNotSolid += s.numNotSolid;
    numNoPathFound  += s.numNoPathFound;
    numMultiplePaths += s.numMultiplePaths;
		numErrorAtEnd += s.numErrorAtEnd;
    return *this;
	}
  Stats &operator+=(Stats &s) {
		Append(s);
		return *this;
  }
	friend std::ostream &operator<<(std::ostream &strm, const Stats &s) {
		strm << s.numIns << "\t" 
				 << s.numDel << "\t" 
				 << s.numMut << "\t" 
				 << s.numNotSolid << "\t" 
				 << s.numNoPathFound << "\t" 
				 << s.numMultiplePaths << std::endl;
				return strm;
	}
	static void PrintHeader(std::ostream &strm) {
		strm << "ins\tdel\tnumMut\tbad\tnoPath\tmultiplePath" << std::endl;
	}
};


#endif
