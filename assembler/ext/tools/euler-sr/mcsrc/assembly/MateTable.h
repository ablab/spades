/***************************************************************************
 * Title:          MateTable.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MATE_TABLE_H_
#define MATE_TABLE_H_
#include <string>
#include <vector>

#include "compatibility.h"


class Clone {
public:
	_SSZT_ ai, bi;
	ssize_t type;
};

class ReadMate {
 public:
	_SSZT_ mateIndex;
	ssize_t mateType;
	char marked;
	ReadMate() {
		marked = 0;
		mateIndex = -1;
		mateType = -1;
	}

	//////// TODO: Here, marked is character and 0=unmarked, 1=marked
	//////// Whereas for GraphVertex and GraphEdge,
	//////// it's a one-bit enum, 0=Marked, 1=NotMarked

};

typedef std::pair<std::string, Clone> NameClonePair;
typedef std::vector<ReadMate> ReadMateList;


#endif
