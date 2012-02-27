/***************************************************************************
 * Title:          ContigMap.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef CONTIG_MAP_H_
#define CONTIG_MAP_H_
#include <iostream>
#include <vector>

class ContigMap {
 public:
	ssize_t index;
	ssize_t refPos;
	ssize_t refEnd;
	ssize_t qryPos;
	ssize_t qryEnd;
	ssize_t length;
	ContigMap() {
		index = refPos = refEnd = qryPos= qryEnd = length = 0;
	}

	friend std::istream& operator>>(std::istream &in, ContigMap &map) {
		in >> map.index >> map.refPos >> map.refEnd >> map.qryPos >> map.qryEnd >> map.length;
		return in;
	}
	friend std::ostream& operator<<(std::ostream &out, const ContigMap &map){ 
		out << map.index << " " 
				<< map.refPos << " " 
				<< map.refEnd << " " 
				<< map.qryPos << " " 
				<< map.qryEnd << " "
				<< map.length;
		return out;
	}
};


void ReadUniqueContigMap(std::istream &in, std::vector<ContigMap> &map);


#endif
