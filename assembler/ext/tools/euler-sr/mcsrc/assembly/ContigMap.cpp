/***************************************************************************
 * Title:          ContigMap.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ContigMap.h"

void ReadUniqueContigMap(std::istream &in, std::vector<ContigMap> &map) {
	ContigMap m, empty;
	//UNUSED// ssize_t i;
	ssize_t index, count;
	while(in) {
		in >> index >> count;
		if (count >0) {
			in >> m;
			map.push_back(m);
		}
		else {
			map.push_back(empty);
		}
	}
}
