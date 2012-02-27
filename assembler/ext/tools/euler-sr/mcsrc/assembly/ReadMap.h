/***************************************************************************
 * Title:          ReadMap.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "utils.h"
#include <vector>

class ReadMap {
 public:
	std::string name;
	//	std::vector<int> map;
	ssize_t start;
	ssize_t end;
};

typedef std::vector<ReadMap> ReadMapList;

void ReadReadMapFile(std::string &mapFileName, 
										 ReadMapList &readMap);

ssize_t ReadMappedReads(std::string readFileName, 
										ReadMapList &readMap,
										DNASequenceList &reads);


