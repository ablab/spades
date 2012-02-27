/***************************************************************************
 * Title:          Contig.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef CONTIG_H_
#define CONTIG_H_

#include "ContigMap.h"
#include "MUMCluster.h"

class Contig {
public:
  ContigMap intervals;
  ssize_t length;
  std::string name;
  ssize_t covered;
  std::vector<MUMCluster*> alignedClusters;
  void FindIntervals();
  ssize_t CalculateCoverage();
  ssize_t SumNumAlignBlocks();
};




#endif
