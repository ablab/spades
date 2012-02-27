/***************************************************************************
 * Title:          MUMCluster.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MUM_CLUSTER_H_
#define MUM_CLUSTER_H_

#include <string>
#include <vector>

#include "MUMAlignBlock.h"

class MUMCluster {
public:
  std::string refName, qryName;
  ssize_t refLen, qryLen;
  std::vector<MUMAlignBlock*> blocks;
  ssize_t size() { return blocks.size(); }
  ssize_t totalSize() {
    int ts = 0;
    ssize_t b;
    for (b = 0; b < blocks.size(); b++ ){
      ts+= blocks[b]->size();
      //      std::cout << "ts: " << ts << std::endl;
    }
    return ts;
  }
  ssize_t CalculateCoverage() {
    ssize_t a, b;
    ssize_t coverage = 0;
    for (b = 0; b < blocks.size(); b++ ) {
      for (a = 0; a < blocks[b]->size(); a++ ) {
	coverage += blocks[b]->length[a];
      }
    }
    return coverage;
  }    
};

#endif
