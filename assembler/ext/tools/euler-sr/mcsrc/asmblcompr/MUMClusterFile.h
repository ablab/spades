/***************************************************************************
 * Title:          MUMClusterFile.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MUM_CLUSTER_FILE
#define MUM_CLUSTER_FILE

#include <vector>
#include <string>
#include "MUMCluster.h"

class MUMClusterFile {
public:
  std::string refFileName, qryFileName, method;
  std::vector<MUMCluster*> clusters;
  ssize_t size() { 
    return clusters.size(); 
  }
};

#endif
