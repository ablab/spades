/***************************************************************************
 * Title:          MUMClusterParser.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MUM_CLUSTER_PARSER_H_
#define MUM_CLUSTER_PARSER_H_

#include "MUMClusterFile.h"
#include "MUMCluster.h"
#include "MUMAlignBlock.h"

class MUMClusterParser {
public:
  static ssize_t ParseAlignBlocks(std::ifstream &inFile, MUMCluster *cluster);
  static ssize_t ParseMUMClusterFile(std::ifstream &inFile,
				 MUMClusterFile &clusterFile);
  static ssize_t ParseCluster(std::ifstream &inFile, MUMCluster *&cluster);
};

#endif
