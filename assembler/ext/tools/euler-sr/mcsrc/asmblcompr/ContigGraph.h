/***************************************************************************
 * Title:          ContigGraph.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef CONGIG_GRAPH_H_
#define CONGIG_GRAPH_H_

#include <vector>
#include <string>
#include <map>
#include <set>

#include "Contig.h"
#include "MUMCluster.h"
#include "MUMClusterParser.h"

typedef std::map<std::string, ssize_t> ContigNameToIndex;
typedef std::set<ssize_t> ConnectSet;
typedef std::vector<ConnectSet> ForestSet;

ssize_t FindTreeIndex(ForestSet &forest, ssize_t index);
void BuildForest(std::vector<Contig> &refContigs,
		 std::vector<Contig> &qryContigs,
		 ContigNameToIndex &qryNameToIndex,
		 ForestSet &forest);
ssize_t ValidTreeSize(ForestSet &forest,
		  ssize_t index,
		  ssize_t minTreeSize);

void CollateClusters(MUMClusterFile &clusterFile,
		     std::vector<Contig> &refContigs,
		     std::vector<Contig> &qryContigs,
		     ContigNameToIndex &refNameToIndex,
		     ContigNameToIndex &qryNameToIndex);

void ChopTree(ConnectSet &tree,
	     ssize_t lastRef,
	     std::vector<ssize_t> &refSet,
	     std::vector<ssize_t> &qrySet);
#endif
