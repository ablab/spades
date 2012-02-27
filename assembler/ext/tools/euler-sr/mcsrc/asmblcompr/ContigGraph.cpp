/***************************************************************************
 * Title:          ContigGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ContigGraph.h"

ssize_t FindTreeIndex(ForestSet &forest, ssize_t index) {
  ssize_t t;
  for (t = 0; t < forest.size(); t++ ){
    if (forest[t].find(index) != forest[t].end()) {
      /*            std::cout << " found " << index << " in " << t << ": ";
										ssize_t t2;
										std::set<ssize_t>::iterator sit, send;
										for (sit = forest[t].begin(); sit != forest[t].end(); sit++) {
										std::cout << " " << *sit;
										}
										std::cout << std::endl;*/
      return t;
    }
  }
  return -1;
}


// Function BuildForest
// Purpose: Each ref contig is aligned to one or more query contigs.
// Consider all alignments, and build a forest where each connected
// component (tree) contains the indices of all contigs connected by
// an alignment.
//
void BuildForest(std::vector<Contig> &refContigs,
								 std::vector<Contig> &qryContigs,
								 ContigNameToIndex &qryNameToIndex,
								 ForestSet &forest) {
  ssize_t r, q;
  ssize_t qryIndex, treeRefIndex, treeQryIndex, treeIndex;
  // Create a list of sets that have the indices of contigs in each set.
  // For each contig on the reference genome, if it is already in a
  // tree, assign everything it is connected to to that tree.  
  // Otherwise, take the first query contig that is assigned a tree,
  // and place all other connected components in that tree.
  // If neither the reference nor any of the query contigs it is
  // aligned to is in a tree, create a new tree.

  //  std::cout << "buiding a forest " << std::endl;
  for (r = 0; r < refContigs.size(); r++ ) {
		// Step 1, find a tree index
    treeIndex = -1;
    treeRefIndex = FindTreeIndex(forest, r);
    //    std::cout << "index for ref: " << r << " " << treeRefIndex << std::endl;
    if (treeRefIndex >= 0) {
      treeIndex = treeRefIndex;
    }
    for (q = 0; q < refContigs[r].alignedClusters.size() and treeIndex != -1; q++) {
      qryIndex = qryNameToIndex[refContigs[r].alignedClusters[q]->qryName] + refContigs.size();
      treeQryIndex = FindTreeIndex(forest, qryIndex);
      //      std::cout << "qry: " << qryIndex << " has treeQryIndex: " << treeQryIndex << std::endl;
      if (treeQryIndex >= 0) {
				treeIndex = treeQryIndex;
      }
    }
    //    std::cout << "tree index: " << treeIndex << std::endl;
    // Step 2. If a tree contains one of the aligned contig indice,
    // place everything else in that tree.
    
    if (treeIndex == -1) {
      treeIndex = forest.size();
      forest.resize(treeIndex+1);
      //      std::cout << "creating a new tree " << treeIndex << std::endl;
    }
    // 2.1 if the treeRefIndex >= 0, then we are using the reference
    // index.  Otherwise we are using a query index.

    if (treeRefIndex == -1) {
      treeRefIndex = treeIndex;
      forest[treeIndex].insert(r);
      //      std::cout << "adding ref " << r << " to " << treeIndex << std::endl;
    }
    for (q = 0; q < refContigs[r].alignedClusters.size() and treeIndex != -1; q++) {
      qryIndex = qryNameToIndex[refContigs[r].alignedClusters[q]->qryName] + refContigs.size();
      treeQryIndex = FindTreeIndex(forest, qryIndex);
      //      std::cout  << "looking to add " << qryIndex << " to " << treeQryIndex << std::endl;
      if (treeQryIndex == -1) {
				// This query contig has not been placed in a tree (connected component)
				forest[treeIndex].insert(qryIndex);
				//	std::cout << " adding qry " << qryIndex << " to " << treeIndex << std::endl;
      }
      else if (treeQryIndex >= 0 and treeQryIndex != treeIndex) {
				// This query contig is in a tree, but in a different tree.
				// This is the case when the reference is aligned to a query
				// node that has already been aligned to something else. 
				// Add everything in the other tree to this one, and remove
				// the other tree.
				forest[treeIndex].insert(forest[treeQryIndex].begin(),
																 forest[treeQryIndex].end());
				/*
					std::cout << "adding everything in " << treeQryIndex << " (" 
					<< forest[treeQryIndex].size() << ") to " 
					<< treeIndex << std::endl; 
				*/
				// erase this tree
				forest[treeQryIndex].clear();
      }
      // else this query contig is already in a tree, and we are
      // adding everyhing else to this tree, so nothing needs to be
      // done for this contig alone.
    }
  }   
}

void ChopTree(ConnectSet &tree,
						 ssize_t lastRef,
						 std::vector<ssize_t> &refSet,
						 std::vector<ssize_t> &qrySet) {
	//UNUSED//  int i;
  ConnectSet::iterator cit;
  for (cit = tree.begin(); cit != tree.end(); cit++) {
    if (*cit < lastRef)
      refSet.push_back(*cit);
    else
      qrySet.push_back(*cit - lastRef);
  }
}

ssize_t ValidTreeSize(ForestSet &forest,
									ssize_t index,
									ssize_t minTreeSize) {
  if (minTreeSize == 0) 
    return 1;
  ssize_t treeIndex;
  treeIndex = FindTreeIndex(forest, index);
  // Every ref contig should be in a connected component
  assert(treeIndex != -1);
  if (forest[treeIndex].size() > minTreeSize) {
    return 1;
  }
  else {
    return 0;
  }
}

void  CollateClusters(MUMClusterFile &clusterFile,
											std::vector<Contig> &refContigs,
											std::vector<Contig> &qryContigs,
											ContigNameToIndex &refNameToIndex,
											ContigNameToIndex &qryNameToIndex) {
    
  ssize_t i;
  ssize_t refContigIndex, qryContigIndex;
  for (i = 0; i < clusterFile.size(); i++) {
    // Collate the reference clusters.
    if (refNameToIndex.find(clusterFile.clusters[i]->refName) == 
				refNameToIndex.end()) {
      refNameToIndex[clusterFile.clusters[i]->refName] = refContigs.size();
      refContigIndex = refContigs.size();
      refContigs.resize(refContigs.size() + 1);
      refContigs[refContigs.size() - 1].length = clusterFile.clusters[i]->refLen;
      refContigs[refContigs.size() - 1].name = clusterFile.clusters[i]->refName;
    }
    else 
      refContigIndex = refNameToIndex[clusterFile.clusters[i]->refName];
    refContigs[refContigIndex].alignedClusters.push_back(clusterFile.clusters[i]);
      
    // Collate the query clusters
    if (qryNameToIndex.find(clusterFile.clusters[i]->qryName) == 
				qryNameToIndex.end()) {
      qryNameToIndex[clusterFile.clusters[i]->qryName] = qryContigs.size();
      qryContigIndex = qryContigs.size();
      qryContigs.resize(qryContigs.size() + 1);
      qryContigs[qryContigs.size() - 1].length = clusterFile.clusters[i]->qryLen;
      qryContigs[qryContigs.size() - 1].name   = clusterFile.clusters[i]->qryName;
    }
    else 
      qryContigIndex = qryNameToIndex[clusterFile.clusters[i]->qryName];
    qryContigs[qryContigIndex].alignedClusters.push_back(clusterFile.clusters[i]);
  }

  for (i = 0; i < refContigs.size(); i++ ) {
    refContigs[i].FindIntervals();
    refContigs[i].CalculateCoverage();
  }

  for (i = 0; i < qryContigs.size(); i++ ) {
    qryContigs[i].FindIntervals();
    qryContigs[i].CalculateCoverage();
  }
}


