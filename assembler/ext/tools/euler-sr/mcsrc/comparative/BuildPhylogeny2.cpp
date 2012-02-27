/***************************************************************************
 * Title:          BuildPhylogeny2.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <sstream>

#include "mysql/mysql.h"

#include "mctypes.h"
#include "utils.h"
#include "blocks/BlockDB.h"
#include "net/NetDB.h"
#include "blocks/dbOrthoPos.h"
#include "CharTree.h"
#include "InversionBins.h"
#include "InversionUtils.h"
#include "FourGameteUtils.h"
#include <set>

typedef std::set<ssize_t> IntSet;
typedef std::vector<CharNode*> NodeVector;

class Collection {
public:
  std::vector<IntSet> sets;
  std::vector<ssize_t> scores;
  std::vector<ssize_t> range;
  void Print() {
    ssize_t s;
    IntSet::iterator sit, send;
    for (s = 0; s < sets.size(); s++ ) {
      std::cout << "collection " << s << " : ";
      sit = sets[s].begin();
      send = sets[s].end();
      for (; sit != send; ++sit) {
	std::cout << *sit << " ";
      }
      std::cout << std::endl;
    }
  }
  void RemoveSet(ssize_t index) {
    if (index > sets.size()) 
      assert(sets.size() == 0);
    else {
      sets.erase(sets.begin() + index);
    }
  }    

  void Replace(ssize_t index, IntVector &values) {
    assert(index < sets.size());
    sets[index].clear();
    ssize_t i;
    for (i = 0; i < values.size(); i++ ){
      sets[index].insert(values[i]);
    }
  }

  ssize_t size() { return sets.size(); }
  void Reset() {
    sets.clear();
    scores.clear();
  }
  void SetToVector(ssize_t setNumber, std::vector<ssize_t> &vect ) {
    IntSet::iterator intit;
    for (intit = sets[setNumber].begin(); 
	 intit != sets[setNumber].end(); ++intit) {
      vect.push_back(*intit);
    }
  }
  void StoreRange(IntVector &vect) {
    ssize_t i;
    range.clear();
    for (i = 0; i < vect.size(); i++) { 
      range.push_back(vect[i]);
    }
    std::sort(range.begin(), range.end());
  }

  ssize_t FindSpanningBranch(std::vector<ssize_t> &query) {
    ssize_t s;
    IntVector difference, intersection;
    assert(query.size() > 0);
    std::sort(query.begin(), query.end());

    // the query must be contained in the range.

    for (s = 0; s < sets.size(); s+=2) {
      // compare this 
      std::set_difference(sets[s].begin(), sets[s].end(),
			  query.begin(), query.end(),
			  std::back_insert_iterator<IntVector>(difference));

      if (difference.size() < query.size() and
	  difference.size() > 0) {
	// A branch was found that spans one of the resolvable
	// subsets. It's ok to use then.
	return s;
      }
    }
    return -1;
  }

  ssize_t CreateDividingSet(std::vector<ssize_t> &query) {
    // assume range is sorted.
    // don't assume the query is sorted.

    // preconditions: query has some elements, and they are all contained in 
    // the range of range.
    assert(query.size() > 0);
    std::sort(query.begin(), query.end());
    ssize_t i; 

    IntVector difference, intersect;
    std::set_intersection(range.begin(), range.end(),
			  query.begin(), query.end(), 
			  std::back_insert_iterator<IntVector>(intersect));
    assert(intersect.size() > 0);

    ssize_t setFound = 0;
    ssize_t setIndex = 0;
    ssize_t s;
    for (s = 0; s < sets.size() and !setFound; s+= 2 ) {
      // Compare this 
      if (sets[s].size() == query.size()) {
	difference.clear();
	std::set_difference(sets[s].begin(), sets[s].end(),
			    query.begin(), query.end(), 
			    std::back_insert_iterator<IntVector>(difference));
	if (difference.size() == 0) {
	  setFound = 1;
	  setIndex = s;
	}
      }
      // Compare against the complement
      if (sets[s+1].size() == query.size()) {
	difference.clear();
	std::set_difference(sets[s+1].begin(), sets[s+1].end(),
			    query.begin(), query.end(), 
			    std::back_insert_iterator<IntVector>(difference));
	if (difference.size() == 0) {
	  setFound = 1;
	  setIndex = s+1;
	}
      }
    }

    if (setFound == 0) {
      sets.resize(sets.size()+2);
      VectorToSet(query, sets[sets.size()-2]);
      // Add the complement
      difference.clear();
      std::set_difference(range.begin(), range.end(),
			  query.begin(), query.end(), 
			  std::back_insert_iterator<IntVector>(difference));
      VectorToSet(difference, sets[sets.size()-1]);
      setIndex = sets.size()-2;
    }
    return setIndex;
  }
  
  void VectorToSet(std::vector<ssize_t> &vect, std::set<ssize_t> &set) {
    ssize_t i;
    for (i = 0; i < vect.size(); i++) {
      set.insert(vect[i]);
    }
  }

  void Intersect(std::vector<ssize_t> &query) {
    ssize_t intersectFound = 0;
    IntVector isectIndices;
    ssize_t i, j, k;
    for (i = 0; i < sets.size() ; i++) {
      intersectFound = 0;
      for (j = 0; j < query.size() and  !intersectFound; j++) {
	if (sets[i].find(query[j]) != sets[i].end()) {
	  intersectFound = 1;
	  isectIndices.push_back(i);
	}
      }
    }
    if (isectIndices.size() > 0) {
      // join all other sets to the first one.
      ssize_t firstSet = isectIndices[0];
      for (i = isectIndices.size()-1; i > 0; i-- ) {
	IntSet::iterator sit, send;
	sit  = sets[isectIndices[i]].begin();
	send = sets[isectIndices[i]].end();
	for (;  sit != send; ++sit) {
	  if (sets[firstSet].find(*sit) == sets[firstSet].end())
	    sets[firstSet].insert(*sit);
	}
	sets.erase(sets.begin() + isectIndices[i]);
      }
      for (i = 0; i < query.size(); i++) {
	if (sets[firstSet].find(query[i]) == sets[firstSet].end())
	  sets[firstSet].insert(query[i]);
      }
    }
    else {
      // searched through all sets and didn't find one.  Add it here
      ssize_t lastSet = sets.size();
      sets.resize(sets.size()+1);
      for (i = 0; i < query.size(); i++) {
	sets[lastSet].insert(query[i]);
      }
    }
  }

  ssize_t Increment(std::vector<ssize_t> &query, ssize_t inc=1) {
    ssize_t i, j;
    ssize_t setMatched = 0;
    ssize_t outsiderFound = 1;

    for (i = 0; i < sets.size() and outsiderFound; !outsiderFound or i++) {
      outsiderFound = 0;
      for (j = 0; j < query.size() and  !outsiderFound; j++) {
	if (sets[i].find(query[j]) == sets[i].end()) {
	  outsiderFound = 1;
	}
      }
    }
    if (outsiderFound == 0) {
      assert(i < scores.size());
      scores[i] += inc;
    }
    else {
      ssize_t lastSet = sets.size();
      sets.resize(sets.size()+1);
      for (i = 0; i < query.size(); i++) {
	sets[lastSet].insert(query[i]);
      }
      scores.push_back(inc);
    }
  }
  ssize_t GetMaxScore(ssize_t &maxScore, ssize_t &maxScoreIndex ) {
    maxScore = -1;
    ssize_t i;
    for (i = 0; i < scores.size(); i++) {
      if (scores[i] > maxScore) {
	maxScore = scores[i];
	maxScoreIndex = i;
      }
    }
    return maxScore > -1;
  }

  ssize_t GetMaxScoreIndices(ssize_t &maxScore, 
			 std::vector<ssize_t> &searchIndices, 
			 std::vector<ssize_t> &maxIndices) {
    ssize_t i,j;
    maxScore = -1;
    ssize_t maxScoreIndex = -1;
    for (i = 0; i < sets.size(); i++) {
      for (j = 0; j < searchIndices.size(); j++) 
	if (sets[i].find(searchIndices[j]) != sets[i].end())
	  if (scores[i] > maxScore) {
	    maxScore = scores[i];
	    maxScoreIndex = i;
	  }
    }
    if (maxScore >= 0) {
      SetToVector(maxScoreIndex, maxIndices);
      return 1;
    }
    else {
      return 0;
    }
  }

  void  GetMaxScoreIndices(ssize_t &maxScore, std::vector<ssize_t> &maxIndices) {
    ssize_t maxScoreIndex;
    GetMaxScore(maxScore, maxScoreIndex);
    if (maxScore <= 0) 
      return;
    ssize_t i, j;
    IntSet maxIndexSet;
    IntSet::iterator isIt, isEnd;
    // Initialize the first set
    for (i = 0; i < sets.size(); i++ ){
      if (scores[i] == maxScore) {
	for (isIt = sets[i].begin(); isIt != sets[i].end(); ++isIt) {
	  maxIndexSet.insert(*isIt);
	  maxIndices.push_back(*isIt);
	}
	break;
      }
    }
    // Look for other sets that overlap with this.
    IntSet::iterator maxIsIt, maxIsEnd;
    ssize_t overlapFound;
    for (; i < sets.size(); i++ ){
      if (scores[i] == maxScore) {
	overlapFound = 0;
	for (isIt = sets[i].begin(); isIt != sets[i].end(); ++isIt) {
	  // look to see if this set overlaps with the current max set.
	  if (maxIndexSet.find(*isIt) != maxIndexSet.end()) {
	    ssize_t ii, jj;
	    IntSet::iterator iii;
	    std::cout << "set: ";
	    for (iii = sets[i].begin(); iii != sets[i].end(); ++iii)
	      std::cout << *iii << " ";
	    std::cout << " overlaps with: " << std::endl;
	    std::cout << " max set: ";
	    for (iii = maxIndexSet.begin(); iii != maxIndexSet.end(); ++iii)
	      std::cout << *iii << " ";
	    std::cout << std::endl;

	    overlapFound = 1;
	    break;
	  }
	}
	if (overlapFound) {
	  // Add the set difference (slow set difference)
	  for (isIt = sets[i].begin(); isIt != sets[i].end(); ++isIt) {
	    // look to see if this set overlaps with the current max set.
	    if (maxIndexSet.find(*isIt) == maxIndexSet.end()) {
	      maxIndexSet.insert(*isIt);
	      maxIndices.push_back(*isIt);
	    }
	  }
	} 
      } // End checking to see if the two sets overlap.
    } // End checkign through all further sets
  }
};

void PrintLabels(StringVector &labels, ssize_t iter,
		 std::ofstream &out) {
  out << "labels_" << iter << " = { ";
  ssize_t i;
  for (i = 0; i < labels.size() - 1; i++ ) {
    out << "'" << labels[i] << "'" << ",";
  }
  if (i < labels.size()) {
    out << "'" << labels[i] << "'";
  }
  out << "};" << std::endl;
}

void PrintOrdering(IntVector &order, ssize_t iter,
		   std::ofstream &out) {
  ssize_t i;
  out << "order_" << iter << " = [ ";
  for (i = 0; i < order.size() - 1; i++ ) {
    out << order[i] << "," ;
  }
  if (i < order.size())
    out << order[i];
  out << "];" << std::endl;
}

void InitializeLabelOrder(std::string &orderFile, 
			  StringVector &species, 
			  StringVector &labels, 
			  IntVector &order) {
  
  std::ifstream orderIn;
  openck(orderFile, orderIn);
  ReadSpeciesLine(orderIn, labels);
  orderIn.close();
  ssize_t l, s;
  for (l = 0; l < labels.size(); l++ ) {
    for (s = 0; s < species.size(); s++ ) {
      if (labels[l] == species[s]) {
	order.push_back(s+1);
      }
    }
  }
}

void PrintCharacters(InversionList &invList,
		     StringVector &labels,
		     IntVector &order,
		     ssize_t iter) {
  std::ofstream charOutFile;
  std::stringstream charOutStream;
  IntMatrix charMat;
  charOutStream.str("");
  charOutStream << "savedchars_";
  charOutStream << iter << ".m";
  openck(charOutStream.str(), charOutFile, std::ios::out);
  BuildCharMat(invList, charMat);
  charOutFile << "chars_" << iter << " = [" << std::endl;
  PrintMatrix(charMat, charOutFile);
  charOutFile << "];" << std::endl;

  if (labels.size() > 0) {
    PrintLabels(labels, iter, charOutFile);
  }
  if (order.size() > 0) {
    PrintOrdering(order, iter, charOutFile);
  }
  charOutFile.close();
  charOutFile.clear();
}


void RemoveSpeciesLabelOrder(std::string &species, 
			     StringVector &labels, 
			     IntVector &order) {
  ssize_t l;
  ssize_t pos;
  for (l = 0; l < labels.size(); l++ ) {
    if (labels[l] == species) 
      break;
  }
  if (l < labels.size()) {
    pos = order[l];
    labels.erase(labels.begin() + l);
    order.erase(order.begin() + l);
    for (l = 0; l < order.size(); l++ ) {
      if (order[l] > pos) 
	order[l]--;
    }
  }
}



void FindLeastConflictingParents(InversionList &invList,
				 IntMatrix &parentGraph, 
				 IntMatrix &conflictGraph,
				 IntVector &minConflictingList,
				 ssize_t &minNumConflicts);

void FindInformativeChars(InversionList &invList,
			  IntVector &consideredRows,
			  IntVector &informativeIndices);

void CalculateSharedMatrix(InversionList &invList,
			   StringVector &species,
			   IntMatrix &simMatrix) ;

void CalculateDifferingMatrix(InversionList &invList,
			      StringVector &species,
			      IntMatrix &simMatrix);

void Merge(Collection &branches, 
	   NodeVector &nodes,
	   StringVector &species,
	   InversionList &invList,
	   StringVector &labels,
	   IntVector &order);

ssize_t NotSame(InversionList &invList, ssize_t index1, ssize_t index2, 
	    double unknownRatio, ssize_t &numNotSame);

ssize_t Same(InversionList &invList, ssize_t index1, ssize_t index2, 
	 double unknownRatio, ssize_t &same, ssize_t &numSame);

ssize_t Same(InversionList &invList, IntVector indices, 
	 double unknownRatio, ssize_t &same);


void FindSame(InversionList &invList, 
	      double unknownRatio,
	      StringVector &species,
	      Collection &sameIndices );

void MergeSame(InversionList &invList,  
	       double unknownRatio,
	       StringVector &species,
	       Collection &sameIndices,
	       Collection &conflictingIndices);

ssize_t BranchResolvable(InversionList &invList,
		     IntVector  &branchingSpecies,
		     Collection &resolvableSets);

void GetMaxScore(IntMatrix &graph, ssize_t &specA, ssize_t &specB, ssize_t &score);

void FindParents(InversionList &invList, IntMatrix &graph, Collection &branches, ssize_t branchingSize = 2);

void RemoveUninformativeRows(std::vector<ssize_t> &startPos,
			     std::vector<ssize_t> &endPos,
			     InversionList &consensus,
			     std::vector<ssize_t> &mult,
			     StringVector &species,
			     ssize_t iteration = 0);

void GetColumns(ValidatedInversion *valInv, std::vector<ssize_t> &cols);
void MergeColumns(InversionList &consensus, std::vector<ssize_t> &cols);
void DeleteColumns(InversionList &consensus, std::vector<ssize_t> &toDelete);
void PrintColumns(InversionList &consensus, std::vector<ssize_t> &indices);

void PrintInversions(std::vector<ssize_t> &startPos,
		     std::vector<ssize_t> &endPos, 
		     InversionList &consensus,
		     std::vector<ssize_t>* mult = NULL);

void RemoveRow(std::vector<ssize_t> &startPos,
	       std::vector<ssize_t> &endPos, 
	       InversionList &consensus,
	       std::vector<ssize_t> &mult, 
	       ssize_t pos);

template<typename t>
void PrintSimMat(t &scores) {
  ssize_t i, j;
  for (i = 0; i < scores.size()-1; i++) {
    for (j = 0; j < i; j++) 
      std::cout << "      ";
    for (j = i + 1; j < scores.size(); j++) {
      std::cout.width(6);
      std::cout.precision(2);
      std::cout << scores[i][j];
    }
    std::cout << std::endl;
  }
}


ssize_t graphOrder; // the order each vertex was created.

int main(int argc, char* argv[]) {

  std::string validInvFileName, binFileName, seqName;
  std::string graphName;
  std::string orderFile;
  validInvFileName = "";
  binFileName = "";
  orderFile   = "";
  
  if (argc < 3) {
    std::cout << "usage: buildphylo binFile graphName  " << std::endl;
    exit(1);
  }
  binFileName = argv[1];
  graphName = argv[2];
  if (argc > 3) {
    int argi;
    argi = 3;
    while (argi < argc) {
      if (strcmp(argv[argi], "-order") == 0) {
	++argi;
	orderFile = argv[argi];
      }
      ++argi;
    }
  }
  StringVector labels;
  IntVector    order;
  StringVector species;
  InversionList invList;
  BinMap binnedInversions;
  std::vector<ssize_t> startPos, endPos;
  ssize_t deleted;
  // Read in the inversions and where they belong.
  ReadBinFile(binFileName, species, binnedInversions);

  // Now try to map everything back to human, and collate orthologous inversions.
  if (orderFile != "") {
    InitializeLabelOrder(orderFile, species, labels, order);
  }

  std::set<std::string> tempTables;

  InversionList::iterator invIt;
  ssize_t humanStartPos, humanEndPos;
  std::string tableName, tempTableName, chainTableName, netTableName;
  ssize_t refLen;
  std::string seqTableName = "sequences";


  ssize_t bin,speciesi;
  
  
  for (speciesi = 0; speciesi < species.size() -1; speciesi++) 
    std::cout << species[speciesi] << " ";
  if (species.size()) 
    std::cout << species[species.size()-1] << std::endl;

  ssize_t n0, n1, n2;
  ssize_t inv;
  std::vector<ssize_t> n;
  // Determine the major orientation of every inversion
 
  if (validInvFileName != "" ) {
    // had to read in the validated inversions and not a bin, print the bins
    binFileName = validInvFileName + ".bins";
    std::ofstream out;
    openck(binFileName, out);
    PrintBins(species, binnedInversions, out);
  }
  
  InversionList consensus;

  std::vector<ssize_t> mult;

  startPos.clear();
  endPos.clear();
  BinMap::iterator binIt, endIt;
  for (binIt = binnedInversions.begin();
       binIt != binnedInversions.end(); ++binIt) {
    InversionList::iterator invIt;
    for (invIt = (*binIt).second.begin(); 
	 invIt != (*binIt).second.end(); 
	 ++invIt) {
      //      OneMajorCharacter(*invIt);
      consensus.push_back(*invIt);
      startPos.push_back((*invIt)->startPos);
      endPos.push_back((*invIt)->endPos);
      mult.push_back(1);
    }
  }

  IntMatrix graph, charMat;
  CreateMatrix(graph, species.size(), species.size());
  ssize_t iter = 0;

  PrintCharacters(consensus, labels, order, iter);

  Collection branches, commonBranches;

  std::ofstream consensusOut;
  consensusOut.open("consensus.out");
  PrintInversionList(consensus,consensusOut);
  consensusOut.close();

  // Do preprocessing and error removal on the data set.
  // 
  std::cout << "Initial row count " << consensus.size() << std::endl;
  RemoveUninformativeRows(startPos, endPos, consensus, mult, species, 0);
  std::cout << "rows after removing uninformative " << consensus.size() << std::endl;
  RemoveMinimumConflicts(consensus, 4);

  std::cout << "after removing conflicts: " << consensus.size() << std::endl;

  ssize_t i, j;

  PrintMatrix(graph, std::cout, 4);

  std::cout << "original matrix: " << std::endl;
  PrintInversions(startPos, endPos, consensus, &mult);
  std::cout << "merging species " << std::endl;

  bin = 0;
  ssize_t s;

  // The initial graph is a collection of n unconnected vertices.
  NodeVector nodes;
  CharNode *nodePtr;
  for (s = 0; s < species.size(); s++) {
    nodePtr = new CharNode;
    nodePtr->data = species[s];
    nodes.push_back(nodePtr);
  }

  ssize_t brMax, brMaxIndex;
  IntSet::iterator intIt, intEnd;

  std::cout << "consensus size: " << consensus.size() 
	    << " first row: " << consensus[0]->inversions.size() 
	    << std::endl;

  ssize_t first = 1;
  while (first or (consensus.size() > 0 and 
		   consensus[0]->inversions.size() > 2 and 
		   branches.size() > 0) ) {
    first = 0;
    // Merge all distinct and defined branches
    ++iter;
    IntVector toRemove;
    IntVector branchIndices;
    //    std::vector<IntVector> mergedColumns;
    Collection resolvableBranches, conflictingBranches;
    
    // Look for some pair of vertices that are the same.
    branches.Reset();

    Collection parentCollection;
    
    IntMatrix simMatrix, diffMatrix, subMatrix;
    CalculateSharedMatrix(consensus, species, simMatrix) ;
    CalculateDifferingMatrix(consensus, species, diffMatrix);
    FindParents(consensus, graph, parentCollection);
    std::cout << "similarities: " << std::endl;
    PrintMatrix(simMatrix, std::cout);
    std::cout << std::endl;
    std::cout << "differences " << std::endl;
    PrintMatrix(diffMatrix, std::cout);
    std::cout << "parents: " << std::endl;
    PrintMatrix(graph, std::cout);
    std::cout << std::endl;
    ssize_t s1,s2;
    ssize_t score;
    GetMaxScore(graph, s1, s2, score);
    IntVector leastConflictingParents;
    ssize_t minNumConflicts;
    FindLeastConflictingParents(consensus, graph, 
				diffMatrix, 
				leastConflictingParents, minNumConflicts);
    // Print the result for error checking
    ssize_t ls;
    std::cout << "got least conflicting parents: ";
    for (ls = 0; ls < leastConflictingParents.size(); ls++) {
      std::cout << leastConflictingParents[ls] << " ";
    }
    std::cout << "  numconflicts: " << minNumConflicts << std::endl;
    std::cout << "max from matrix: " << s1 << " " << s2 << " " << score << std::endl;
    std::vector<ssize_t> pair;
    pair.resize(2);
    pair[0] = s1;
    pair[1] = s2;

    branches.Intersect(leastConflictingParents);

    Merge(branches, nodes, species, consensus, labels, order);




    RemoveUninformativeRows(startPos, endPos, consensus, mult, species, iter);
      
    for (s = 0; s < species.size(); s++) 
      std::cout << species[s] << " ";
      
    std::cout << std::endl;
    std::cout << "inversions now: " << consensus.size() << std::endl;
    PrintInversions(startPos, endPos, consensus, &mult);
    std::cout << std::endl;

    if (consensus.size() > 0 ) {
      PrintCharacters(consensus, labels, order, iter);
    }
  }

  PrintInversions(startPos, endPos, consensus);

  CharTree phylogeny;
  ssize_t node;
  phylogeny.root = new CharNode;
  for (node = 0; node < nodes.size(); node++) {
    phylogeny.root->AddChild(nodes[node]);
  }
  std::ofstream phyloout;
  openck(graphName, phyloout, std::ios::out);
  phylogeny.NewickPrint(phyloout);
  phyloout.close();
}

void RemoveRow(std::vector<ssize_t> &startPos,
	       std::vector<ssize_t> &endPos, 
	       InversionList &consensus,
	       std::vector<ssize_t> &mult, 
	       ssize_t pos) {

  assert(pos < startPos.size());
  
  startPos.erase(startPos.begin() + pos);
  endPos.erase(endPos.begin() + pos);
  consensus.erase(consensus.begin() + pos);
  mult.erase(mult.begin() + pos);
}

void PrintInversions(std::vector<ssize_t> &startPos,
		     std::vector<ssize_t> &endPos, 
		     InversionList &consensus,
		     std::vector<ssize_t>* mult) {
    
  ssize_t bin = 0;
  InversionList::iterator invIt;
    
  for (invIt = consensus.begin(); 
       invIt != consensus.end(); ++invIt) {
    ssize_t inv;
    (*invIt)->startPos = startPos[bin];
    (*invIt)->endPos   = endPos[bin];
    std::cout << " ";
    std::cout.width(10);
    std::cout << startPos[bin];
    std::cout.width(10);
    std::cout << endPos[bin];
    if (mult != NULL) {
      std::cout.width(4); std::cout << (*mult)[bin] <<"\t";
    }
    ++bin;
    for (inv = 0; inv < (*invIt)->inversions.size(); ++inv) {
      std::cout << " " << (*invIt)->inversions[inv];
    }
    std::cout << std::endl; 
  }
}




  
void RemoveUninformativeRows(std::vector<ssize_t> &startPos,
			     std::vector<ssize_t> &endPos,
			     InversionList &consensus,
			     std::vector<ssize_t> &mult,
			     StringVector &species,
			     ssize_t iteration) {
  std::vector<ssize_t>::iterator startIt, endIt, multIt;
  InversionList::iterator consensusIt;
  
  std::ofstream uninformativeOut;
  ssize_t searchChar, uninformativeIndex, charNumber;
  char outName[256];
  sprintf(outName, "uninformative_%d.txt", iteration);
  openck(outName, uninformativeOut, std::ios::out);
  ssize_t s;
  for (s = 0; s < species.size(); s++ ) {
    uninformativeOut << species[s] << " ";
  }
  uninformativeOut << std::endl;
  startIt = startPos.begin();
  endIt   = endPos.begin();
  multIt  = mult.begin();
  consensusIt = consensus.begin(); 
  ssize_t i;
  std::vector<ssize_t> n;
  charNumber = 0;
  while (consensusIt != consensus.end()) {
    CountNs(*consensusIt, n);
    if (n[1] < 2 or n[0] < 2) {
      if (n[1] == 1) {
	searchChar = 1;
      }
      else if (n[0] == 1) {
	searchChar = 0;
      }

      uninformativeIndex = -1;
      for (i = 0; i < (*consensusIt)->inversions.size(); i++) {
	if ((*consensusIt)->inversions[i] == searchChar) {
	  uninformativeIndex = i;
	  break;
	}
      }
      if (uninformativeIndex != -1) {
	uninformativeOut << uninformativeIndex << "\t" << charNumber << std::endl;
      }


      // Only one N in this row, remove it.
      startIt = startPos.erase(startIt);
      endIt   = endPos.erase(endIt);
      multIt  = mult.erase(multIt);
      consensusIt = consensus.erase(consensusIt);
      // Find the species that has only 1, and output it to a file.
      searchChar = -1;
    }
    else {
      ++startIt;
      ++endIt;
      ++multIt;
      ++consensusIt;
    }
    ++charNumber;
  }
}


void GetColumns(ValidatedInversion *valInv, std::vector<ssize_t> &cols) {
  ssize_t inv;
  for (inv = 0; inv < valInv->inversions.size(); inv++) {
    if (valInv->inversions[inv] == 1)
      cols.push_back(inv);
  }
}

void PrintColumns(InversionList &consensus, std::vector<ssize_t> &indices) {
  ssize_t i, j;
  ssize_t lineLength = 50;
  char alignString[51];
  alignString[50] = '\0';
  ssize_t c, ci, l;
  char alignChar;
  c = 0;
  while (c < consensus.size()) { 
    for (l = 0; l < lineLength; l++) 
      alignString[l] = ' ';
    for (i = 0; i < indices.size(); i++) {
      ci = c;
      for (l = 0; l < lineLength and ci < consensus.size(); l++, ci++) {
	for (j = i+1; j < indices.size(); j++) {
	  if (consensus[ci]->inversions[indices[i]] != 2 and
	      consensus[ci]->inversions[indices[j]] != 2 and
	      consensus[ci]->inversions[indices[i]] != 
	      consensus[ci]->inversions[indices[j]]) {
	    alignString[l] = 'X';
	  }
	}
	std::cout << consensus[ci]->inversions[indices[i]];
      }
      std::cout << std::endl;
    }
    alignString[l] = '\0';
    std::cout << alignString << std::endl;
    
    std::cout << std::endl;
    
    c = ci;
  }
}

void DeleteColumns(InversionList &consensus, std::vector<ssize_t> &toDelete) {
  ssize_t i;
  InversionList::iterator consensusIt;
  consensusIt = consensus.begin(); 
  while (consensusIt != consensus.end()) {
    for (i = toDelete.size()-1; i >= 0 ; i--) {
      (*consensusIt)->inversions.erase((*consensusIt)->inversions.begin() + toDelete[i]);
    }
    ++consensusIt;
  }
}
 

void MergeColumns(InversionList &consensus, std::vector<ssize_t> &toMerge) {
  assert(toMerge.size() > 1);
  InversionList::iterator consensusIt;
  consensusIt = consensus.begin(); 
  ssize_t c;
  ssize_t invFound, unknown, invSign, delFound;
  ssize_t from, to, fromVal, toVal;
  to   = toMerge[0];
  std::vector<ssize_t> n;
  n.resize(3);
  ssize_t i;
  ssize_t invIndex = 0;
  std::cout << "merging: " << std::endl;
  for (i = 0; i < toMerge.size(); i++ ) {
    std::cout << toMerge[i] << " ";
  }
  std::cout << std::endl;
  PrintColumns(consensus, toMerge);
  while (consensusIt != consensus.end()) {
    // Just iterate the merging conditions
    // Case 1. Both are known
    for (i = 0; i < 3; i++ ) 
      n[i] = 0;
    
    
    for (i = 0; i < toMerge.size(); i++) {
      n[(*consensusIt)->inversions[toMerge[i]]]++;
    }

    // find the consensus character

    if (n[0] > n[1])
      (*consensusIt)->inversions[to] = 0;
    else if (n[1] > n[0])
      (*consensusIt)->inversions[to] = 1;
    else 
      (*consensusIt)->inversions[to] = 2;
    /*
      leave this up to delete columns
      for (i = 1; i < toMerge.size(); i++) {
      (*consensusIt)->inversions.erase((*consensusIt)->inversions.begin() + toMerge[i]-(i-1));
      }
    */
    ++consensusIt;
    ++invIndex;
  }
}

void GetMaxScore(IntMatrix &graph, ssize_t &specA, ssize_t &specB, ssize_t &score) {
  ssize_t maxScore = -100000;
  ssize_t i, j;   
  specA = -1; specB = -1;
  ssize_t nspec = graph.size();
  for (i = 0; i < nspec-1; i++) {
    for (j = i+1; j < nspec; j++ ) {
      if (graph[i][j] > maxScore) {
	specA = i;
	specB = j;
	maxScore = graph[i][j];
      }
    }
  }
  score = maxScore;
  assert(specA != -1 and specB != -1); 
}

void FindParents(InversionList &invList, 
		 IntMatrix &graph, 
		 Collection &branches, ssize_t branchingSize) {
  ssize_t aIndex, bIndex;
  std::vector<ssize_t> na, nb;
  
  // Perform consistency check.  This will remove some chars from
  // being considered as parents since they have question marks that 
  // cause conflicts.
  CreateMatrix(graph, 
	       invList[0]->inversions.size(),
	       invList[0]->inversions.size());
  std::vector<ssize_t> cladeIndices;
  cladeIndices.resize(branchingSize);
  ssize_t cladeIndex;
  for (aIndex = 0; aIndex < invList.size()-1; aIndex++) {
    CountNs(invList[aIndex], na);
    //    if (na[0] == 2) {
    // char a at aIndex is a potentially valid column. 
    // Check its consistency with other potentially valid columns
    // Use information in other potentially valid parents to decide if 
    /*
      for (bIndex = aIndex+1; bIndex < invList.size(); bIndex++) {
      CountNs(invList[bIndex], nb);
    */
      
    ssize_t s1, s2, i;
    ssize_t branchOrientation = -1;
    ssize_t branchSupport;
    if (na[0] == branchingSize) {
      branchOrientation = 0;
      branchSupport = na[1];
    }
    else if (na[1] == branchingSize) {
      branchOrientation = 1;
      branchSupport = na[0];
    }

    if (branchOrientation >= 0) {
      s1 = s2 = -1;
      cladeIndex = 0;
      for (i = 0; i < invList[aIndex]->inversions.size() and cladeIndex < branchingSize; i++) {
	if (invList[aIndex]->inversions[i] == branchOrientation) {
	  cladeIndices[cladeIndex] = i;
	  cladeIndex++;
	}
      }

      // Update the counts of all clade indices
      ssize_t j;
      for (i = 0; i < branchingSize; i++) {
	for (j = 0; j < branchingSize; j++) {
	  if (i != j) {
	    graph[cladeIndices[i]][cladeIndices[j]]+=1;//branchSupport;
	    graph[cladeIndices[j]][cladeIndices[i]]+=1;//branchSupport;
	  }
	}
      }
      //      branches.Increment(cladeIndices, branchSupport);
      branches.Increment(cladeIndices, 1);
    } // End check to make sure there was a branch
  } // end counting branches.
} // .

ssize_t BranchResolvable(InversionList &invList,
		     IntVector  &branchingSpecies,
		     Collection &resolvableSets) {

  // The check to see if this branch can be resolved involves 
  // looking to see if there are any distinguishing columns in 
  // this set.

  ssize_t inv, spec;
  // Make sure that the indices are in order.
  ssize_t n[4];
  std::sort(branchingSpecies.begin(), branchingSpecies.end());
  ssize_t resolvable = 0;
  ssize_t sidx;
  IntVector branchSplit;
  for (inv = 0; inv < invList.size()-1; inv++) {
    n[0] = n[1] = n[2] = n[3] = 0;
    for (spec = 0; spec < branchingSpecies.size(); spec++ ) {
      sidx = branchingSpecies[spec];
      if (invList[inv]->inversions[sidx] == 0) {
	n[0]++;
      }
      else if (invList[inv]->inversions[sidx] == 1)
	n[1]++;
      else if (invList[inv]->inversions[sidx] == 2)
	n[2]++;
      else 
	n[3]++;
    }
    if (n[1] > 0 and  n[0] > 0) {
      resolvable = 1;
      branchSplit.clear();
      //      std::cout << "found more than one n: " << n[1] << " " << n[0] << " ";
      for (spec = 0; spec < branchingSpecies.size(); spec++ ) {
	sidx = branchingSpecies[spec];
	//	std::cout << invList[inv]->inversions[sidx] << " ";
	if (invList[inv]->inversions[sidx] == 0 ||
	    invList[inv]->inversions[sidx] == 1)
	  branchSplit.push_back(sidx);
      }
      //      std::cout << std::endl;
      
      resolvableSets.CreateDividingSet(branchSplit);
    }
  }
  return resolvable;
}

void GVZPrintGraph(IntMatrix &graph, 
		   std::ostream &out, 
		   std::vector<std::string> vertexNames) {

  out << "graph G { " << std::endl;
  out << "size=\"8,10\";" << std::endl;
  
  out << "\"]; " << std::endl;

  out << " }" << std::endl;
}

ssize_t Same(InversionList &invList, IntVector indices, double unknownRatio, ssize_t &same) {
  if (invList.size() <= 0) {
    return 0;
  }
  ssize_t i;
  // Sanity check on the dimensions of the indices
  for (i = 0; i < indices.size(); i++ )
    assert(indices[i] < invList[0]->inversions.size());

  // Check to make sure that everything in the entire set is the same
  ssize_t inv;
  ssize_t firstKnown;
  ssize_t numUnknown = 0;
  for (inv = 0; inv < invList.size(); inv++ ) {
    firstKnown = -1;
    for (i = 0; i < indices.size(); i++) {
      if (invList[inv]->inversions[indices[i]] != 2 and 
	  firstKnown == -1) {
	firstKnown = invList[inv]->inversions[indices[i]];
      }
      else if (firstKnown != -1 and 
	       invList[inv]->inversions[indices[i]] != 2 and
	       invList[inv]->inversions[indices[i]] != firstKnown) {
	same = 0;
      }
    }
    for (i = 0; i < indices.size(); i++) {
      if (invList[inv]->inversions[i] == 2) {
	++numUnknown;
	break;
      }
    }
  }
}

ssize_t Same(InversionList &invList, ssize_t index1, ssize_t index2, 
	 double unknownRatio, ssize_t &same, ssize_t &numSame) {
  ssize_t inv;
  if (invList.size() <= 0) {
    return 0;
  }
  assert(index1 < invList[0]->inversions.size());
  assert(index2 < invList[0]->inversions.size());
  ssize_t numUnknown;
  numUnknown = 0;
  numSame = 0;
  same = 1;
  for (inv = 0; inv < invList.size() and same; inv++) {
    if (invList[inv]->inversions[index1] == 2 or
	invList[inv]->inversions[index2] == 2)
      ++numUnknown;
    
    if (invList[inv]->inversions[index1] != 2 and
	invList[inv]->inversions[index1] == 
	invList[inv]->inversions[index2])
      ++numSame;

    if (!(invList[inv]->inversions[index1] == 2 or
	  invList[inv]->inversions[index2] == 2 or 
	  invList[inv]->inversions[index1] ==
	  invList[inv]->inversions[index2] ) )
      same = 0;
  }
  //  std::cout << "nu: " << numUnknown << " ils: " << invList.size() << std::endl;
  if ( double(numUnknown) / invList.size() > unknownRatio) 
    same = 0;

  return 1;
}

ssize_t NotSame(InversionList &invList, ssize_t index1, ssize_t index2, 
	    double unknownRatio, ssize_t &numNotSame) {
  ssize_t inv;
  if (invList.size() <= 0) {
    return 0;
  }
  assert(index1 < invList[0]->inversions.size());
  assert(index2 < invList[0]->inversions.size());
  ssize_t numUnknown;
  numUnknown = 0;
  numNotSame = 0;
  for (inv = 0; inv < invList.size(); inv++) {
    if (invList[inv]->inversions[index1] != 2 and
	invList[inv]->inversions[index2] != 2 and
	invList[inv]->inversions[index1] != 
	invList[inv]->inversions[index2])
      ++numNotSame;
  }
  return 1;
}

void CalculateSharedMatrix(InversionList &invList,
			   StringVector &species,
			   IntMatrix &simMatrix) {

  ssize_t s1, s2;
  ssize_t same, numSame;
  CreateMatrix(simMatrix, species.size(), species.size());
  for (s1 = 0; s1 < species.size() - 1; s1++) {
    for (s2 = s1+1; s2 < species.size(); s2++ ) {
      Same(invList, s1, s2, 0, same, numSame);
      simMatrix[s1][s2] = numSame;
      simMatrix[s2][s1] = numSame;
    }
  }
}

void CalculateDifferingMatrix(InversionList &invList,
			      StringVector &species,
			      IntMatrix &simMatrix) {

  ssize_t s1, s2;
  ssize_t numNotSame;
  CreateMatrix(simMatrix, species.size(), species.size());
  for (s1 = 0; s1 < species.size() - 1; s1++) {
    for (s2 = s1+1; s2 < species.size(); s2++ ) {
      NotSame(invList, s1, s2, 0, numNotSame);
      simMatrix[s1][s2] = numNotSame;
      simMatrix[s2][s1] = numNotSame;
    }
  }
}


void FindSame(InversionList &invList, 
	      double unknownRatio,
	      StringVector &species,
	      Collection &sameIndices ) {
  std::cout << "inside merge same " << std::endl;
  ssize_t s1, s2;

  // Look for species that have the same inversions, and merge
  // them.
  
  std::vector<IntVector> mergedColumns;  
  ssize_t same, numSame;
  IntVector pair;
  pair.resize(2);
  
  for (s1 = 0; s1 < species.size() - 1; s1++) {
    for (s2 = s1+1; s2 < species.size(); s2++ ) {
      Same(invList, s1, s2, unknownRatio, same, numSame);
      if (same) { 
	std::cout << species[s1] << " is the exact same as " << species[s2] 
		  << " " << numSame << std::endl;

	pair[0] = s1;
	pair[1] = s2;
	sameIndices.Intersect(pair);
      }
    }
  }
}

void MergeSame(InversionList &invList, 
	       double unknownRatio,
	       StringVector &species,
	       Collection &sameIndices,
	       Collection &conflictingIndices) {

  // Merge all collections of equivalent characters
  ssize_t c;
  ssize_t setSame;
  ssize_t same, numSame;

  for (c = sameIndices.size() - 1; c >= 0 ; c-- ){
    IntSet::iterator sit, send;
    if (sameIndices.sets[c].size() > 1) {
      std::cout << "got the same indices: " << sameIndices.sets[c].size() << " ";
      if (sameIndices.sets[c].size() > 2) {
	// Try a consistency check for a same set larger than 2
	std::cout << "checking set of size " << sameIndices.sets[c].size() << std::endl;
	IntVector sameIndexVect;
	sameIndices.SetToVector(c, sameIndexVect);
	Same(invList, sameIndexVect, 1, same);
	if (same) {
	  std::cout << "all checked out consistent " << std::endl;
	}
	else {
	  std::cout << "had at least one inconsistency " << std::endl;
	  IntVector conflictingSet;
	  sameIndices.SetToVector(c, conflictingSet);
	  ssize_t is; 
	  for (is = 0; is < conflictingSet.size(); is++ ) {
	    std::cout << conflictingSet[is] << " ";
	  }
	  std::cout << std::endl;
	  conflictingIndices.Intersect(conflictingSet);
	  std::cout << "there are : " << conflictingIndices.size() 
		    << " setf of ci's " << std::endl;
	  sameIndices.RemoveSet(c);
	}
      }
      for (sit = sameIndices.sets[c].begin();
	   sit != sameIndices.sets[c].end();
	   ++sit) {
	std::cout << *sit << " ";
      }
      std::cout << std::endl;
    }
  }
}



void Merge(Collection &branches, 
	   NodeVector &nodes,
	   StringVector &species,
	   InversionList &invList,
	   StringVector &labels,
	   IntVector &order) {
  // Pretty print what's getting merged
  std::cout << "merging: ";
      
  ssize_t b;
  ssize_t bi;
  IntVector branchIndices, toRemove;
  for (b = 0; b < branches.sets.size(); b++ ) {
    branchIndices.clear();
    if (branches.sets[b].size() == 1) 
      continue;
	
    branches.SetToVector(b, branchIndices);
	  
    for (bi = 0; bi < branchIndices.size(); bi++) {
      std::cout << branchIndices[bi] << " " 
		<<  species[branchIndices[bi]] << " ";
    }
    std::cout << std::endl;

    // Create the ancestral vertex
    CharNode *nodePtr;
    nodePtr = new CharNode;
    for (bi = 0; bi < branchIndices.size(); bi++) {
      nodePtr->AddChild(nodes[branchIndices[bi]]);
    }
      
    // Link the branches back into the graph
    nodes[branchIndices[0]] = nodePtr;
    std::stringstream sstr;
    sstr << graphOrder;
    nodePtr->data = sstr.str();
    ++graphOrder;

    // Store what to remove
    for (bi = 1; bi < branchIndices.size(); bi++) {
      toRemove.push_back(branchIndices[bi]);
    }
  }
  ssize_t i;
  // Merge columns corresponding to merged nodes.
  IntVector toMerge;
  for (b = 0; b < branches.size(); b++) {
    if (branches.sets[b].size() <= 1) 
      continue;
    toMerge.clear();
    branches.SetToVector(b, toMerge);
    MergeColumns(invList, toMerge); 
  }
      
  // Remove merged nodes
  std::sort(toRemove.begin(), toRemove.end());
  for (i = toRemove.size()-1; i >= 0 ; i--) {
    RemoveSpeciesLabelOrder(species[toRemove[i]], labels, order);
    nodes.erase(nodes.begin() + toRemove[i]);
    species.erase(species.begin() + toRemove[i]);
  }
    
  // Delete unneeded columns
  DeleteColumns(invList, toRemove);
}


void FindLeastConflictingParents(InversionList &invList,
				 IntMatrix &parentGraph, 
				 IntMatrix &conflictGraph,
				 IntVector &minConflictingList,
				 ssize_t &minNumConflicts) {
  Collection parents;
  FindParents(invList, parentGraph, parents, 2);  
  ssize_t p;
  IntVector parentList;
  ssize_t conflicts;
  minNumConflicts = -1;
  minConflictingList.clear();
  for (p = 0; p < parents.size(); p++) {
    parentList.clear();
    parents.SetToVector(p, parentList);
    conflicts = 0;
    ssize_t pi;
    for (pi = 0; pi < parentList.size()-1; pi++ ) {
      conflicts += conflictGraph[parentList[pi]][parentList[pi+1]];
    }
    if (minNumConflicts < 0 or minNumConflicts > conflicts) {
      minNumConflicts = conflicts;
      minConflictingList = parentList;
    }	  
  }
}
				
void FindInformativeChars(InversionList &invList,
			  IntVector &consideredRows,
			  IntVector &informativeIndices) {
  ssize_t inv;
  ssize_t spec, specIdx;
  IntVector n;
  n.resize(4);
  std::cout << "finding informative chars on: ";
  ssize_t i;
  for (i = 0; i < consideredRows.size(); i++ )
    std::cout << consideredRows[i] << " ";
  std::cout << std::endl;
  for (inv = 0; inv < invList.size(); ++inv) {
    n[0] = n[1] = n[2] = n[3] = 0;
    CountNs(invList[inv], consideredRows, n);
    if (n[0] > 1 and n[1] > 1) {
      std::cout << "inv: " << inv << " no: " << n[0] << " " << n[1] << std::endl;
      informativeIndices.push_back(inv);
    }
  }
}



