/***************************************************************************
 * Title:          BuildPhylogeny.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <string>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "mctypes.h"
#include "utils.h"
#include "CharTree.h"
#include "InversionBins.h"
#include "InversionUtils.h"
#include "FourGameteUtils.h"
#include <set>


typedef std::set<ssize_t> IntSet;
typedef std::vector<CharNode*> NodeVector;


ssize_t EQUAL_CHAR = 0;
ssize_t CONFLICTING_CHAR = 1;

class Collection {
public:
  std::vector<IntSet> sets;
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
  }

  void SetToVector(ssize_t setNumber, std::vector<ssize_t> &vect ) {
    IntSet::iterator intit;
    for (intit = sets[setNumber].begin(); 
	 intit != sets[setNumber].end(); ++intit) {
      vect.push_back(*intit);
    }
  }

  void VectorToSet(std::vector<ssize_t> &vect, std::set<ssize_t> &set) {
    ssize_t i;
    for (i = 0; i < vect.size(); i++) {
      set.insert(vect[i]);
    }
  }

  void Intersect(ssize_t query) {
    ssize_t intersectFound = 0;
    IntVector isectIndices;
    ssize_t i, j, k;
    for (i = 0; i < sets.size() ; i++) {
      intersectFound = 0;
      if (sets[i].find(query) != sets[i].end()) {
	intersectFound = 1;
	isectIndices.push_back(i);
      }
    }
    // If no intersection was found, create a new set 
    // that includes this.  Don't need to join other sets
    // because this only involves one element.
    if (isectIndices.size() == 0) {
      ssize_t lastSet = sets.size();
      sets.resize(sets.size()+1);
      sets[lastSet].insert(query);
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
    if (outsiderFound > 0) {
      ssize_t lastSet = sets.size();
      sets.resize(sets.size()+1);
      for (i = 0; i < query.size(); i++) {
	sets[lastSet].insert(query[i]);
      }
    }
  }
};

void JoinMaximalCliques(ssize_t maxClique, 
			IntMatrix &maxCliques, IntVector &maxCliqueSizes,
			IntMatrix &similarity,
			IntVector &cliqueVertices);

ssize_t IsClique(IntVector &vertexSubset,
	     ssize_t numVertices,
	     IntVector &verticex,
	     IntMatrix &graph);

ssize_t FindMaximalCliques(IntVector &vertices,
		       IntMatrix &conflictGraph,
		       IntMatrix &maximalCliques,
		       IntVector &cliqueSizes,
		       InversionMatrix &invesions);

void CreateConflictGraph(IntMatrix &simMatrix,
			 IntMatrix &diffMatrix,
			 IntMatrix &conflictGraph);

void BreakClique(IntVector &clique,
		 ssize_t cliqueSize,
		 IntVector &vertices,
		 IntMatrix &conflictGraph,
		 Collection &brokenCliques,
		 IntVector &setClique);


void CreateConflictGraph(IntMatrix &simMatrix,
			 IntMatrix &diffMatrix,
			 IntMatrix &conflictGraph);



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

void PrintCharacters(InversionMatrix &invMat,
		     StringVector &labels,
		     IntVector &order,
		     ssize_t &iter) {
  ++iter;
  std::ofstream charOutFile;
  std::stringstream charOutStream;
  IntMatrix charMat;

  charOutStream.str("");
  charOutStream << "savedchars_";
  charOutStream << iter << ".m";
  openck(charOutStream.str(), charOutFile, std::ios::out);
  charOutFile << "chars_" << iter << " = [" << std::endl;
  invMat.PrintMiserly(charOutFile);
  charOutFile << "];" << std::endl;
  if (invMat.size() > 0) {
    charOutFile << "inv_" << iter <<  " = [" ;
    ssize_t i;
    for (i = 0; i < invMat.size()-1; i++ ) {
      charOutFile << invMat.loci[i]->species << ",";
    }
    charOutFile << invMat.loci[i]->species;
    charOutFile << "];" << std::endl;
  }
  
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

void FindLeastConflictingParents(InversionMatrix &inversions,
				 IntMatrix &parentGraph, 
				 IntMatrix &conflictGraph,
				 IntVector &minConflictingList,
				 ssize_t &minNumConflicts);

void FindInformativeChars(InversionMatrix &inversions,
			  IntVector &consideredRows,
			  IntVector &informativeIndices);

void CalculateSharedMatrix(InversionMatrix &inversions,
			   StringVector &species,
			   IntMatrix &simMatrix, double unknownRatio = 0.0) ;

void CalculateDifferingMatrix(InversionMatrix &inversions,
			      StringVector &species,
			      IntMatrix &simMatrix);


void DeleteRemoved(IntVector &toRemove,
		   NodeVector &nodes,
		   StringVector &species,
		   InversionMatrix &inversions,
		   StringVector &labels,
		   IntVector &order);

void Merge(Collection &branches, 
	   NodeVector &nodes,
	   StringVector &species,
	   InversionMatrix &inversions,
	   IntVector &toRemove);

ssize_t NotSame(InversionMatrix &inversions, IntVector &specIndices,
	    double unknownRatio, ssize_t &numNotSame);

ssize_t NotSame(InversionMatrix &inversions, ssize_t index1, ssize_t index2, ssize_t &numNotSame);

ssize_t CheckIfSame(InversionMatrix &inversions, ssize_t index1, ssize_t index2, 
	 double unknownRatio, ssize_t &same, ssize_t &numSame);

ssize_t CheckIfSame(InversionMatrix &inversions, IntVector indices, 
		double unknownRatio, ssize_t &same);

void CountCertificates(InversionMatrix &inversions,
		       IntVector &branchSpecies,
		       IntVector &certificates);

void JoinSame(Collection &samePairs,
	      Collection &sameSets);

void FindSame(InversionMatrix &invList, 
	      double unknownRatio,
	      StringVector &species,
	      Collection &sameIndices );

void GetMaxScore(IntMatrix &graph, ssize_t &specA, ssize_t &specB, ssize_t &score);
//void FindParents(InversionMatrix &invList, IntMatrix &graph, Collection &branches, int branchingSize = 2);
void GetColumns(ValidatedInversion *valInv, std::vector<ssize_t> &cols);
void MergeColumns(InversionMatrix &inversions, std::vector<ssize_t> &cols);
void DeleteColumns(InversionMatrix &inversions, std::vector<ssize_t> &toDelete);
void PrintColumns(InversionMatrix &inversions, std::vector<ssize_t> &indices);

void PrintInversions(InversionMatrix &inversions);
void RemoveRow(InversionMatrix &inversions, ssize_t pos);

void FindMostSimilar(IntMatrix &similarityMatrix, IntVector &set, 
		     ssize_t &indexA, ssize_t &indexB);
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

  std::string validInvFileName, 
    inversionFileName, 
    seqName,
    speciesFileName;
  std::string graphName;
  std::string orderFile;
  ssize_t doFourGameteCheck;
  doFourGameteCheck = 0;
  validInvFileName = "";
  inversionFileName = "";
  orderFile   = "";
  
  if (argc < 4) {
    std::cout << "usage: buildphylo binFile speciesFile graphName [-order] [-4] [-delayunknown val]" 
	      << std::endl;
    std::cout << "  -order specifies an order of species for printing labels "
	      << std::endl;
    std::cout << "  -4 : perform four-gamete check" << std::endl;
    std::cout << "  -delayunknown [value] : Do not attempt to merge two species if the fraction of known values is less than " 
	      << std::endl
	      << "       value*the number of unmerged characters. " << std::endl; 
      
    exit(1);
  }
  inversionFileName = argv[1];
  speciesFileName = argv[2];
  graphName = argv[3];
  double unknownDelay = 0.5;
  if (argc > 4) {
    int argi;
    argi = 4;
    while (argi < argc) {
      if (strcmp(argv[argi], "-delayunknown") == 0) {
	++argi;
	unknownDelay = atof(argv[argi]);
      }
      if (strcmp(argv[argi], "-order") == 0) {
	++argi;
	orderFile = argv[argi];
      }
      if (strcmp(argv[argi], "-4") == 0) {
	doFourGameteCheck = 1;
      }
      ++argi;
    }
  }
  StringVector labels, fullLabels;
  IntVector    order, fullOrder;
  StringVector species;
  InversionList invList;
  ssize_t deleted;
  // Read in the inversions and where they belong.
  ReadInversionFile(inversionFileName, species, invList);

  ReadSpeciesFile(speciesFileName, species);
  // Now try to map everything back to human, and collate orthologous inversions.
  if (orderFile != "") {
    InitializeLabelOrder(orderFile, species, labels, order);
    fullLabels = labels;
    fullOrder  = order;
  }
  else {
    ssize_t specIndex;
    fullLabels.resize(species.size());
    fullOrder.resize(species.size());
    for (specIndex = 0; specIndex < species.size(); specIndex++) {
      fullLabels[specIndex] = species[specIndex];
      fullOrder[specIndex] = specIndex + 1;
    }
  }

  std::string tableName, tempTableName, chainTableName, netTableName;
  ssize_t refLen;
  std::string seqTableName = "sequences";


  ssize_t bin,speciesi;
  /*
    for (speciesi = 0; speciesi < species.size() -1; speciesi++) 
    std::cout << species[speciesi] << " ";
    if (species.size()) 
    std::cout << species[species.size()-1] << std::endl;
  */
  ssize_t n0, n1, n2;
  ssize_t inv;
  std::vector<ssize_t> n;
 
  InversionMatrix inversions;
  ListToInvMatrix(invList, inversions);
  InversionMatrix origInversions(inversions);

  IntMatrix graph, charMat;
  ssize_t iter = 0;
  ssize_t i, j;

  CreateMatrix(graph, species.size(), species.size());

  // Print the initial char set to a file.
  //  PrintCharacters(inversions, labels, order, iter);
  PrintCharacters(origInversions, fullLabels, fullOrder, iter);

  // Branches contains the sets of species that are joined in an iteration
  Collection branches;

  bin = 0;
  ssize_t s;

  // The initial graph is a collection of n unconnected vertices.
  NodeVector nodes;
  CharNode *nodePtr;
  for (s = 0; s < species.size(); s++) { nodePtr = new CharNode;nodePtr->data = species[s]; nodes.push_back(nodePtr);  }

  ssize_t brMax, brMaxIndex;
  IntSet::iterator intIt, intEnd;


  ssize_t first = 1;

  // My method, in comments. 
  // While there remains characters to merge on and possible merges have been found
  // 1. Look for all species that have a common parent
  // 2. Join together all common parents that are linked (this creates a branch point)
  // 3. Remove any inconsistent branch points.
  // 4. Merge if possible.
  
  Collection sameIndexPairs, sameIndexSets;
  IntMatrix equalMatrix, conflictMatrix, conflictGraph;
  IntVector uniqueInversionCount;
  while (inversions.size() > 0 and  // there are character left
	 inversions.numSpecies() > 2) {

    RemoveUninformativeRows(inversions, origInversions, uniqueInversionCount);
    std::cout << "inversions size: " << inversions.size() 
	      << " on: " << inversions.numSpecies() << " species "
	      << std::endl;
    /*
      for (inv = 0; inv < inversions.size(); inv++ ){  
      std::cout << "remaining: " << inversions.loci[inv]->number << " ";
      inversions.loci[inv]->PrintOrientations(std::cout);
      std::cout << std::endl;
    }
    */
    if (doFourGameteCheck) {
      ssize_t numRemoved = RemoveMinimumConflicts(inversions, 2);
      std::cout << "Removed: " << numRemoved << std::endl;
      doFourGameteCheck = 0;
    }
    
    std::cout << "unique inversion count " << std::endl;
    ssize_t spec;
    for (spec = 0; spec < species.size(); spec++ ) {
      std::cout << species[spec] << " " << uniqueInversionCount[spec] << std::endl;;
    }
    std::cout << std::endl;

    IntersectLists(inversions,origInversions);

    // Save the current set of characters to a matlab file.
    if (inversions.size() > 0 ) {
      //      PrintInversions(inversions);
      PrintCharacters(origInversions, fullLabels, fullOrder, iter);
    }


    // Merge all distinct and defined branches
    // Look for some pair of vertices that are the same.
    sameIndexPairs.Reset();
    sameIndexSets.Reset();

    // Create graphs for resolving non-transitively equal characters
    ClearMatrix(equalMatrix);
    ClearMatrix(conflictMatrix);
    CreateMatrix(equalMatrix, species.size(), species.size());
    CreateMatrix(conflictMatrix, species.size(), species.size());

    CalculateSharedMatrix(inversions, species, equalMatrix, unknownDelay);
    CalculateDifferingMatrix(inversions, species, conflictMatrix);

    CreateConflictGraph(equalMatrix, conflictMatrix, conflictGraph);

    FindSame(inversions, unknownDelay, species, sameIndexPairs);

    JoinSame(sameIndexPairs, sameIndexSets);

    ssize_t set;
    IntVector setIndices, certificates;
    IntMatrix toMerge;
    IntVector toRemove;
    branches.Reset();
    ssize_t numNotSame;
    for (set = 0; set < sameIndexSets.size() ; set++ ) {
      std::cout << "for set: " << set << " (" ;
      setIndices.clear();
      toRemove.clear();
      sameIndexSets.SetToVector(set, setIndices);
      ssize_t c;
      for (c = 0; c < setIndices.size();c++) {
	std::cout << " " << species[setIndices[c]];
      }
      std::cout << " ) ";
      CountCertificates(inversions, setIndices, certificates);
      std::cout << "got " << certificates.size() << " certificates :";
	for (c = 0; c < certificates.size(); c++ ) {
	std::cout << " " << certificates[c];
      }
      std::cout << std::endl;
      
      // Change to see if one may merge into vertices with 0
      // inversions on the edge first.
      if (certificates.size() >= 0 ) {
	NotSame(inversions, setIndices, 0, numNotSame);
	if (numNotSame == 0 ) {
	  branches.Intersect(setIndices);
	}
	else {
	    //	  This becomes more important on larger data sets.
	  std::cout << "branch " << set << " has " << numNotSame 
		    << " inconsistencies " << std::endl;
	  // Try to find out some clique information for the set.
	  IntMatrix maximalCliques;
	  IntVector maximalCliqueSizes;
	  ssize_t maxClique;
	  Collection unbrokenCliques;
	  maxClique = FindMaximalCliques(setIndices, 
					 conflictMatrix, 
					 maximalCliques, 
					 maximalCliqueSizes,
					 inversions);
	  IntVector cliqueVertices;
	  JoinMaximalCliques(maxClique, maximalCliques, maximalCliqueSizes, equalMatrix, cliqueVertices);
	  ssize_t cl;
	  std::cout << "got clique vertices: ";
	  for (cl = 0; cl < cliqueVertices.size(); cl++) { std::cout << cliqueVertices[cl] << " "; }
	  std::cout << std::endl;
	  if (cliqueVertices.size() > 0) {
	    branches.Intersect(cliqueVertices);
	    CountCertificates(inversions, cliqueVertices, certificates);
	  }
	  std::cout << "merged clique: ";
	  for (cl = 0; cl < cliqueVertices.size(); cl++ ) {
	    std::cout << species[cliqueVertices[cl]] << " ";
	  }
	  std::cout << "got " << certificates.size() << " certificates ";
	  for (cl = 0; cl < certificates.size(); cl++) {
	    std::cout << certificates[cl] << " ";
	  }
	  std::cout << std::endl;
	}
      }
    }

    if (branches.size() > 0) {
      Merge(branches, nodes, species, inversions, toRemove);
    }
    else {
      // No more processing to do, bail out.
      break;
    }
    DeleteRemoved(toRemove, nodes, species, inversions, labels, order);
    if (inversions.size() > 0) {
      PrintCharacters(origInversions, fullLabels, fullOrder, iter);
    }
  }

  //  PrintInversions(inversions);

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

void JoinMaximalCliques(ssize_t maxClique, 
			IntMatrix &maxCliques, IntVector &maxCliqueSizes,
			IntMatrix &equalMatrix,
			IntVector &cliqueVertices) {
  ssize_t maxCliqueSize = maxCliqueSizes[maxClique];
  Collection connectedCliques;
  ssize_t c;
  ssize_t cl;
  IntVector tmpVect;
  ssize_t tmp;
  ssize_t nEqual;
  ssize_t maxNEqual;
  ssize_t maxNEqualClique;
  ssize_t s1, s2;
  maxNEqual = -1;
  maxNEqualClique = -1;
  ssize_t nMaximal = 0;
  cliqueVertices.clear();

  for (c = 0; c < maxCliques.size(); c++) {
    if (maxCliqueSizes[c] == maxCliqueSize) {
      // Step 1. Count the number of characters that this clique 
      // shares.
      nMaximal++;
      nEqual = 0;
      for (s1 = 0; s1 < maxCliqueSize-1; s1++ ) {
	for (s2 = 1; s2 < maxCliqueSize; s2++ ) {
	  nEqual += equalMatrix[maxCliques[c][s1]][maxCliques[c][s2]];
	}
      }
      std::cout << " clique: " << c << " size: " << nEqual << std::endl;
      if (nEqual > maxNEqual) {
	maxNEqual = nEqual;
	maxNEqualClique = c;
      }

      tmpVect.resize(maxCliqueSizes[c]);
      std::cout << " joining clique: ";
      for (cl = 0; cl < maxCliqueSizes[c]; cl++) { 
	std::cout << " " << maxCliques[c][cl];
	tmpVect[cl] = maxCliques[c][cl];
      }
      std::cout << std::endl;
      connectedCliques.Intersect(tmpVect);
      connectedCliques.SetToVector(0, tmpVect);
    }
  }
 
  connectedCliques.SetToVector(0, cliqueVertices);
  std::cout << "first clique size: " << connectedCliques.sets[0].size() << " of " << connectedCliques.sets.size() << std::endl;
  for (tmp = 0; tmp < cliqueVertices.size(); tmp++) 
    std::cout << cliqueVertices[tmp] << " ";
  std::cout << std::endl;

  if (cliqueVertices.size() >  maxCliqueSize) {
    std::cout << "merged conflicting cliques, break the merging" << std::endl;
    cliqueVertices.clear();
  }
}

void RemoveRow(InversionMatrix &inversions, ssize_t pos) {
  assert(pos < inversions.size());
  inversions.Erase(pos);
}

void PrintInversions(InversionMatrix &inversions) {
  ssize_t inv;
  for (inv = 0; inv < inversions.size(); inv++) {
    inversions.loci[inv]->Print(std::cout);
  }
}

void GetColumns(ValidatedInversion *valInv, std::vector<ssize_t> &cols) {
  ssize_t spec;
  for (spec = 0; spec < valInv->size(); spec++) {
    if ((*valInv)[spec] == 1)
      cols.push_back(spec);
  }
}

void PrintColumns(InversionMatrix &inversions, std::vector<ssize_t> &indices) {
  ssize_t i, j;
  ssize_t lineLength = 50;
  char alignString[51];
  alignString[50] = '\0';
  ssize_t c, ci, l;
  char alignChar;
  c = 0;
  while (c < inversions.size()) { 
    for (l = 0; l < lineLength; l++) 
      alignString[l] = ' ';
    for (i = 0; i < indices.size(); i++) {
      ci = c;
      for (l = 0; l < lineLength and ci < inversions.size(); l++, ci++) {
	for (j = i+1; j < indices.size(); j++) {
	  if (inversions.loci[ci]->orient[indices[i]] != 2 and
	      inversions.loci[ci]->orient[indices[j]] != 2 and
	      inversions.loci[ci]->orient[indices[i]] != 
	      inversions.loci[ci]->orient[indices[j]]) {
	    alignString[l] = 'X';
	  }
	}
	std::cout << inversions.loci[ci]->orient[indices[i]];
      }
      std::cout << std::endl;
    }
    alignString[l] = '\0';
    std::cout << alignString << std::endl;
    std::cout << std::endl;
    c = ci;
  }
}

void DeleteColumns(InversionMatrix &inversions, std::vector<ssize_t> &toDelete) {
  ssize_t inv,i;
  i = inversions.size() -1;
  ssize_t d = toDelete.size() - 1;
  for (inv = 0; inv < inversions.size(); inv++) {
    for (i = toDelete.size()-1; i >= 0 ; i--) {
      inversions.loci[inv]->Erase(toDelete[i]);
    }
  }
}
 

void MergeColumns(InversionMatrix &inversions, std::vector<ssize_t> &toMerge) {
  assert(toMerge.size() > 1);
  ssize_t c;
  ssize_t invFound, unknown, invSign, delFound;
  ssize_t from, to, fromVal, toVal;
  if (toMerge.size() == 0) {
    std::cout << "Nothing specified to merge " << std::endl;
    return;
  }
  to = toMerge[0];
  std::vector<ssize_t> n;
  n.resize(3);
  ssize_t i;
  /*
    std::cout << "merging: " << std::endl;
    for (i = 0; i < toMerge.size(); i++ ) {
    std::cout << toMerge[i] << " ";
    }
    std::cout << std::endl;
    PrintColumns(inversions, toMerge);
  */
  ssize_t inv;
  for (inv = 0; inv < inversions.size(); inv++) {
    // Just iterate the merging conditions
    // Case 1. Both are known
    for (i = 0; i < 3; i++ ) 
      n[i] = 0;
    
    for (i = 0; i < toMerge.size(); i++) {
      n[inversions.loci[inv]->orient[toMerge[i]]]++;
    }
    // find the inversions character
    if (n[0] > n[1])
      inversions.loci[inv]->orient[to] = 0;
    else if (n[1] > n[0])
      inversions.loci[inv]->orient[to] = 1;
    else 
      inversions.loci[inv]->orient[to] = 2;
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

void CountCertificates(InversionMatrix &inversions,
		       IntVector &branchSpecies,
		       IntVector &certificates) {
  // Find the number of times a character appears 
  // in the list of certificates for the parent.
  certificates.clear();

  // Make sure the indices in species should be sorted.
  std::sort(branchSpecies.begin(), branchSpecies.end());
  ssize_t inv, spec, b;
  for (inv = 0; inv < inversions.size(); inv++ ) {
    b = 0;
    // Look to see if inversion inv is a certificate.
    // Step 1.  All species in the branch must be in the same orientation.
    //          AND known. 
    ssize_t branchOrientation = -1;
    ssize_t numUnknown = 0;
    ssize_t numKnown = 0;
    for (b = 0; b < branchSpecies.size() ; b++ ) {
      spec = branchSpecies[b];
      if (inversions.loci[inv]->orient[spec] == 2)
	++numUnknown;
      else 
	++numKnown;
      
      if (branchOrientation == -1 and inversions.loci[inv]->orient[spec] != 2){
	branchOrientation = inversions.loci[inv]->orient[spec];
      }
      else {
	if (branchOrientation != -1 and 
	    inversions.loci[inv]->orient[spec] != 2 and 
	    inversions.loci[inv]->orient[spec] != branchOrientation) {
	  // Not all the sequences in this set are in the same 
	  // orientation, don't bother checking further.
	  branchOrientation = -1;
	  break;
	}
      }
    }
    // If the orientation of this branch cannot be determined, continue
    if (branchOrientation == -1) 
      continue;
    // The orientation of the branch may be determined.  All other species
    // must be in the opposite orientation for this to be a certificate.

    ssize_t validCertificate = 1;
    b = 0;
    for (spec = 0; spec < inversions.loci[inv]->size() and 
	   validCertificate; spec++ ) {
      if (spec != branchSpecies[b]) {
	if (inversions.loci[inv]->orient[spec] == branchOrientation) 
	  // Found something outside the branch that is in the same orientation
	  // This means that when merged the character is not yet unique, so 
	  // the species cannot be merged.  
	  validCertificate = 0;
      }
      else {
	b++;
      }
    }
    if (validCertificate and numKnown > numUnknown ) {
      certificates.push_back(inversions.loci[inv]->number);
    }
  }
}


void GVZPrintGraph(IntMatrix &graph, 
		   std::ostream &out, 
		   std::vector<std::string> vertexNames) {

  out << "graph G { " << std::endl;
  out << "size=\"8,10\";" << std::endl;
  
  out << "\"]; " << std::endl;

  out << " }" << std::endl;
}

ssize_t CheckIfSame(InversionMatrix &inversions,
		IntVector indices, 
		double unknownRatio, 
		ssize_t &same) {
  same = 0;
  if (inversions.size() == 0)
    return 0;

  ssize_t i;
  // Sanity check on the dimensions of the indices
  for (i = 0; i < indices.size(); i++ )
    assert(indices[i] < inversions[0]->size());

  // Check to make sure that everything in the entire set is the same
  ssize_t inv;
  ssize_t firstKnown;
  ssize_t numUnknown = 0;
  for (inv = 0; inv < inversions.size(); inv++ ) {
    firstKnown = -1;
    for (i = 0; i < indices.size(); i++) {
      if ((*inversions[inv])[indices[i]] != 2 and 
	  firstKnown == -1) {
	firstKnown = (*inversions[inv])[indices[i]];
      }
      else if (firstKnown != -1 and 
	       inversions[inv]->orient[indices[i]] != 2 and
	       inversions[inv]->orient[indices[i]] != firstKnown) {
	same = 0;
      }
    }
    for (i = 0; i < indices.size(); i++) {
      if (inversions[inv]->orient[indices[i]] == 2) {
	++numUnknown;
	break;
      }
    }
  }
  // erase the result if there aren't many known species
  //  std::cout << numUnknown << " " << numUnknown *1.0 / inversions.size() << " " << unknownRatio << std::endl;
  if (double(numUnknown) / inversions.size() < unknownRatio) {
    same = 0;
  }
}

ssize_t CheckIfSame(InversionMatrix &inversions, ssize_t index1, ssize_t index2, 
		double unknownRatio, ssize_t &same, ssize_t &numSame) {
  ssize_t inv;
  numSame = 0;
  if (inversions.size() <= 0) {
    return 0;
  }
  assert(index1 < inversions.numSpecies());
  assert(index2 < inversions.numSpecies());
  ssize_t numUnknown;
  numUnknown = 0;
  numSame = 0;
  same = 1;

  for (inv = 0; inv < inversions.size() and same; inv++) {
    if (inversions.loci[inv]->orient[index1] == 2 or
	inversions.loci[inv]->orient[index2] == 2)
      ++numUnknown;
    
    if (inversions.loci[inv]->orient[index1] != 2 and
	inversions.loci[inv]->orient[index1] == 
	inversions.loci[inv]->orient[index2])
      ++numSame;

    if (inversions.loci[inv]->orient[index1] != 2 and
	inversions.loci[inv]->orient[index2] != 2 and
	inversions.loci[inv]->orient[index1] !=
	inversions.loci[inv]->orient[index2] ) {
      same = 0;
      return 1;
    }
  }
  if (numUnknown >  inversions.size() * unknownRatio) {
    same = 0;
  }
  numSame = numUnknown;
  return 1;
}

ssize_t NotSame(InversionMatrix &inversions, IntVector &specIndices,
	    double unknownRatio, ssize_t &numNotSame) {
  ssize_t inv;
  ssize_t charValue;
  ssize_t index, orientation;
  ssize_t spec;
  numNotSame = 0;
  for (inv = 0; inv < inversions.size(); inv++ ) {
    orientation = 2;
    for (index = 0; index < specIndices.size() and orientation == 2; index++) {
      spec = specIndices[index];
      orientation = inversions.loci[inv]->orient[spec];
    }
    for (; index < specIndices.size(); index++ ) {
      spec = specIndices[index];
      if (inversions.loci[inv]->orient[spec] != 2 and
	  inversions.loci[inv]->orient[spec] != orientation) {
	// Found a character that is not the same in the 
	// potential ancestral sequence.  
	// Increment the number of errors, and 
	++numNotSame;
      }
    }
  }
}


ssize_t NotSame(InversionMatrix &inversions, 
	    ssize_t index1, ssize_t index2, ssize_t &numNotSame) {
  ssize_t inv;
  if (inversions.size() <= 0) {
    return 0;
  }
  assert(index1 < inversions.numSpecies());
  assert(index2 < inversions.numSpecies());
  ssize_t numUnknown;
  numNotSame = 0;
  for (inv = 0; inv < inversions.size(); inv++) {
    if (inversions.loci[inv]->orient[index1] != 2 and
	inversions.loci[inv]->orient[index2] != 2 and
	inversions.loci[inv]->orient[index1] != 
	inversions.loci[inv]->orient[index2])
      ++numNotSame;
  }
  return 1;
}

void CalculateSharedMatrix(InversionMatrix &inversions,
			   StringVector &species,
			   IntMatrix &simMatrix,
			   double unknownRatio) {

  ssize_t s1, s2;
  ssize_t same, numSame;
  CreateMatrix(simMatrix, species.size(), species.size());
  for (s1 = 0; s1 < species.size() - 1; s1++) {
    for (s2 = s1+1; s2 < species.size(); s2++ ) {
      CheckIfSame(inversions, s1, s2, unknownRatio, same, numSame);
      if (same) {
	simMatrix[s1][s2] = numSame;
	simMatrix[s2][s1] = numSame;
      }
      else {
	simMatrix[s1][s2] = 0;
	simMatrix[s2][s1] = 0;
      }
    }
  }
}

void CalculateDifferingMatrix(InversionMatrix &inversions,
			      StringVector &species,
			      IntMatrix &simMatrix) {

  ssize_t s1, s2;
  ssize_t numNotSame;
  CreateMatrix(simMatrix, species.size(), species.size());
  for (s1 = 0; s1 < species.size() - 1; s1++) {
    for (s2 = s1+1; s2 < species.size(); s2++ ) {
      NotSame(inversions, s1, s2, numNotSame);
      simMatrix[s1][s2] = numNotSame;
      simMatrix[s2][s1] = numNotSame;
    }
  }
}



void FindSame(InversionMatrix &inversions,
	      double unknownRatio,
	      StringVector &species,
	      Collection &sameIndices ) {
  ssize_t s1, s2;

  // Look for species that have the same inversions, and merge
  // them.
  
  std::vector<IntVector> mergedColumns;  
  ssize_t same, numSame;
  IntVector pair;
  pair.resize(2);
  std::cout << "starting find same with ur: " << unknownRatio << std::endl;
  for (s1 = 0; s1 < species.size() - 1; s1++) {
    for (s2 = s1+1; s2 < species.size(); s2++ ) {
      CheckIfSame(inversions, s1, s2, unknownRatio, same, numSame);
      if (same) { 
	std::cout << species[s1] << " is the exact same as " << species[s2] 
		  << " " << numSame << " " << unknownRatio * inversions.size() << std::endl;
	pair[0] = s1;
	pair[1] = s2;
	sameIndices.Intersect(pair);
      }
    }
  }
}

void JoinSame(Collection &samePairs,
	      Collection &sameSets) {
  sameSets.Reset();
  IntVector pair;
  ssize_t p;
  for (p = 0; p < samePairs.size(); p++ ) { 
    pair.clear();
    samePairs.SetToVector(p, pair);
    sameSets.Intersect(pair);
  }
}

void Merge(Collection &branches, 
	   NodeVector &nodes,
	   StringVector &species,
	   InversionMatrix &inversions,
	   IntVector &toRemove) {

  // Merge all collections of species listed by branches
  // at once.  In order to keep consistency, rather than
  // merging one collection at a time, they are processed
  // in phases.  Otherwise the indices specified by the 
  // second branch are invalidated by merging the species
  // specified by the first branch.
  
  // Pretty print what's getting merged

  // Step 1.  Create the internal vertices and link the
  // leaves to them.
  ssize_t b;
  ssize_t bi;
  IntVector branchIndices;
  for (b = 0; b < branches.sets.size(); b++ ) {
    branchIndices.clear();
    if (branches.sets[b].size() == 1) 
      continue;
	
    branches.SetToVector(b, branchIndices);
    std::cout << "merging: ";
    for (bi = 0; bi < branchIndices.size(); bi++) {
      std::cout <<  species[branchIndices[bi]] << ", ";
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
 
  // Step 2.  Merge columns into the leftmost column
  // without deleting the remaining columns.

  ssize_t i;
  // Merge columns corresponding to merged nodes.
  IntVector toMerge;
  for (b = 0; b < branches.size(); b++) {
    if (branches.sets[b].size() <= 1) 
      continue;
    toMerge.clear();
    branches.SetToVector(b, toMerge);
    MergeColumns(inversions, toMerge); 
  }
}

void DeleteRemoved(IntVector &toRemove,
		   NodeVector &nodes,
		   StringVector &species,
		   InversionMatrix &inversions,
		   StringVector &labels,
		   IntVector &order) {
  // Remove merged nodes 
  std::sort(toRemove.begin(), toRemove.end());
  ssize_t i;
  for (i = toRemove.size()-1; i >= 0 ; i--) {
    RemoveSpeciesLabelOrder(species[toRemove[i]], labels, order);
    nodes.erase(nodes.begin() + toRemove[i]);
    species.erase(species.begin() + toRemove[i]);
  }
  // Delete the non-leftmost columns
  DeleteColumns(inversions, toRemove);
}




void FindLeastConflictingParents(InversionMatrix &inversions,
				 IntMatrix &parentGraph, 
				 IntMatrix &conflictGraph,
				 IntVector &minConflictingList,
				 ssize_t &minNumConflicts) {
  std::cout << "not implemented" << std::endl;
  assert(0);
  ssize_t p;
  IntVector parentList;
  ssize_t conflicts;
  minNumConflicts = -1;
}
				
void FindInformativeChars(InversionMatrix &inversions,
			  IntVector &consideredRows,
			  IntVector &informativeIndices) {
  ssize_t inv;
  ssize_t spec, specIdx;
  IntVector n;
  n.resize(4);
  //  std::cout << "finding informative chars on: ";
  ssize_t i;
  /*
    for (i = 0; i < consideredRows.size(); i++ )
    std::cout << consideredRows[i] << " ";
    std::cout << std::endl;
  */
  for (inv = 0; inv < inversions.size(); ++inv) {
    n[0] = n[1] = n[2] = n[3] = 0;
    CountNs(inversions.loci[inv], consideredRows, n);
    if (n[0] > 1 and n[1] > 1) {
      /*
	std::cout << "inv: " << inv 
	<< " no: " << n[0] << " " << n[1] << std::endl;
      */
      informativeIndices.push_back(inv);
    }
  }
}


void CreateConflictGraph(IntMatrix &simMatrix,
			 IntMatrix &diffMatrix,
			 IntMatrix &conflictGraph) {
  
  ssize_t i, j;
  CreateMatrix(conflictGraph, simMatrix.size(), simMatrix.size());
  for (i = 0; i < simMatrix.size() - 1; i++ ) {
    conflictGraph[i][i] = EQUAL_CHAR;
    for (j = i + 1; j < simMatrix.size(); j++ ) {
      if (diffMatrix[i][j] == 0) {
	conflictGraph[i][j] = EQUAL_CHAR;
	conflictGraph[j][i] = EQUAL_CHAR;
      }
      else {
	conflictGraph[i][j] = CONFLICTING_CHAR;
	conflictGraph[j][i] = CONFLICTING_CHAR;
      }
    }
  }
  //  std::cout << "conflictGraph: " << std::endl;
  //  PrintMatrix(conflictGraph, std::cout, 2);
}


ssize_t IsClique(IntVector &vertexSubset,
	     ssize_t numVertices,
	     IntVector &vertices,
	     IntMatrix &graph) {
  ssize_t u, v;
  ssize_t i, j;
  for (i = 0; i < numVertices - 1; i++ ) {
    u = vertices[vertexSubset[i]];
    for ( j = i + 1; j < numVertices; j++ ) {
      v = vertices[vertexSubset[j]];
      /*      std::cout << "gs: " << graph.size() 
		<< " vs: size: " << vertices.size() 
		<< " u: " << u << " v: " << v
		<< " vsi: " << vertexSubset[i] 
		<< " vsj: " << vertexSubset[j] << std::endl;
      */
      if (graph[u][v] != 0) {
	return 0;
      }
    }
  }
  return 1;
}

ssize_t FindMaximalCliques(IntVector &vertices,
		       IntMatrix &conflictGraph,
		       IntMatrix &maximalCliques,
		       IntVector &cliqueSizes,
		       InversionMatrix &inversions ) {
  // Input:
  // vertices: a list of vertices to look for a clique in
  // conflictGraph: a graph of either 0 for match, or nonzero for conflict
  // 
  // Output:
  // maximalCliques: the list of other vertices in the maximal clique
  //                 each vertex is a part of.  The indices are relative
  //                 to the list of vertices, not the outside world (the vertices 
  //                 the graph.)  In other words, if vertices = [2,5,7,9] and the 
  //                 maximal clique involves 5,7,and 9, maximalCliques[1] = 1,2,3.
  // 
  // cliqueSizes   : an array of the size of every maximal clique for each vertex


  // Find the maximal nonconflicting clique for each vertex.
  
  // I'm pretty sure this will stay on small data sets, so use an
  // exact solution.  We'll see if this becomes a problem later.
  ssize_t max;
  max = 1 << vertices.size();
  
  IntVector vertexSubset;
  ssize_t binaryVertexSet, mutableVertexSet;
  ssize_t numVertices;
  ssize_t s, v, c, t, p;
  vertexSubset.resize(vertices.size());
  IntVector tempClique;
  ssize_t tempCliqueSize;
  ssize_t maxClique;

  maxClique = 0;
  CreateMatrix(maximalCliques, vertices.size(), vertices.size());
  cliqueSizes.resize(vertices.size());
  for (c = 0; c < cliqueSizes.size(); c++) {
    cliqueSizes[c] = 1;
    maximalCliques[c][0] = c;
  }
  for (binaryVertexSet = 1; binaryVertexSet < max; binaryVertexSet++) {
    // Transform the # into a set of vertices.
    numVertices = 0;
    mutableVertexSet = binaryVertexSet;
    for (s = 0; s < vertices.size(); s++ ) {
      v = mutableVertexSet & 1;
      mutableVertexSet >>= 1;
      if (v == 1) {
	vertexSubset[numVertices] = s;
	++numVertices;
      }
    }
    
    // Check to see if the set of vertices is a clique
    
    if (IsClique(vertexSubset, numVertices, vertices, conflictGraph)) {
      // Look to see if this clique beats out all known clique sizes
      // so far for the vertices in the clique.
      for (c = 0; c < cliqueSizes.size() and cliqueSizes[c] <= numVertices; c++);
    
      if (c == cliqueSizes.size()) {
	// This is the best clique for the vertices found so far
	// Remove the other cliques found for these vertices.
	for (v = 0; v < numVertices; v++ ) {
	  tempClique = maximalCliques[vertexSubset[v]];
	  tempCliqueSize = cliqueSizes[vertexSubset[v]];
	  if (tempCliqueSize < numVertices) {
	    for (t = 0; t < tempCliqueSize; t++) {
	      // Get rid of the stored clique
	      ssize_t t2;
	      maximalCliques[vertexSubset[t]][0] = tempClique[t];
	      for (t2 = 1; t2 < tempCliqueSize; t2++) {
		maximalCliques[vertexSubset[t]][t2] = 0;
	      }
	      cliqueSizes[tempClique[t]] = 1;
	    }
	  }
	}
	if (numVertices >= cliqueSizes[maxClique]) {
	  maxClique = vertexSubset[0];
	}
	// Now add the clique 
	for (v = 0; v < numVertices; v++ ) {
	  // Just in case, zero out the clique
	  for (t = 0; t < vertices.size(); t++ ) {
	    maximalCliques[vertexSubset[v]][t] = 0;
	  }
	  /*
	  std::cout << "maximal cliques size: " 
	  << maximalCliques.size() << std::endl;
	  */
	  for (t = 0; t < numVertices; t++ ) {
	    maximalCliques[vertexSubset[v]][t] = vertexSubset[t];
	  }
	  cliqueSizes[vertexSubset[v]] = numVertices;
	}
      }
    }
  }
  
  // Done finding cliques. For now, just report them.

  ssize_t maxCliqueSize = 0;
  //    std::cout << "found cliques " << std::endl;
  for (v = 0; v < vertices.size(); v++ ) {
    if (maxCliqueSize < cliqueSizes[v] )
      maxCliqueSize = cliqueSizes[v];
	
    for (t = 0; t < cliqueSizes[v]; t++ ) {
      maximalCliques[v][t] = vertices[maximalCliques[v][t]];
      //      std::cout << maximalCliques[v][t] << " ";
    }
    //    std::cout << std::endl;
    }
    //  std::cout << std::endl;
 
  if (maxCliqueSize == 2) {
    // problem! the largest clique size is 2. This makes it arbitrary to decide
    // which ones to merge, and for now I'm not merging any. 
    // One possible method is to count the number of certificates each clique
    // is good for, and removing any that are 0.
    IntVector tmpVector, certificates;
    tmpVector.resize(2);
    for (v = 0; v < vertices.size(); v++ ) {
      if (cliqueSizes[v] == 2) {
	for (t = 0; t < cliqueSizes[v]; t++ ) 
	  tmpVector[t] = maximalCliques[v][t];
	certificates.clear();
	CountCertificates(inversions, tmpVector, certificates);	
	std::cout << "clique for " << v << " has " << certificates.size() << " certificates " << std::endl;
	if (certificates.size() == 0) {
	  // there are no certificates on this clique, remove it
	  cliqueSizes[v] = 1;
	  maximalCliques[v].resize(1);
	}
      }
    }
  }
  return maxClique;
}


void BreakClique(IntVector &clique,
		 ssize_t cliqueSize,
		 IntVector &vertices,
		 IntMatrix &conflictGraph,
		 Collection &brokenCliques,
		 IntVector &remainingCliques) {

  // The clique must have some vertex that has a conflicting
  // edge, otherwise it would have already been merged. 

  // Split the clique into vertices that are connected to other 
  // vertices with conflict edges.

  ssize_t c, cv, v, cfl;
  IntVector cliqueConflictLabels;
  cliqueConflictLabels.resize(cliqueSize);
  for (c = 0; c < cliqueSize; c++ )
    cliqueConflictLabels[c] = -1;

  for (c = 0; c < cliqueSize; c++ ) {
    // Look for nodes that conflict with c.
    cv = vertices[clique[c]];
    for (v = 0; v < conflictGraph.size(); v++) {
      if (conflictGraph[cv][v] != EQUAL_CHAR) {
	// A vertex inside the clique has a conflict with another 
	// vertex.  Find what that vertex doesn't have conflicts to
	// that are inside the clique, and these must be separated 
	// from the clique (label them differently)
	for (cfl = 0; cfl < cliqueSize; cfl++ ) {
	  if (conflictGraph[vertices[clique[cfl]]][v] == EQUAL_CHAR) {
	    cliqueConflictLabels[cfl] = v;
	  }
	}
      }
    }
  }
  // Now divide the clique
  ssize_t nRemoved = 0;
  std::vector<ssize_t> unbroken;
  std::cout << "conflict labels: " << cliqueConflictLabels.size() << " ";
  for (c = 0; c < cliqueSize; c++ ) {
    std::cout << cliqueConflictLabels[c] << " ";
  }
  std::cout << std::endl;
  for (c = 0; c < cliqueSize; c++ ) {
    if (cliqueConflictLabels[c] == -1) {
      unbroken.push_back(c);
      if (remainingCliques.size() == 0)
	remainingCliques.push_back(vertices[clique[c]]);
    }
  }
  for (c = 0; c < unbroken.size(); c++ ) {
    unbroken[c]  = vertices[clique[c]];
  }
  brokenCliques.Intersect(unbroken);
  std::cout << "cliqueSize: " << cliqueSize << std::endl;
  for (c = 0; c < cliqueSize; c++ ) {
    if (cliqueConflictLabels[c] != -1) {
      brokenCliques.Intersect(vertices[clique[c]]);
      remainingCliques.push_back(vertices[clique[c]]);
    }
  }

  std::cout << "broke cliques: " << std::endl;
  brokenCliques.Print();
}


void FindMostSimilar(IntMatrix &similarityMatrix, IntVector &set, 
		     ssize_t &indexA, ssize_t &indexB) {

  indexA = -1;
  indexB = -1;
  if (set.size() < 2) {
    return;
  }

  ssize_t i, j;
  ssize_t maxSim = 0;
  for (i = 0; i < set.size() - 1; i++ ) {
    for (j = i+1; j < set.size(); j++ ){ 
      if (similarityMatrix[set[i]][set[j]] > maxSim) {
	indexA = i;
	indexB = j;
	maxSim = similarityMatrix[set[i]][set[j]];
      }
    }
  }
}
