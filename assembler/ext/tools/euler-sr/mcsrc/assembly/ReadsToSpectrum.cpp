/***************************************************************************
 * Title:          ReadsToSpectrum.cpp 
 * Author:         Mark Chaisson, Glenn Tesler
 * Created:        2008
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"
#include "IntegralTupleStatic.h"
#include "utils.h"

#include <vector>
#include <string>

class CompareWithCombine {
public:
	ssize_t operator()(CountedIntegralTuple &a, CountedIntegralTuple &b) const {
		if (a == b) { // compare tuples
		//		if (a.tuple == b.tuple) {
			//			a.count = b.count = a.count + b.count;
			a.AdvanceMult(b.count);
			b.count = a.count;
		}
		//		return a.tuple < b.tuple;
		return a < b; // compare tuples
	}
};

class CompareCountedIntegralTuples {
public:
	ssize_t operator()(const CountedIntegralTuple &a, const CountedIntegralTuple &b) const {
		//		return a.tuple < b.tuple;
		return a < b; // compare tuples
	}
};

using namespace std;

typedef std::vector<CountedIntegralTuple> TupleVector;

// debugging only
void DumpTupleList(TupleVector &tupleList) {
	std::string tupStr;
	for (size_t j=0; j < tupleList.size(); j++) {
		tupleList[j].ToString(tupStr);
		std::cout << "j=" << j << "  count=" << tupleList[j].count << "  tuple=" << tupStr << std::endl;
	}
}

void CompressReverseComplements(TupleVector &tupleList) {
	//	std::cout << "CompressReverseComplements - before sort" << std::endl; // DEBUG
	//	DumpTupleList(tupleList);

	std::sort(tupleList.begin(), tupleList.end(), CompareCountedIntegralTuples());

	//	std::cout << "CompressReverseComplements - after sort" << std::endl; // DEBUG
	//	DumpTupleList(tupleList); // DEBUG
	//	std::cout << std::endl; // DEBUG

	ssize_t i;
	CountedIntegralTuple rc;
	i = 0;

	// Combine rev comp multiplicities.
	for (i = 0; i < tupleList.size(); i++ ) {
		if (tupleList[i].count < 0) {
			continue;
		}
		tupleList[i].MakeRC(rc);



		//		if (rc.tuple > tupleList[i].tuple) {
		if (rc > tupleList[i]) { // compare tuples
			ssize_t rcTupleIndex;
			rcTupleIndex = LookupBinaryTuple(&tupleList[i+1], tupleList.size() - i - 1, rc);
			if (rcTupleIndex != -1) {
				//				tupleList[i].count = tupleList[i].count + tupleList[i+1+rcTupleIndex].count;
				tupleList[i].AdvanceMult(tupleList[i+1+rcTupleIndex].count);
				tupleList[rcTupleIndex + i + 1].count = -1;
			}
		}
	}

	ssize_t cur;
	for (i = 0, cur = 0; i < tupleList.size(); i++) {
		if (tupleList[i].count >= 0) {
			tupleList[cur] = tupleList[i];

			//			tupleList[cur].CopyTuple(tupleList[i]);
			//			tupleList[cur].count = tupleList[i].count;
			
			//			tupleList[cur].tuple = tupleList[i].tuple;
			//			tupleList[cur].count = tupleList[i].count;
			cur++;
		}
	}
	tupleList.resize(cur);
}

void Compress(TupleVector &tupleList) {
	std::sort(tupleList.begin(), tupleList.end(), CompareCountedIntegralTuples());
	//UNUSED+// ssize_t n;
	ssize_t i, c ;
	i = 0; c = 0;
	ssize_t numRemoved = 0;
	ssize_t indexOfLast = tupleList.size() - 1;
	ssize_t listSize = tupleList.size();
	ssize_t curI;
	ssize_t count;
	for (i = 0; i < listSize ; i++ ){ 
		count = 1;
		curI = i;
		while (i < indexOfLast and 
					 tupleList[curI] == tupleList[i+1]) {
			i++;
			count++;
		}
		tupleList[c] = tupleList[i];
		//		tupleList[c].count = count;
		tupleList[c].SetMult(count);
		c++;
		numRemoved += (count - 1);
	}
	
	//	cout << "Removed: " << numRemoved << " during compression." << endl;
	tupleList.resize(c);
}


void Join(TupleVector &srcA, 
					TupleVector &srcB, 
					TupleVector &dest) {
	ssize_t a = 0, b = 0;

	ssize_t nA = srcA.size(), nB = srcB.size();
	//UNUSED// ssize_t d;

	//	back_insert_iterator<TupleVector> mergeOutput(dest);
	//	set_union(srcA.begin(), srcA.end(), srcB.begin(), srcB.end(), mergeOutput, CompareWithCombine() );
	while (a < nA and b < nB) {
		if (srcA[a] < srcB[b]) {
			dest.push_back(srcA[a]);
			a++;
		}
		else if (srcB[b] < srcA[a]) {
			dest.push_back(srcB[b]);
			b++;
		}
		else {
			assert(srcA[a] == srcB[b]); // compare tuples
			//			assert(srcA[a].tuple == srcB[b].tuple);
			//UNUSED// int prevA = srcA[a].count;
			//			srcA[a].count = srcA[a].count + srcB[b].count;
			srcA[a].AdvanceMult(srcB[b].count);
			dest.push_back(srcA[a]);
			//UNUSED// int lastd = dest.size() - 1;
			a++;
			b++;
		}
	}
	while (a < nA) {
		dest.push_back(srcA[a]);
		a++;
	}
	while (b < nB) {
		dest.push_back(srcB[b]);
		b++;
	}

	//	CompressReverseComplements(dest);
}

void Update(TupleVector &src,
						TupleVector &dest) {
	//UNUSED// ssize_t nSrc = src.size(), nDest = dest.size();
	
	TupleVector tmp;
	Compress(src);
	Join(src,dest,tmp);
	dest = tmp;
	/*	cout << "Update: " << nSrc << " " << nDest << " " << src.size() 
			 << " " << dest.size() << endl;
	*/
}


// Define data structures and operations for a balanced binomial sorted list tree

class BinomialVertexNode {
public:
	ssize_t totalSize;
	BinomialVertexNode *left, *right, *parent;
	ssize_t leftSize, rightSize;
	TupleVector tupleList;
	string pagedOutFileName;
	ssize_t isPaged;
	static ssize_t pageIndex;
	static ssize_t pageSize;
	BinomialVertexNode() {
		leftSize = rightSize = 0;
		left = right = parent = 0;
		isPaged = 0;
	}
	void AddList(TupleVector &newList) {
		tupleList = newList;
		Compress(tupleList);
		//CompressReverseComplements(tupleList);
		totalSize = tupleList.size();
	}

	void MergeChildren() {
		assert(left != NULL);
		assert(right != NULL);
		Join(left->tupleList, right->tupleList, tupleList);
		delete left;
		delete right;
		left = right = NULL;
	}
	
	void PageOut(string pageFileName,
							 std::ostream &report = std::cout) {
		cout << "NOT DONE" << endl;
		exit(0);
		std::ofstream pageOut;
		openck(pageFileName, pageOut, std::ios::out|std::ios::binary, report);
		pageOut.write((const char *) &tupleList[0], sizeof(CountedIntegralTuple) * tupleList.size());
		//UNUSED// ssize_t listSize = tupleList.size();
		isPaged = 1;
		pageOut.close();
	}

	void MapPage(string pageFileName) {
		cout << "MapPage not done" << endl;
		exit(0);
		//		int fildes = open(pageFileName.c_str(), );
		/*
		CountedIntegralTuple *mappedList;
		mappedList = (CountedIntegralTuple*) mmap(0, 
																							sizeof(CountedIntegralTuple) * listSize , 
																							PROT_READ, MAP_PRIVATE, 
																						fildes, 0);
		*/
	}
};


ssize_t BinomialVertexNode::pageIndex = 0;
ssize_t BinomialVertexNode::pageSize  = 0;

BinomialVertexNode* AddNewList(BinomialVertexNode *&root, TupleVector &newList) {
	if (root == NULL) {
		root = new BinomialVertexNode;
		root->AddList(newList);
		return root;
	}
	else {
		if (root->left == NULL) {

			// Nodes are always added as a sibling to the current node. 
			// A new parent is created, and the current parent is added 
			// as the right child to the current node.
			// The current node is demoted to a child

			// Make the current node children.
			BinomialVertexNode *newRoot = new BinomialVertexNode, *leftChild, *rightChild;

			// the left child is new.
			leftChild = newRoot->left   = new BinomialVertexNode;
			leftChild->parent = newRoot;
			leftChild->AddList(newList);

			// the right child is the root, demoted.
			rightChild = newRoot->right = root;

		 	// Link the new root back up.
			newRoot->parent = root->parent;

			// Move the root down.
			root->parent = newRoot;
			newRoot->totalSize = root->totalSize + newRoot->left->totalSize;


			// Reassign whatever was pointing to the root to be the new root.
			root = newRoot;
			return leftChild;
		}
		else {
			// left is not null
			assert(root->left != NULL);
			assert(root->right != NULL);
			if (root->left->totalSize < root->right->totalSize) {
				return AddNewList(root->left, newList);
			}
			else {
				return AddNewList(root->right, newList);
			}
		}
	}
}

bool Close(ssize_t a, ssize_t b, double percent) {
	if ( abs(a - b) < percent*a or 
			 abs(a - b) < percent*b) 
		return 1;
	else
		return 0;
}

BinomialVertexNode *BalanceTree(BinomialVertexNode *leaf) {

	// Cannot balance the root, since if we get here, there is just 
	// one vertex.

	BinomialVertexNode *parent = leaf->parent, *root;
	if (parent == NULL)
		root = leaf;
	else
		root = parent;

	while (parent != NULL and 
				 Close(parent->left->totalSize, 
							 parent->right->totalSize, 0.15)) {
		/*		cout << "merging children of size: " << parent->left->totalSize << " " 
				 << parent->right->totalSize << endl;
		*/
		parent->MergeChildren();
		parent->totalSize = parent->tupleList.size();
		root   = parent;
		parent = parent->parent;
	}
	return root;
}

BinomialVertexNode *CollapseEntireTree(BinomialVertexNode *root) {
	if (root->left != NULL) {
		CollapseEntireTree(root->left);
		CollapseEntireTree(root->right);

		root->MergeChildren();
	}
	return root;
}

void PrintUsage() {
	cout << "usage: readsToSpectrum reads tupleSize spectrumOut" << endl;
	cout << "  -minMult M (0) Only output tuples with multiplicity >= M " << endl << endl;
	cout << "  -bufferSize B (1000000) Read tuples into buffers of size B" << endl << endl;
	cout << "  -printCount    Print the number of times each tuple appears" << endl << endl;
	cout << "  The tuple list is output in a binary format, with the tuple first and " <<endl
			 << "  the count next if printCount is given." << endl
			 << "  -usePaging  S  Page out buffers when they pass 'S' in size." << endl;
}

int main(int argc, char* argv[]) {

	//	cout << "sizeof countedintegraltuple: " << sizeof(CountedIntegralTuple) << " "
	//			 << "sizeof it: " << sizeof(IntegralTuple) << endl;

	if (argc < 4) {
		PrintUsage();
		exit(1);
	}
	string readsFileName         = argv[1];
	int tupleSize                = atoi(argv[2]);
	string spectrumOutFileName   = argv[3];
	int argi = 4;
	ssize_t minMult = 0;
	ssize_t bufferSize = 1000000;      // TODO: deal with hard-coded constants
	ssize_t printCount = 0;
	ssize_t maxReadLength = 1000000;   // TODO: deal with hard-coded constants
	while (argi < argc) {
		if (strcmp(argv[argi], "-minMult") == 0) {
			minMult = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-bufferSize") == 0) {
			bufferSize = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-printCount") == 0) {
			printCount = 1;
		}
		else if (strcmp(argv[argi], "-maxReadLength") == 0) {
			maxReadLength = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			cout << " bad option: " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}
	//	CountedIntegralTuple::tupleSize = tupleSize;
	CountedIntegralTuple::SetTupleSize(tupleSize);
	TupleVector tupleList, buffer;
	DNASequence seq;

	string reportFileName = FormReportName(readsFileName);
	std::ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);

	ifstream readsIn;
	openck(readsFileName, readsIn, std::ios::in, report);

	ssize_t r = 0;
	_SZT_ total = 0;
	BinomialVertexNode *listRoot = NULL, *listLeaf = NULL, *balancedRoot;


	std::ofstream spectrumOut;
	openck(spectrumOutFileName, spectrumOut, std::ios::out | std::ios::binary, report);
	buffer.resize(bufferSize);
	ssize_t bufCur = 0;
	while (SeqReader::GetSeq(readsIn, seq, SeqReader::noConvert)) {
		ssize_t i;
		CountedIntegralTuple t;
		i = 0;

		// advance to first valid.
		while(i < seq.length - tupleSize + 1 and
					!t.StringToTuple(&seq.seq[i])) 
			i++;

		char nuc;
		while (i < seq.length - tupleSize + 1) {
			
			if (i < seq.length - tupleSize + 1) {
				buffer[bufCur] = t;
				++bufCur;
				//				buffer.push_back(t);
			}
			
			nuc = seq.seq[i+tupleSize];
			if (nuc != 'A' and
					nuc != 'C' and
					nuc != 'T' and
					nuc != 'G') {
				i+= tupleSize;
				while (i < seq.length - tupleSize + 1 and
							 !t.StringToTuple(&seq.seq[i]))
					i++;
			}
			else {
				if (i+tupleSize < maxReadLength) {
					CountedIntegralTuple next;
					ForwardNuc(t, numeric_nuc_index[seq.seq[i+tupleSize]], next);
					t = next;
				}
			}
			i++;
			total++;
			if (bufCur == bufferSize) {
				listLeaf = AddNewList(listRoot, buffer);
				balancedRoot = BalanceTree(listLeaf);
				bufCur = 0;
			}
		}
				
		/*
		for (i = 0; i < seq.length - tupleSize + 1; i++ ){ 
			if (t.StringToTuple(&seq.seq[i]))
				tupleList.push_back(t);
		}
		*/
		r++;
		if (r % 500000 == 0) {
			cout << r << endl;
		}
	}
	// only one buffer has been created.
	buffer.resize(bufCur);
	listLeaf = AddNewList(listRoot, buffer);
	
	balancedRoot = CollapseEntireTree(listRoot);

	std::sort(listRoot->tupleList.begin(), listRoot->tupleList.end(), CompareCountedIntegralTuples());
	
	CompressReverseComplements(listRoot->tupleList);
	
	_SZT_ i;
	// _SZT_ nTuples = listRoot->tupleList.size();


	// First count the number above threshold.
	ssize_t freqListSize = 0;
	ssize_t numPalindrome = 0;
	CountedIntegralTuple rc;
	for (i = 0; i < listRoot->tupleList.size(); i++) {
		if (listRoot->tupleList[i].count >= minMult) {
			freqListSize++;
			listRoot->tupleList[i].MakeRC(rc);
			//			if (rc.tuple == listRoot->tupleList[i].tuple)
			if (rc == listRoot->tupleList[i]) // compare tuples
				numPalindrome++;
		}
	}	
	freqListSize *=2;
	freqListSize -= numPalindrome;

	ssize_t tupleSize_SSZT = tupleSize;
	spectrumOut.write((const char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));

	spectrumOut.write((const char*) &freqListSize, sizeof(ssize_t));
	for (i = 0; i < listRoot->tupleList.size(); i++) {
		if (listRoot->tupleList[i].count >= minMult) {
			spectrumOut.write((const char*) &listRoot->tupleList[i], sizeof(CountedIntegralTuple));
			listRoot->tupleList[i].MakeRC(rc);
			//			if (listRoot->tupleList[i].tuple != rc.tuple) {
			if (listRoot->tupleList[i] != rc) { // compare tuples
				rc.count = listRoot->tupleList[i].count;
				//UNUSED// LongTuple tuple = rc.tuple;
					spectrumOut.write((const char*) &rc, sizeof(CountedIntegralTuple));		
			}
		}
	}
	spectrumOut.close();
	EndReport(report);
	report.close();
	return 0;
}
