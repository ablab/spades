/***************************************************************************
 * Title:          ReadsToVector.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <vector>
#include <string>
using namespace std;

#include "SeqReader.h"
#include "DNASequence.h"
#include "utils.h"
#include "IntegralTupleStatic.h"

typedef std::vector<CountedIntegralTuple> TupleVector;

void CompressReverseComplements(TupleVector &tupleList) {
	std::sort(tupleList.begin(), tupleList.end());
	ssize_t i;
	CountedIntegralTuple rc;
	i = 0;

	// Combine rev comp multiplicities.
	for (i = 0; i < tupleList.size(); i++ ) {
		if (tupleList[i].GetMult() < 0)
			continue;
		tupleList[i].MakeRC(rc);
		if (rc > tupleList[i]) { // compare tuples
		//		if (rc.tuple > tupleList[i].tuple) {
			ssize_t rcTupleIndex;
			rcTupleIndex = LookupBinaryTuple(&tupleList[i+1], tupleList.size() - i - 1, rc);
			if (rcTupleIndex != -1) {
				tupleList[i].SetMult(tupleList[i].GetMult() + tupleList[rcTupleIndex].GetMult());
				tupleList[rcTupleIndex].SetMult(-1);
			}
		}
	}
	// Remove the rc from the list.
	ssize_t cur = 0;
	for (i = 0, cur = 0; i < tupleList.size(); i++) {
		if (tupleList[i].GetMult() >= 0) {
			tupleList[cur].SetTuple(tupleList[i]);
			//			tupleList[cur].tuple = tupleList[i].tuple;
			tupleList[cur].SetMult(tupleList[i].GetMult());
			cur++;
		}
	}
	cout << "removed " << tupleList.size() - cur << " reverse complements of " << tupleList.size() << endl;
	tupleList.resize(cur);
}

void Compress(TupleVector &tupleList) {
	std::sort(tupleList.begin(), tupleList.end());
	ssize_t i, c, n;
	i = 0; c = 0;
	ssize_t numRemoved = 0;
	ssize_t endI = tupleList.size() - 1;
	ssize_t prevI;
	ssize_t mult;
	for (i = 0; i < endI ; ){ 
		tupleList[c] = tupleList[i];
		mult = 1;
		prevI = i;
		while (tupleList[c] == tupleList[i]) {
			i++;
			mult++;
		}
		tupleList[c].SetMult(mult);
		c++;
		numRemoved += (i - prevI - 1);
	}
	//	cout << "Removed: " << numRemoved << " during compression." << endl;
	tupleList.resize(c);
}


void Join(TupleVector &srcA, 
					TupleVector &srcB, 
					TupleVector &dest) {
	ssize_t a = 0, b = 0;

	ssize_t nA = srcA.size(), nB = srcB.size();
	ssize_t d;
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
			srcA[a].SetMult(srcA[a].GetMult() + srcB[b].GetMult());
			dest.push_back(srcA[a]);
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
	CompressReverseComplements(dest);
}

void Update(TupleVector &src,
						TupleVector &dest) {
	ssize_t nSrc = src.size(), nDest = dest.size();
	
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
	BinomialVertexNode() {
		leftSize = rightSize = 0;
		left = right = parent = 0;
	}
	void AddList(TupleVector &newList) {
		tupleList = newList;
		Compress(tupleList);
		CompressReverseComplements(tupleList);
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

	
};

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

	BinomialVertexNode *parent = leaf->parent;

	while (parent != NULL and 
				 Close(parent->left->totalSize, 
							 parent->right->totalSize, 0.15)) {
		/*		cout << "merging children of size: " << parent->left->totalSize << " " 
				 << parent->right->totalSize << endl;
		*/
		parent->MergeChildren();
		parent->totalSize = parent->tupleList.size();
		parent = parent->parent;
	}
}


int main(int argc, char* argv[]) {

	if (argc < 3) {
		cout << "usage: readsToSpectrum reads tupleSize spectrum" << endl;
		exit(1);
	}
	string readsFileName = argv[1];
	int tupleSize =atoi(argv[2]);
	CountedIntegralTuple::SetTupleSize(tupleSize);
	//	CountedIntegralTuple::tupleSize = tupleSize;

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
	BinomialVertexNode *listRoot = NULL, *listLeaf = NULL;

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
				buffer.push_back(t);
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
				CountedIntegralTuple next;
				ForwardNuc(t.tuple, numeric_nuc_index[seq.seq[i+tupleSize]], next.tuple);
				t = next;
			}
			i++;
			total++;
			if (total %10000000 == 0) {
				listLeaf = AddNewList(listRoot, buffer);
				BalanceTree(listLeaf);
				buffer.clear();
			}
		}
				
		/*
		for (i = 0; i < seq.length - tupleSize + 1; i++ ){ 
			if (t.StringToTuple(&seq.seq[i]))
				tupleList.push_back(t);
		}
		*/
		r++;
		if (r % 50000 == 0) {
			cout << r << endl;
		}
	}

	std::sort(tupleList.begin(), tupleList.end());

	EndReport(report);
	report.close();
	return 0;
}
