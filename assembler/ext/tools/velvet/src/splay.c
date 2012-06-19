/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#include <stdlib.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "globals.h"
#include "recycleBin.h"
#include "kmer.h"
#include "utility.h"

#define CHUNKSIZE 10000

static RecycleBin *treeMemory = NULL;

struct splayNode_st {
	Kmer kmer;
	Coordinate position;
	struct splayNode_st *left;
	struct splayNode_st *right;
	IDnum seqID;
}  ATTRIBUTE_PACKED;

typedef struct splayNode_st SplayNode;
typedef struct splayNode_st SplayTree;


#ifdef _OPENMP
void initSplayTreeMemory(void)
{
	int n;

	n = omp_get_max_threads();
	#pragma omp critical
	if (treeMemory == NULL)
		treeMemory = newRecycleBinArray(n, sizeof(SplayNode), CHUNKSIZE);
}
#endif

static SplayNode *allocateSplayNode()
{
#ifdef _OPENMP
#ifdef DEBUG
	if (treeMemory == NULL)
	{
		velvetLog("The memory for splay trees seems uninitialised, "
			  "this is probably a bug, aborting.\n");
		abort();
	}
#endif
	return allocatePointer(getRecycleBinInArray(treeMemory,
						    omp_get_thread_num()));
#else
	if (treeMemory == NULL)
		treeMemory = newRecycleBin(sizeof(SplayNode), CHUNKSIZE);

	return allocatePointer(treeMemory);
#endif
}

void destroyAllSplayTrees()
{
#ifdef _OPENMP
	destroyRecycleBinArray(treeMemory);
#else
	destroyRecycleBin(treeMemory);
#endif
	treeMemory = NULL;
}

/* This function can be called only if K2 has a left child */
/* Perform a rotate between a node (K2) and its left child */
/* Update heights, then return new root */

static SplayNode *SingleRotateWithLeft(SplayNode * K2)
{
	SplayNode *K1;

	K1 = K2->left;
	K2->left = K1->right;
	K1->right = K2;

	return K1;		/* New root */
}

/* This function can be called only if K1 has a right child */
/* Perform a rotate between a node (K1) and its right child */
/* Update heights, then return new root */

static SplayNode *SingleRotateWithRight(SplayNode * K1)
{
	SplayNode *K2;

	K2 = K1->right;
	K1->right = K2->left;
	K2->left = K1;

	return K2;		/* New root */
}

/* Top-down splay procedure, */
/* not requiring kmer to be in tree */

static SplayTree *Splay(Kmer * kmer, SplayTree * T)
{
	SplayNode Header;
	SplayNode *LeftTreeMax, *RightTreeMin;

	if (T == NULL)
		return NULL;

	Header.left = Header.right = NULL;
	LeftTreeMax = RightTreeMin = &Header;

	while (compareKmers(kmer, &(T->kmer))) {
		if (compareKmers(kmer, &(T->kmer)) < 0) {
			if (T->left == NULL)
				break;
			if (compareKmers(kmer, &(T->left->kmer)) < 0)
				T = SingleRotateWithLeft(T);
			if (T->left == NULL)
				break;
			/* Link right */
			RightTreeMin->left = T;
			RightTreeMin = T;
			T = T->left;
		} else {
			if (T->right == NULL)
				break;
			if (compareKmers(kmer, &(T->right->kmer)) > 0)
				T = SingleRotateWithRight(T);
			if (T->right == NULL)
				break;
			/* Link left */
			LeftTreeMax->right = T;
			LeftTreeMax = T;
			T = T->right;
		}
	}			/* while kmer != T->kmer */

	/* Reassemble */
	LeftTreeMax->right = T->left;
	RightTreeMin->left = T->right;
	T->left = Header.right;
	T->right = Header.left;

	return T;
}

Kmer * findInTree(Kmer * X, SplayTree ** T)
{
	*T = Splay(X, *T);
	return &((*T)->kmer);
}

void insertIntoTree(Kmer * kmer, SplayTree ** T)
{
	SplayNode *newNode;

	if (*T == NULL) {
		newNode = allocateSplayNode();
		copyKmers(&(newNode->kmer), kmer);
		newNode->left = newNode->right = NULL;
		*T = newNode;
		return;
	}

	*T = Splay(kmer, *T);
	if (compareKmers(kmer, &((*T)->kmer)) < 0) {
		newNode = allocateSplayNode();
		copyKmers(&(newNode->kmer), kmer);
		newNode->left = (*T)->left;
		newNode->right = *T;
		(*T)->left = NULL;
		*T = newNode;
	} else if (compareKmers(&((*T)->kmer), kmer) < 0) {
		newNode = allocateSplayNode();
		copyKmers(&(newNode->kmer), kmer);
		newNode->right = (*T)->right;
		newNode->left = *T;
		(*T)->right = NULL;
		*T = newNode;
	}
}

boolean
findOrInsertOccurenceInSplayTree(Kmer * kmer, IDnum * seqID,
				 Coordinate * position, SplayTree ** T)
{
	SplayNode *newNode;

	if (*T == NULL) {
		newNode = allocateSplayNode();
		copyKmers(&(newNode->kmer), kmer);
		newNode->seqID = *seqID;
		newNode->position = *position;

		newNode->left = newNode->right = NULL;

		*T = newNode;

		return false;
	}

	*T = Splay(kmer, *T);
	if (compareKmers(kmer, &((*T)->kmer)) < 0) {
		newNode = allocateSplayNode();
		copyKmers(&(newNode->kmer), kmer);
		newNode->seqID = *seqID;
		newNode->position = *position;

		newNode->left = (*T)->left;
		newNode->right = *T;
		(*T)->left = NULL;

		*T = newNode;

		return false;
	} else if (compareKmers(kmer, &((*T)->kmer)) > 0) {
		newNode = allocateSplayNode();
		copyKmers(&(newNode->kmer), kmer);
		newNode->seqID = *seqID;
		newNode->position = *position;

		newNode->right = (*T)->right;
		newNode->left = *T;
		(*T)->right = NULL;

		*T = newNode;

		return false;
	} else {
		*seqID = (*T)->seqID;
		*position = (*T)->position;

		return true;
	}
}
