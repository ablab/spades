/***************************************************************************
 * Title:          StripGen.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _STRIPGEN_H_
#define _STRIPGEN_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include "EnumUtils.h"

using namespace std;

class Strip {
public:
  long startQryEnum;
  long endQryEnum;
  long startRefEnum;
  long endRefEnum;
  long startRefPos;
  long endRefPos;
  long startQryPos;
  long endQryPos;
  ssize_t numMerged;
  ssize_t hashLength;
  ssize_t size;
  ssize_t sign;
  Strip() {
    numMerged = 0;
    size =0;
  }
};


typedef std::list<Strip> T_Strips;

class SortStripsByQryEnum  {
public:
  bool operator()(T_Strips::iterator a, T_Strips::iterator b) {
    return labs(a->startQryEnum) > labs(b->startQryEnum);
  }
};

class SortStripsByRefEnum {
public:
  bool operator()(T_Strips::iterator a, T_Strips::iterator b) {
    return labs(a->startRefEnum) > labs(b->startRefEnum);
  }
};

class SortStripsByQryPos {
public:
  bool operator()(T_Strips::iterator a, T_Strips::iterator b) {
    return labs(a->startQryPos) > labs(b->startQryPos);
  }
};

typedef  std::priority_queue<T_Strips::iterator, std::vector<T_Strips::iterator>, SortStripsByQryEnum > T_StripQQry;
typedef  std::priority_queue<T_Strips::iterator, std::vector<T_Strips::iterator>, SortStripsByRefEnum > T_StripQRef;

typedef std::priority_queue<T_Strips::iterator, std::vector<T_Strips::iterator>, SortStripsByQryPos > T_StripQQryPos;

class  TreeStrip: public Strip {
public:
  TreeStrip *leftChild;
  TreeStrip *rightChild;
  long enumValue;
  TreeStrip() {
    leftChild = NULL;
    rightChild = NULL;
    enumValue = 0;
  }
  friend bool  operator< ( const  TreeStrip &a, const TreeStrip &b) const {
    return a.startRefEnum > b.startRefEnum;
  }
  bool EnumerationAssigned() { return enumValue != 0; }
};

struct SortTreeStrips {
  public:
  bool operator()(TreeStrip*a, TreeStrip*b) {
    return a->startRefEnum > b->startRefEnum;
  }
};

typedef std::queue<TreeStrip*> Squeue;

void Insert(TreeStrip *&head, Strip newStrip);

TreeStrip * BuildStripRelations(T_Strips strips);


/* 
   Input: 
     strips - a list of strips in order of appearance in a reference
   genome, and enumerations in the query genome.  The enumerations are
   the original enumerations before merging strips or removing small
   strips.  There are n = |strips| strips as input.
   
   Result:
     The strips are enumerated e_i, i = 1..n, e_i \in [i,-1].

*/
void EnumerateStrips(T_Strips &strips);

void PrintEnumerations(T_Strips &strips, std::ostream & out);

ssize_t MergeStrips(T_Strips::iterator prev, T_Strips::iterator cur, ssize_t forceMerge);

ssize_t DiscardSmallStrips(T_Strips& strips, ssize_t distThresh);

ssize_t DiscardAndMerge(T_Strips &strips, ssize_t distThreshold);

void Merge(T_Strips &strips);

void BFSPrintEnumerations(TreeStrip *head, std::ofstream &out);

void AddStripsToPQueue(TreeStrip *head, 
		       std::priority_queue<TreeStrip> &pqueue);

void PQPrintEnumerations(TreeStrip *head, std::ofstream & out);

void DFSAssignEnumerations(TreeStrip *head);

ssize_t GetLineOrig(std::ifstream &in,  long &curEnumRef, long &curEnumComp, long &curPosRef, long &curPosComp, ssize_t withSeq);

ssize_t GetLineSigned(std::ifstream &in, long &enumRef, long &enumComp, 
		long &startPosRef, long &endPosRef, 
		  long &startPosComp, long &endPosComp, long & size, ssize_t& numMerged, ssize_t &sign);

ssize_t GetLineUnsigned(std::ifstream &in, long &enumRef, long &enumComp, 
		    long &startPosRef, long &endPosRef, 
		    long &startPosComp, long &endPosComp, long & size, ssize_t& numMerged);

ssize_t GetLineSecond(std::ifstream &in, long &enumRef, long &enumComp, 
		  long &startPosRef, long &endPosRef, 
		  long &startPosComp, long &endPosComp, long & size, ssize_t& numMerged);

ssize_t GetStrips(std::ifstream &inFile, T_Strips &strips, ssize_t &numRead, ssize_t mergeAdjacent = 1);

ssize_t GetStrips(ssize_t *enumerations, ssize_t *startRefLocations, ssize_t *endRefLocations,
	       ssize_t *startQryLocation, ssize_t *endQryLocation, ssize_t *stripSize, 
	      ssize_t *stripsMerged, ssize_t length, T_Strips &strips, ssize_t hashLength, ssize_t mergeAdjacent = 1);

void ReadStrips(std::string inName, 
		ssize_t *&enumerations, ssize_t *&startRefLocations, ssize_t *&endRefLocations,
		ssize_t *&startQryLocation, ssize_t *&endQryLocation, ssize_t *&stripSize, 
		ssize_t *&stripsMerged, ssize_t &length);

void ReadStrips(std::ifstream &inFile, 
		ssize_t *&enumerations, ssize_t *&startRefLocations, ssize_t *&endRefLocations,
		ssize_t *&startQryLocation, ssize_t *&endQryLocation, ssize_t *&stripSize, 
		ssize_t *&stripsMerged, ssize_t &length);

void CreateStrip(ssize_t stripStart, ssize_t pos, ssize_t enumRefStart, ssize_t enumRefEnd,
		 ssize_t *enumerations, ssize_t *startRefLocations, ssize_t *endRefLocations, 
		 ssize_t *startQryLocations, ssize_t *endQryLocations, ssize_t stripSize, ssize_t numMerged, 
		 Strip &strip, ssize_t hashLength);

ssize_t GetMUMs(ssize_t *enumerations, ssize_t *startRefLocations, ssize_t *endRefLocations,
	    ssize_t *startQryLocations, ssize_t *endQryLocations, ssize_t *stripSizes, 
	    ssize_t *stripsMerged, ssize_t length, T_Strips &strips, ssize_t hashlength);

ssize_t GetStripsMerged(ssize_t *merged, ssize_t pos);
ssize_t GetStripSize(ssize_t *strip, ssize_t pos);

template <typename T>
void BuildStripPqueue(T_Strips &strips, T &queue) {
  // BuildStripPqueue is used for ordering the strips according to
  // their output order.  They are generated in order of appearance in
  // the input file.
  long l;
  T_Strips::iterator it, end = strips.end();
  for (it = strips.begin(); it != end; ++it)
    queue.push(it);
}

template <typename T>
void AssignEnumerations(T &queue) {
  long e = 1;
  T_Strips::iterator it;
  while (queue.size() > 0) {
    it = queue.top();
    queue.pop();
    it->startQryEnum = e++ * it->sign;
    it->endQryEnum = it->startQryEnum;
  }
}

#endif
