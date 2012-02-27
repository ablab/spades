/***************************************************************************
 * Title:          TupleLib.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _TUPLELIB_H_
#define _TUPLELIB_H_
#include <iostream>
#include <string>
#include <fstream>
#include <ostream>

#include <unistd.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "Alignment.h"
#include "../tupal/Params.h"
#include "TupleDef.h"
#include "modify_tuple.h"
#include "DNASequence.h"
 

#define FORWARD 1
#define REVERSE -1

/*
In order.
244, 245, 246, 247, 248, 249
g    a    c    t    g    a 
250, 251, 252, 253, 254, 255}
c    t    g    a    c    t
*/

class Node;
class Node {
public:
  Tuple value;     // Value of the tuple for this node.
  static  long total_nodes;
  static  ssize_t totalLookups;
  static  ssize_t totalReferences;
  Node *match;      // begin by encoding match as !NULL for unique or
                    // NULL not.  Next if the tuple is shared, match
                    // is set to the shared tuple.
  Node *next;       // next node in chained hash table
  ssize_t location;     // Location of this tuple in the genome. Location
                    // also encodes the direction.  Positive location
                    // encodes forward, negative reverse.
  ssize_t enumeration;  // Order of appearance of this node in the list of
                    // unique nodes.

  Node(Tuple value = -1) {
    value = -1;
    total_nodes++;
    enumeration = 0;
    match = this;
  }
  Node(Tuple v, long loc, Node* n) {
    value = v;
    next  = n;
    match = this;  // set to null when this node is not unique, or
                   // its partner is not unique.   
    enumeration = 0;
    location  = loc;
    total_nodes++;
  }
  void setEnumeration(ssize_t e) {
    enumeration = e;
  }
  ssize_t getEnumeration() {
    return enumeration;
  }
  void setLocation(ssize_t l) {
    location = l;
  }
  ssize_t getLocation() {
    return location;
  }
  ssize_t IsUnique() {
    return match == this;
  }
};



ssize_t  GetNextValidPos(DNASequence &seq, ssize_t pos, ssize_t tupleLength, ssize_t masked);
void MaskNotUnique(DNASequence &seq1, DNASequence &seq2,
		   ssize_t dist, ssize_t doIndel,
		   ssize_t hashLength, ssize_t wordLength);
void UnmaskShared(DNASequence &seq1,
		   DNASequence &seq2,
		  ssize_t hashLength, ssize_t wordLength);

void EnumerateUnique(DNASequence &seq1, 
		     DNASequence &seq2, 
		     ssize_t hashLength, ssize_t wordLength,
		     ssize_t *permutedEnumerations,
		     ssize_t *permutedLocations);

Tuple trans_seq(DNASequence &seq, ssize_t pos, ssize_t len);


void StoreGenome(DNASequence &seq,long &startPos, Node* table[], 
		 ssize_t probelen, ssize_t wordlen, ssize_t storerc = 1, ssize_t enumeration = 0);

std::ostream& printTuple(Tuple tup, ssize_t tuple_size, std::ostream &out);

void printGenome(char *s, long len, Node *reference[], Node*
		 permuted[], std::ofstream &out, ssize_t probe_len, ssize_t word_len, 
		 ssize_t outputTuple);

void findCommon(Node *tableA[], Node *tableB[], long tableSize, int
		probe_len, ssize_t word_len);

void store(Node* table[], Tuple val, long pos, char nextNuc, ssize_t word_len,
	   ssize_t enumeration = 0);
void AllocateTable(Node **&table, ssize_t wordLength);
void DeallocateTable(Node **&table, ssize_t wordLength);
void findNode(Node* table[],  Node *&node, Tuple tuple, ssize_t probe_len, ssize_t word_len);

Tuple trans_seq(char *seq, ssize_t len);
Tuple getHashValue(Tuple val, ssize_t probe_len, ssize_t word_len);

long assignEnumerations(Node* head, Node *tableA[], Node *tableB[], ssize_t probe_len, ssize_t word_len, ssize_t nt, ssize_t doIndel, Node*&newHead);
typedef Node* nodeptr;
void PrintUniqueCommon(Node* ahead,  Node* tableB[], ssize_t probe_len, ssize_t word_len, std::ofstream &out);
void FilterUniqueCommon(Node* tableA[], Node* tableB[], long
			table_size, ssize_t probe_len, ssize_t word_len,
			std::ofstream &out);

void search_k(Node* result, ssize_t & numNeighbors, Node* ref[], Tuple tuple, 
	      ssize_t k, ssize_t probe_len, ssize_t word_len, Tuple &match, ssize_t start, char next, ssize_t doIndel);
ssize_t search_to_k(Node* result, ssize_t & numNeighbors, Node* ref[], Tuple tuple, 
		ssize_t k, ssize_t probe_len, ssize_t word_len, Tuple &match, char next, ssize_t doIndel);


void enumerateGenome(char *s,  long length, long &enumeration, long &location, 
		     ssize_t dir, Node *tableA[], Node* tableB[], ssize_t probe_len,
		     ssize_t word_len,
		     ssize_t k, ssize_t doIndel);



void PrintLocations(ssize_t *locations, ssize_t length, std::ostream &out);


template<typename t>
void BinaryPrintArray(t *data, ssize_t &length, std::ostream &out) 
{
	// TODO: may need to update file format to have tupleSize at start, as .spect does
  out.write((char*) &length, (std::streamsize) (sizeof(ssize_t)));
  out.write((char*) data, (std::streamsize) (sizeof(t)*length));
}

template<typename T>
void BinaryReadArray(T *&data, ssize_t &length, std::istream &in){
	// TODO: may need to update file format to have tupleSize at start, as .spect does
  in.read((char*) &length, (std::streamsize) (sizeof(ssize_t)));
  data = new T[length];
  in.read((char*) data, (std::streamsize) (sizeof(T)*length));
}

void PrintCheckpoint(T_Alignment &alignment, T_Params &params, 
		     DNASequence &refSeq, DNASequence &qrySeq, 
		     std::string outFileName);

void ReadCheckpoint(T_Alignment &alignment, T_Params &params, 
		    DNASequence &refSeq, DNASequence &qrySeq,
		    std::string inFilename);

void RotateAlignment(ssize_t *alignment, ssize_t alignLength, ssize_t alignSpan);


#endif
