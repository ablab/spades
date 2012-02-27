/***************************************************************************
 * Title:          TupalLib.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _TUPLE_LIB
#define _TUPLE_LIB
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <iomanip>
#include <algorithm>
#include "pvm3.h"
#include "Params.h"
#include "SeqReader.h"
#include "TupleLib.h"
#include "DNASequence.h"
#include "CommunicationLib.h"
#include "StripGen.h"
#include "alignutils.h"
#include "MaximizeStrip.h"
#include "Alignment.h"
#include "Mapping.h"

class T_FragmentLocation {
public:
  ssize_t contig;
  ssize_t start;
  ssize_t end;
  ssize_t length; // not always start - end+1
  ssize_t strand;
  ssize_t contigsSpanned;
  // deep copy operator
  T_FragmentLocation &operator=(const T_FragmentLocation &c) {
    if (this != &c) {
      contig  = c.contig;
      start   = c.start;         
      end     = c.end;           
      length  = c.length;        
      strand  = c.strand;        
      contigsSpanned = c.contigsSpanned;
		}
    return *this;
  }
  // copy constructor
  T_FragmentLocation(const T_FragmentLocation &c) {
    (*this) = c;
  }
  // default constructor
  T_FragmentLocation() {
    contig   = 0;
    start    = 0;
    end		   = 0;
    length	 = 0;
    strand	 = 0;
    contigsSpanned = 0;
  }
};

// Define types
typedef std::vector<DNASequence*> T_Sequences;
typedef std::vector<T_Sequences*> T_GenomeVect;
typedef std::vector<T_FragmentLocation> T_Locations;
typedef std::vector<ssize_t*> T_IntegralValueVect;

// Functions

void ReadParamFile(std::string paramFileName, 
		   T_Params &parameters, 
		   T_ParallelParams &parallelParams);

// Take a sequence and split it into numFragments fragments.
// Record where in the genomne all of the fragments were taken from.
void FragmentOneGenome(DNASequence &sequence, 
		       ssize_t numFragments, ssize_t hashLength,
		       T_Sequences &fragments);

// Given a bunch of fragments that were created from sequence, combine
// them back into sequence.

void UnfragmentGenome( T_Sequences &fragments, 
		       ssize_t hashLength,
		       DNASequence &sequence );

void CombineEnumerations( T_IntegralValueVect
			  &fragmentEnumerations, 
			  T_Sequences &fragments,
			  ssize_t wordLength,
			  ssize_t *enumerations, ssize_t totalLength );

void PrintEnumeration( DNASequence &seq, 
		       ssize_t hashLength,
		       ssize_t *enumerations,
		       ssize_t *locations,
		       std::ostream &out );



void CombineVectors( T_IntegralValueVect
		     &fragmentVectors, 
		     T_Sequences &fragments,
		     ssize_t wordLength,
		     ssize_t *vector, ssize_t totalLength );

void JoinStrings( std::vector<char*> &strings, 
		  std::vector<ssize_t> &lengths, 
		  char *&seq,
		  ssize_t  &len );

void DistributeMasking(T_Sequences &fragments, 
		       T_Params &params,
		       T_ParallelParams &parInfo);

void DistributeUnmaskShared(T_GenomeVect &genomeFragments, 
			    T_Params &params,
			    T_ParallelParams &parInfo);

void DistributeEnumerateGenome(T_GenomeVect &genomeFragments, 
			       T_IntegralValueVect &permutedEnumerations,
			       T_IntegralValueVect &permutedLocations,
			       ssize_t hashLength, ssize_t wordLength,
			       ssize_t myid,
			       ssize_t numChildren, ssize_t *childIds);


void FragmentGenomes(T_GenomeVect &genomes,
		     T_GenomeVect &genomeFragments, 
		     ssize_t numFragments,
		     T_Params &params);

void MaskFragments(T_GenomeVect &genomeFragments, 
		   T_Params params);

void MergeVectors(ssize_t *vect1, ssize_t *vect2, ssize_t len, ssize_t unsetValue=0);

void DistributeMaskingGenomes(T_GenomeVect genomeFragments, 
			      T_Params &mainParams, 
			      T_ParallelParams &parInfo);

void InitializeNotCommon(T_GenomeVect &genomeFragments, T_Params params);

void UnmaskSharedGenomes(T_GenomeVect &genomes, T_Params &params);

ssize_t EnumerateSequences(DNASequence &refSeq, DNASequence &querySeq, 
		       T_Params &mainParams,
		       T_ParallelParams &parInfo,
		       ssize_t *&enumerations,
		       ssize_t *&locations);

ssize_t EnumerateSequences(DNASequence &refSeq, DNASequence &querySeq, 
		       T_Params &mainParams,
		       T_ParallelParams &parInfo,
		       T_Alignment *&alignment);

void CreateStripArrays(DNASequence &refSeq,  // description of first seq
		       // permutation of the second sequence
		       T_Params params,
		       ssize_t *& seqEnumerations, ssize_t *& seqLocations, 
		       ssize_t & numTuples,
		       // Output
                       //    condensed enumerations
		       ssize_t *& enumerations,  
		       //    Condensed locations
		       ssize_t *& refLocations, ssize_t *&qryLocations
		       );


ssize_t FindGap(DNASequence &refSeq, 
	     ssize_t *enumerations, ssize_t *locations, 
	     ssize_t minDist, ssize_t maxDist, 
	     ssize_t &refPos, 
	     ssize_t &refStartPos, ssize_t &refEndPos );


void BuildScaffold(DNASequence &refSeq, DNASequence &qrySeq,
		   T_Alignment *&alignment,
		   T_Strips *&strips,
		   T_Params params, T_ParallelParams &parInfo,
		   std::string padding);

void EnumerateRecursively(DNASequence &refSeq, DNASequence &qrySeq,
			  T_Alignment *&alignment,
			  ssize_t dir,
			  T_Params params, T_ParallelParams &parInfo,
			  std::string padding);

T_Alignment *AlignSubsequences(DNASequence &refSeq, DNASequence &qrySeq,
			       ssize_t startRefPos, ssize_t endRefPos,
			       ssize_t startQryPos, ssize_t endQryPos,
			       T_Alignment *alignment,
			       T_Params params, T_ParallelParams parInfo, 
			       std::string padding);

void AssignAlignmentLocations(DNASequence &seq,
			      T_Alignment *alignment, 
			      ssize_t refOffsset, 
			      ssize_t qryOffset,
			      ssize_t *locations,
			      ssize_t *enumerations,
			      ssize_t lenLocations, 
			      ssize_t *qryLocations, ssize_t lenLocations, ssize_t level);



void LocationsToEnumeration(ssize_t *locations, ssize_t *enumerations, ssize_t length, 
			    ssize_t *enumeration, ssize_t enumLength);

void StoreEnumerationInAlignment(ssize_t *locations,
				 ssize_t *enumerations, //T_Alignment &alignment,
				 ssize_t length,
				 ssize_t *qryEnumerations,
				 ssize_t enumLength);
void CopyAlignment(DNASequence &seq,
		   T_Alignment *alignment, ssize_t parentRefOffset, ssize_t parentQryOffset, 
		   ssize_t *locations, ssize_t *enumerations, ssize_t lenLocations,
		   ssize_t *qryLocations, ssize_t qryLength);

void FreeAlignments(T_Alignment *alignment);

void AffineAlign2Seq(DNASequence &refSeq, DNASequence &qrySeq, 
		     T_Alignment*&alignment,
		     T_Params params, std::string padding,
		     ssize_t dir);

double EndAlign2Seq(DNASequence &longSeq, DNASequence &shortSeq,
		   T_Alignment *&alignment, ssize_t dir);

void ReannotateAlignment(DNASequence &seq,
			 T_Strips *strips,
			 T_Alignment *alignment);

ssize_t ExtendMatch(DNASequence &refSeq, DNASequence &qrySeq, 
		ssize_t refPos, ssize_t qryPos, 
		ssize_t dir, double perror);

void AlignStrips(DNASequence &refSeq,
			 DNASequence &qrySeq,
			 T_Alignment *alignment, 
			 T_Alignment *&curAlignment,
			 T_Strips *&strips,
			 T_Params params,
			 T_ParallelParams &parInfo,
			 std::string padding);

void AlignBetweenStrips(DNASequence  &refSeq,
			DNASequence  &qrySeq,
			T_Alignment *&alignment,
			T_Alignment *&curAlignment,
			T_Strips    *&strips,
			T_Params        params,
			T_ParallelParams &parInfo,
			std::string padding);


#endif
