/***************************************************************************
 * Title:          TupalLib.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "TupalLib.h"
#include <cmath>
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

void PrintEnumeration( DNASequence &seq, 
		       ssize_t hashLength,
		       ssize_t *enumerations,
		       ssize_t *locations,
		       std::ostream &out ) {
  ssize_t i;
  ssize_t refEnumeration;
  refEnumeration = 0;
  ssize_t numTuples  = 0;
  // Count the number of positions.
  for (i = 0; i < seq.length - hashLength; i++)
    if (seq.IsACTG(i) && !seq.IsPosMasked(i))
      numTuples++;

  out << numTuples << std::endl;
  for (i = 0; i < seq.length - hashLength + 1; i++) {
    if ((!seq.IsACTG(i) && locations[i] != -1) || !seq.IsPosMasked(i)) {
      if (seq.IsACTG(i)) 
	assert(enumerations[i] != 0 && locations[i] != -1);
      //      refEnumeration++;
      out << enumerations[i] << "\t" << i << "\t" << i << "\t"
	  << locations[i] << "\t" << locations[i] << "\t1\t0" 
	  << std::endl; 
    }
    else {
      assert(locations[i] == -1);
      assert(enumerations[i] == 0);
    }
  }
}

void FragmentOneGenome(DNASequence &genome, 
		       ssize_t numFragments, ssize_t hashLength,
		       T_Sequences &fragments) {
  // Break a string up into numFragments fragments, and remember where
  // they came from.
  ssize_t fragmentLength;
  ssize_t fullFragmentLength;
  fragmentLength = (genome.length - hashLength+1) / numFragments;
  ssize_t curPos = 0;
  ssize_t useForward = 1;
  ssize_t i;
  DNASequence *fragment;
  // determine how long this fragment will be
  fullFragmentLength = fragmentLength + hashLength-1;
  for (i = 0; i < numFragments-1; i++) {
    fragment = new DNASequence(fullFragmentLength);
    // Copy the fragment including the last tuple
    memcpy(fragment->seq, &(genome.seq[curPos]), fullFragmentLength);
    fragment->length = fullFragmentLength;
    fragment->startPosition = curPos;

    // what do we have so far?
    curPos += fragmentLength;
    fragments.push_back(fragment);
  }
	//  fragmentLength = (int) ceil(double(genome.length - hashLength+1) / numFragments);
  fragmentLength = (genome.length - hashLength + numFragments) / numFragments;
  fullFragmentLength = fragmentLength + hashLength - 1;
  // Save some space (especially when numchildren = 1)  by not copying
  // the last fragment. 
  fragment = new DNASequence;
  fragment->seq = &genome.seq[curPos];
  fragment->length = genome.length - curPos;
  fragment->startPosition = curPos;
  /*  fragment->PrintSeq(std::cout);*/
  fragments.push_back(fragment);
}


void UnfragmentGenome(T_Sequences &fragments, 
		      ssize_t hashLength,
		      DNASequence &sequence) {
  // Combine all fragments into sequence

  // assume the last fragments is stored in sequence, and so it does
  // need to be copied.
  
  ssize_t i, j, pos;
  DNASequence *fragment;
  pos = 0;
  for (i = 0; i < (fragments.size() - 1); i++) {
    fragment = fragments[i];
    for (j = 0; j < (fragment->length - hashLength+1); j++) {
      sequence.seq[pos] = fragment->seq[j];
      pos++;
    }
  }
}

void CombineVectors(T_IntegralValueVect &fragmentVectors, 
		    T_Sequences &fragments,
		    ssize_t hashLength,
		    ssize_t *fullVector, ssize_t totalLength) {
  ssize_t i, j, pos;
  ssize_t *vector;
  DNASequence *fragment;
  pos = 0;
  for (i = 0; i < fragmentVectors.size(); i++) {
    fragment = fragments[i];
    vector = fragmentVectors[i];
    assert(i < fragments.size());
    if (fragment->length > 0)
      assert(pos < totalLength);
    for (j = 0; j < (fragment->length - hashLength+1); j++) {
      fullVector[pos] = vector[j];
      pos++;
    }
  }
}


void JoinStrings(std::vector<char*> &strings, 
		 std::vector<ssize_t> &lengths, 
		 char *&seq,
		 ssize_t  &len) {

  ssize_t i;

  // Find out the total length of the strings
  ssize_t totalLength = 0;
  for (i = 0; i < lengths.size(); i++) {
    totalLength += lengths[i];
  }

  totalLength += lengths.size() - 1;


  seq = (char*) new char[totalLength];
  ssize_t curPos = 0;
  for (i = 0; i < strings.size()-1; i++ ){ 
    strcpy(&(seq[curPos]), strings[i]);
    curPos += lengths[i];
    seq[curPos] = 5; // Create break in sequence;
    curPos++;
  }

  // Copy last sequence, don't need a terminating N
  strcpy(&(seq[curPos]), strings[i]);
}

void CreateMaskingIndices(T_Sequences &fragments,
			  std::vector<ssize_t> genomeIndices,
			  std::vector<ssize_t> frag1Indices,
			  std::vector<ssize_t> frag2Indices) {

  // Output a vector containing all the indices to mask
  
}

void DistributeMasking(T_Sequences &fragments, 
		       T_Params &params,
		       T_ParallelParams &parInfo) {

  ssize_t i, j;
  // In order to mask a sequence, each fragment needs to be 
  // masked against all other fragments, for nC2 comparisons.  
  //
  std::vector<ssize_t> idx1, idx2;
  ssize_t index1, index2;
  DNASequence resultSeq1, resultSeq2;
  ssize_t resultIndex1, resultIndex2;

  // Hold the indices of the resulting fragmentsn
  ssize_t childId;
  // Create a queue of jobs to run.
  for (i = 0; i < fragments.size() - 1; i++) 
    for (j = i+1; j < fragments.size(); j++ ) {
      idx1.push_back(i);
      idx2.push_back(j);
    }
  
  ssize_t numPendingJobs;
  ssize_t jobIndex;
  
  numPendingJobs = 0;
  jobIndex       = 0;
  // Submit as many jobs as possible to children
  for (i = 0; i < parInfo.numChildren && jobIndex < idx1.size(); i++) {
    childId = parInfo.childIds[i];
    index1 = idx1[jobIndex];
    index2 = idx2[jobIndex];
    jobIndex++;

    std::cout << childId << " sending " << index1 << " " << index2 << std::endl;
    // Submit a job when I figure out how to do this.
    CommunicationLib::InitiateTask(childId, doUniqueMask);
    CommunicationLib::SendMaskingInformation(childId, params.neighborDist,
					     params.doIndel, params.hashLength, params.wordLength);  
    CommunicationLib::SendSequences(parInfo.mytid, 
				    index1, index2,
				    *fragments[index1],
				    *fragments[index2],
				    childId,
				    CommunicationLib::msgUniqueMaskCompute);
    numPendingJobs++;
  }

  while (numPendingJobs > 0) {
    // Wait for a job to finish when I figure out how.
    childId = CommunicationLib::ReceiveSequences(CommunicationLib::idAny,
						 CommunicationLib::msgUniqueMaskResult,
						 &index1, &index2,
						 resultSeq1,
						 resultSeq2);
    std::cout << childId << " received " << index1 << " " << index2 << std::endl;
    ssize_t numLookups, numReferences;
    CommunicationLib::ReceivePerformanceReport(childId, 
					       &numLookups, &numReferences);
    
    numPendingJobs--;

    
    // Merge the results into the known sequence
    fragments[index1]->MergeSequence(resultSeq1);
    fragments[index2]->MergeSequence(resultSeq2);
    
    if (jobIndex < idx1.size()) {
      numPendingJobs++;
      index1 = idx1[jobIndex];
      index2 = idx2[jobIndex];
      jobIndex++;
      std::cout << childId << " sending " << index1 << " " << index2 << std::endl;
      CommunicationLib::InitiateTask(childId, doUniqueMask);
      CommunicationLib::SendMaskingInformation(childId, 
					       params.neighborDist, 
					       params.doIndel, 
					       params.hashLength, 
					       params.wordLength);
      CommunicationLib::SendSequences(parInfo.mytid, 
				      index1, index2,
				      *fragments[index1],
				      *fragments[index2],
				      childId,
				      CommunicationLib::msgUniqueMaskCompute);
    }
  }    
}


void DistributeUnmaskShared(T_GenomeVect &genomeFragments, 
			    T_Params &params,
			    T_ParallelParams &parInfo) {
  ssize_t i, j, g;
  // In order to mask a sequence, each fragment needs to be 
  // masked against all other fragments, for nC2 comparisons.  
  //
  std::vector<ssize_t> idx0, idx1;
  ssize_t index0, index1;
  DNASequence resultSeq0, resultSeq1;
  ssize_t resultIndex1, resultIndex2;
  T_Sequences *fragments;
  // Hold the indices of the resulting fragmentsn
  ssize_t childId;
  // Create a queue of jobs to run.
  assert ("can only share between two fragments " &&
	  genomeFragments.size() == 2); 

  for (i = 0; i < genomeFragments[0]->size(); i++) 
    for (j = 0; j < genomeFragments[1]->size(); j++) {
      idx0.push_back(i); // index into genomeFragments[0]
      idx1.push_back(j); // index into genomeFragments[1]
    }

  ssize_t numPendingJobs;
  ssize_t jobIndex;  
  numPendingJobs = 0;
  jobIndex       = 0;
  // Submit as many jobs as possible to children
  for (i = 0; i < parInfo.numChildren && jobIndex < idx1.size(); i++,jobIndex++) {
    // get information about the type of job to do
    childId = parInfo.childIds[i];
    index0 = idx0[jobIndex];
    index1 = idx1[jobIndex];

    // Submit a jobp
    CommunicationLib::InitiateTask(childId, doCommonMask);
    CommunicationLib::SendCommonInformation(childId, params.hashLength, params.wordLength);  
    CommunicationLib::SendSequences(parInfo.mytid, 
				    index0, index1,
				    *(*genomeFragments[0])[index0],
				    *(*genomeFragments[1])[index1],
				    childId,
				    CommunicationLib::msgCommonMaskCompute);

    // Remember how many jobs are out
    numPendingJobs++;
  }
  
  while (numPendingJobs > 0) {
    // Wait for a job to finish
    childId = CommunicationLib::ReceiveSequences(CommunicationLib::idAny,
						 CommunicationLib::msgCommonMaskResult,
						 &index0, &index1,
						 resultSeq0,
						 resultSeq1);

    numPendingJobs--;
    
    // Merge the results into the known sequence
    (*genomeFragments[0])[index0]->MergeSequence(resultSeq0, 0);
    (*genomeFragments[1])[index1]->MergeSequence(resultSeq1, 0);

    if (jobIndex < idx1.size()) {
      numPendingJobs++;
      index0 = idx0[jobIndex];
      index1 = idx1[jobIndex];
      jobIndex++;
      CommunicationLib::InitiateTask(childId, doCommonMask);
      CommunicationLib::SendCommonInformation(childId, params.hashLength, params.wordLength);
      CommunicationLib::SendSequences(parInfo.mytid, 
				      index0, index1,
				      *(*genomeFragments[0])[index0],
				      *(*genomeFragments[1])[index1],
				      childId,
				      CommunicationLib::msgCommonMaskCompute);
    }
  }    
}

void DistributeEnumerateGenome(T_GenomeVect &genomeFragments, 
			       T_IntegralValueVect &permutedEnumerations,
			       T_IntegralValueVect &permutedLocations,
			       ssize_t hashLength, ssize_t wordLength,
			       ssize_t myid,
			       ssize_t numChildren, ssize_t *childIds) { 

  ssize_t i, j;
  std::vector<ssize_t> idx0, idx1;
  ssize_t index0, index1;
  DNASequence resultSeq0, resultSeq1;
  ssize_t childId;
  ssize_t *receivedEnumerations, *receivedLocations; 
  ssize_t lenEnum0, lenEnum1;
  DNASequence *seq0, *seq1;

  ssize_t numPendingJobs;
  ssize_t jobIndex;
  for (i = 0; i < genomeFragments[0]->size(); i++) 
    for (j = 0; j < genomeFragments[1]->size(); j++) {
      idx0.push_back(i);
      idx1.push_back(j);
    }

  numPendingJobs = 0;
  jobIndex = 0;
  for (i = 0; i < numChildren && jobIndex < idx1.size(); i++,jobIndex++) {
    // get information about the type of job to do
    childId = childIds[i];
    index0 = idx0[jobIndex];
    index1 = idx1[jobIndex];

    // Submit a job

    CommunicationLib::InitiateTask(childId, doEnumerate);
    CommunicationLib::SendCommonInformation(childId, hashLength, wordLength);  
    CommunicationLib::SendSequences(myid, 
				    index0, index1,
				    *(*genomeFragments[0])[index0],
				    *(*genomeFragments[1])[index1],
				    childId,
				    CommunicationLib::msgEnumerateCompute);

    // Remember how many jobs are out
    numPendingJobs++;
  }
  
  while (numPendingJobs > 0) {
    // Wait for a job to finish
    childId = CommunicationLib::ReceiveEnumeration(CommunicationLib::idAny,
						   receivedEnumerations,
						   receivedLocations,
						   &index0, &index1);
    
    std::cout << "received: " << index0 << " " << index1 << std::endl;

    MergeVectors(permutedEnumerations[index0], 
		 receivedEnumerations, 
		 (*genomeFragments[0])[index0]->length);

    MergeVectors(permutedLocations[index0], 
		 receivedLocations, 
		 (*genomeFragments[0])[index0]->length, -1);
    
    delete receivedEnumerations;
    delete receivedLocations;
    
    numPendingJobs--;
    
    // Check out the result
    // Merge the results into the known sequence

    if (jobIndex < idx1.size()) {
      numPendingJobs++;

      index0 = idx0[jobIndex];
      index1 = idx1[jobIndex];

      jobIndex++;

      CommunicationLib::InitiateTask(childId, doEnumerate);
      CommunicationLib::SendCommonInformation(childId, hashLength, wordLength);
      CommunicationLib::SendSequences(myid, 
				      index0, index1,
				      *(*genomeFragments[0])[index0],
				      *(*genomeFragments[1])[index1],
				      childId,
				      CommunicationLib::msgEnumerateCompute); 
    }
  }    
}
    
void MergeVectors(ssize_t *vect1, ssize_t *vect2, ssize_t len, ssize_t unsetValue) {
  ssize_t i;
  for (i = 0; i < len; i++) {
    if (vect2[i] != unsetValue)
      vect1[i] = vect2[i];
  }
}

void FragmentGenomes(T_GenomeVect &genomes,
		     T_GenomeVect &genomeFragments, 
		     ssize_t numFragments,
		     T_Params &params) {
  ssize_t i;
  T_Sequences *fragments;
  // The blank fragment allows for masking a sequence for repeats within iteslf.
  for (i = 0; i < genomes.size(); i++) {
    // Fragment each genome into children fragments.
    DNASequence *blankFragment = new DNASequence;
    fragments = (T_Sequences*) new T_Sequences;
    FragmentOneGenome(*(*genomes[i])[0],
		      numFragments, 
		      params.hashLength, 
		      *fragments  // store individual fragments here
		      );
    fragments->push_back(blankFragment);
    genomeFragments.push_back(fragments);
  }
}

void MaskFragments(T_GenomeVect &genomeFragments, 
		   T_Params params) {
  ssize_t g, i, j;
  T_Sequences *fragments;
  for (g = 0; g < genomeFragments.size(); g++ ){
    fragments = genomeFragments[g];
    //  Compare all pairs of fragments
    for (i = 0; i < fragments->size() -1; i++) {
      for (j = i + 1; j < fragments->size(); j++) {
	MaskNotUnique(*(*fragments)[i], *(*fragments)[j], 
		      params.neighborDist, params.doIndel, 
		      params.hashLength, params.wordLength);
      }
    }
  } 
}

void InitializeNotCommon(T_GenomeVect &genomeFragments, T_Params params) {
  T_Sequences *fragments;
  DNASequence  *fragment;
  ssize_t g, i, j;
  for (g = 0; g < genomeFragments.size(); g++) {
    fragments = genomeFragments[g];
    for (i = 0; i < fragments->size(); i++) {
      fragment = (*fragments)[i];
      for (j = 0; j < fragment->length - params.hashLength + 1; j++) {
	if ( ! fragment->IsPosMasked(j) )
	  fragment->MarkPosNotCommon(j);
      }
    }
  }
}

void UnmaskSharedGenomes(T_GenomeVect &genomeFragments, T_Params &params) {
   ssize_t i, j;
   for (i = 0; i < genomeFragments[0]->size(); i++) 
     for (j = 0; j < genomeFragments[1]->size(); j++) 
       UnmaskShared(*((*genomeFragments[0])[i]), *((*genomeFragments[1]))[j], 
		    params.hashLength, params.wordLength);
}

void DistributeMaskingGenomes(T_GenomeVect genomeFragments, 
			      T_Params &mainParams, 
			      T_ParallelParams &parInfo) {
  ssize_t g;
  for (g = 0; g < genomeFragments.size(); g++ ){
    DistributeMasking(*(genomeFragments[g]), 
		      mainParams,
		      parInfo);
  }
}

ssize_t EnumerateSequences(DNASequence &refSeq, DNASequence &qrySeq, 
			T_Params &params,
			T_ParallelParams &parInfo,
			T_Alignment *&alignment) {
  alignment->length = refSeq.length;
  return EnumerateSequences(refSeq, qrySeq, params, parInfo,
			    alignment->enumerations, alignment->locations);
}


ssize_t EnumerateSequences(DNASequence &refSeq, DNASequence &querySeq, 
			T_Params &mainParams,
			T_ParallelParams &parInfo,
			ssize_t *&enumerations,
			ssize_t *&locations) {
  ssize_t g, i, j;
  // Initiailize the genoms vect.  This is because I didn't plan on 
  // having just two sequences, but for now there are.

  T_Sequences refGenome, queryGenome;
  refGenome.push_back(&refSeq);
  queryGenome.push_back(&querySeq);

  T_GenomeVect genomes;
  genomes.push_back(&refGenome);
  genomes.push_back(&queryGenome);

  // Temporary variables for dealing with many small fragments
  T_GenomeVect genomeFragments;
  T_Sequences *fragments;
  T_Locations *fragmentLocations;

  // Step 0.1  Create the fragments that will be processed.
  FragmentGenomes(genomes, genomeFragments, parInfo.numFragments, mainParams);
  
  if (parInfo.numChildren == 0) {
    //  Step 1. Mask out all nucleotides that correspond to the start of
    //  a non-unique tuple.
    MaskFragments(genomeFragments, mainParams);
      
    // Step 2. Find all common tuples between all input sequences (just do 2 for now).  
    //   Step 2.1 Initialize all positoins that are not already masked to be 
    //            not common.  They will be marked as common later.
    InitializeNotCommon(genomeFragments, mainParams);

    //   Step 2.2 Search through all fragments of each genome to find the
    //   common nucleotides.  
    UnmaskSharedGenomes(genomeFragments, mainParams);
  }
  else {
    // Parallel version of the above.
    DistributeMaskingGenomes(genomeFragments, mainParams, parInfo);
    InitializeNotCommon(genomeFragments, mainParams);
    DistributeUnmaskShared(genomeFragments, mainParams, parInfo);
  }

  // Step 3.4.1. Get the new sequence
  for (i = 0; i < genomeFragments.size(); i++) {
    UnfragmentGenome(*(genomeFragments[i]), mainParams.hashLength,
		     *((*genomes[i])[0]));
  }


  // Step 3.  Enumerate all shared common tuples in each fragment.
  //    Step 3.1 Count the number of unique (non-masked) tuples in each
  //    fragment.   
  DNASequence *fragment;
    
  ssize_t numUnmasked;
  for (g = 0; g < genomeFragments.size(); g++ ) {
    fragments= genomeFragments[g];
    numUnmasked = 0;
    for (i = 0; i < fragments->size(); i++) {
      fragment = (*fragments)[i];
      fragment->startEnumeration = numUnmasked;
      for (j = 0; j < ((*fragments)[i]->length - mainParams.hashLength + 1); j++) {
	if (! DNASequence::IsMasked(fragment->seq[j])) {
	  fragment->numUnmasked++;
	  numUnmasked++;
	}
      }
    } // end looping over all fragments
  } // end looping over all genomes.
    

  // Step 3.2 Allocate arrays that will store the enumerations of the
  // permuted sequence.
    
  ssize_t *permutedEnumeration;
  ssize_t *permutedLocation;
  T_IntegralValueVect permutedEnumerations;
  T_IntegralValueVect permutedLocations;
  fragments =  genomeFragments[0];
  for (i = 0; i < fragments->size(); i++) {
    permutedEnumeration = new ssize_t[(*fragments)[i]->length];
    permutedLocation    = new ssize_t[(*fragments)[i]->length];
    permutedEnumerations.push_back(permutedEnumeration);
    permutedLocations.push_back(permutedLocation);
    for (j = 0; j < (*fragments)[i]->length; j++) {
      permutedEnumeration[j] = 0;
      permutedLocation[j] = -1;
    }
  }

  // Step 3.3 Fill the permuted enumerations array by comparing all
  // the fragments
  ssize_t l;
  if (parInfo.numChildren <= 1) {
    for (i = 0; i < genomeFragments[0]->size(); i++) 
      for (j = 0; j < genomeFragments[1]->size(); j++) {
	EnumerateUnique(*((*genomeFragments[0])[i]),
			*((*genomeFragments[1]))[j], 
			mainParams.hashLength, mainParams.wordLength, 
			permutedEnumerations[i],
			permutedLocations[i]
			);
      }
  }
  else {
    // send out the sequences to workers
    DistributeEnumerateGenome(genomeFragments, 
			      permutedEnumerations,
			      permutedLocations,
			      mainParams.hashLength, mainParams.wordLength,
			      parInfo.mytid,
			      parInfo.numChildren,  parInfo.childIds);
  }

  
  // Step 3.4 Combine the results of all permuted enumerations into
  // one.


  // Step 3.4.1. Get the new sequence
  for (i = 0; i < genomeFragments.size(); i++)
    UnfragmentGenome(*(genomeFragments[i]), mainParams.hashLength,
		     *((*genomes[i])[0]));
    

  // Step 3.4.2. Combine enumerations.

  ssize_t length;
  length = (*genomes[0])[0]->length;
  enumerations = new ssize_t[length];
  locations    = new ssize_t[length];

  for (i = 0; i < length; i++) {
    locations[i] = -1;
    enumerations[i] = 0;
  }
  
  CombineVectors(permutedEnumerations,
		 *genomeFragments[0],
		 mainParams.hashLength,
		 enumerations, (*genomes[0])[0]->length);
  
  CombineVectors(permutedLocations,
		 *genomeFragments[0],
		 mainParams.hashLength,
		 locations, (*genomes[0])[0]->length);
  ssize_t countEnumerations = 0;
  for (i = 0; i < length; i++) {
    if (locations[i] != -1) {
      ++countEnumerations;
    }
  }

  ssize_t f;
  for (g = 0; g < genomeFragments.size(); g++) {
    for (f = 0; f < (*genomeFragments[g]).size(); f++) {
      delete (*genomeFragments[g])[f];
    }
    delete genomeFragments[g];
  }

  for (i = 0; i < permutedEnumerations.size(); i++){
    delete[] permutedEnumerations[i];
  }
  for (i = 0; i < permutedLocations.size(); i++){
    delete[] permutedLocations[i];
  }

  return countEnumerations;
}

void FindNextGap( DNASequence &seq, 
		  T_Params &params,
		  ssize_t *enumerations,
		  ssize_t *locations,
		  ssize_t minGapLength,
		  double minGapRatio,
		  ssize_t curRefPos,         // index into reference sequence
		  ssize_t curQryPos,         // index into query enumeration and location arrays
		  ssize_t &nextGapRefStart,  // position in the reference genome where there is a gap between 
                                         // anchors. 
		  ssize_t &nextGapQryStart   // position in the
		                         // enumerations and locations of the next gap in  
		                         //   query genome.  
		  ) {
  ssize_t i;
  ssize_t refEnumeration;
  refEnumeration = 0;
  
  ssize_t curEnumeration = enumerations[curQryPos];
  ssize_t curPosition    = locations[curQryPos];

  // Scan forward in the reference sequence and consider the locations 
  // of the anchors.  When the difference between consecutive 
  for (curRefPos; curRefPos < (seq.length - params.hashLength + 1); curRefPos++) {
    if (seq.IsACTG(i) && !seq.IsPosMasked(i)) {
      /* Sanity check. */
      assert("should not have an enumeration at a masked pos " 
	     && enumerations[curRefPos] != 0);
      assert("should have assigned a location at unmasked pos "
	     && locations[curRefPos] != -1);

      refEnumeration++;
    }
  }
}


void CreateStripArrays(DNASequence &refSeq,  // description of first seq
		       // permutation of the second sequence
		       T_Params params,
		       ssize_t *& seqEnumerations, ssize_t *& seqLocations, 
		       // Output
		       ssize_t & numTuples,
                       //    condensed enumerations
		       ssize_t *& enumerations,  
		       //    Condensed locations
		       ssize_t *& refLocations, ssize_t *&qryLocations
		       ) {
  ssize_t i, e;
  // Count the number of positions.
  numTuples = 0;
  for (i = 0; i < (refSeq.length - params.hashLength + 1); i++)
    if (refSeq.IsACTG(i) && !refSeq.IsPosMasked(i))
      numTuples++;

  ssize_t numLocations = 0;
  for (i = 0; i < (refSeq.length - params.hashLength + 1); i++)
    if (seqLocations[i] != -1) 
      numLocations++;

  assert(numLocations == numTuples);

  enumerations = new ssize_t[numTuples];
  refLocations = new ssize_t[numTuples];
  qryLocations = new ssize_t[numTuples];

  e = 0;
  for (i = 0; i < (refSeq.length - params.hashLength + 1); i++) {
    if (refSeq.IsACTG(i) && !refSeq.IsPosMasked(i)) {
      if(seqEnumerations[i] == 0 || seqLocations[i] == -1) {
	std::cout << "bailing out " << seqEnumerations[i] << " " 
		  << seqLocations[i] << std::endl;
	pvm_exit();
	exit(0);
      }
      enumerations[e] = seqEnumerations[i];
      refLocations[e] = i;
      qryLocations[e] = seqLocations[i];
      e++;
    }
  }
}

ssize_t FindBetweenProbes(DNASequence &refSeq,
		      ssize_t *enumerations, ssize_t *locations,
		      ssize_t minDist, ssize_t maxDist,
		      ssize_t &refPos, ssize_t &refStartPos, ssize_t &refEndPos) 
{
  ssize_t foundGap;
  foundGap = 0;
  while (! foundGap && refPos < refSeq.length) {
    while (refPos < refSeq.length && locations[refPos] == -1)
      refPos++;
    refEndPos = refStartPos = refPos;
  
    refPos++;
    while (refPos < refSeq.length && locations[refPos] == -1)
      refPos++;
    
    refEndPos = refPos;
    
    if (refEndPos - refStartPos > minDist && refPos < refSeq.length)
      foundGap = 1;
  }
  return foundGap;
}
  

ssize_t FindGap(DNASequence &refSeq, 
	     ssize_t *enumerations, ssize_t *locations, 
	     ssize_t minDist, ssize_t maxDist, 
	     ssize_t &refPos, 
	     ssize_t &refStartPos, ssize_t &refEndPos ) {
  // Given a list of enumerations and locations, look for two mapped
  // locations that have adjacent enumerations but are separated by some 
  // number of nucleotides.

  ssize_t startEnumeration = -1;
  ssize_t endEnumeration;
  // Located the boundaries of the next consecutively enumerated
  // nucleotides.
  ssize_t foundGap = 0;
  while( ! foundGap && refPos < refSeq.length) {
    while (refPos < refSeq.length && (!refSeq.IsACTG(refPos) || refSeq.IsPosMasked(refPos))) 
      refPos++;
    
    refStartPos = refPos;
    refPos++;
    
    while (refPos < refSeq.length && (!refSeq.IsACTG(refPos) || refSeq.IsPosMasked(refPos))) 
      refPos++;

    // Can't say there is a gap when there was no matching
    // tuple at the end of the sequence.
    if (refPos == refSeq.length)
      return 0;

    refEndPos = refPos;

    startEnumeration = enumerations[refStartPos];
    endEnumeration   = enumerations[refEndPos];
    
    if (startEnumeration + 1 == endEnumeration)
      if (locations[refEndPos] - locations[refStartPos] > minDist  &&
	  refEndPos - refStartPos > minDist)
	foundGap = 1;
  }
  refPos++;
  return foundGap;
}


/* 
   Add all of this later.
  // Condense the enumerations and locations.
  CreateStripArrays(*genomes[0][0][0],  // description of first seq
		    // permutation of the second sequence
		    enumerations, locations, 
		    qrySeqArray,
		    refLocations, qryLocations
		    );

  GetMUMs(*refSeq, qrySeqArray, refLocations, qryLocations, 
	  refSeq->length, strips);

*/

void AddAlignment(T_Alignment *&parent, T_Alignment *&current, T_Alignment *recAlignment) {
  // Store the alignment in a chain.
  if (recAlignment == NULL) 
    return;
  assert(current != recAlignment);
  if (parent->child == NULL) {
    // If this is the first alignment on a new level, the
    // parent alignment needs a pointer to  the first child
    // alignment. The child is then chained to all of it's
    // subsequent alignment (stripIt).
    parent->child = recAlignment;
  }
  else {
    assert("recAlignment is NULL??? " && recAlignment != NULL);
    current->next = recAlignment;
  }
  current = recAlignment;
  current->parent = parent;
}

double EndAlign2Seq(DNASequence &refSeq, DNASequence &qrySeq,
		  T_Alignment *&alignment, ssize_t dir) {
  // End alignment means finding the optimal alignment of the shorter
  // sequence in the longer sequence such that the ends of the sequences
  // align together, and there is one large gap allowed in the 
  // long sequence.
  ssize_t i, j;
  ssize_t band = 10;
  ssize_t *forLocations, *revLocations;
  double *forScores /* and seven years ago */, *revScores;
  DNASequence rcQrySeq;
  if (dir == -1) {
    MakeRC(qrySeq, rcQrySeq);
    if (alignment != NULL)
      delete alignment;
    alignment = NULL;
    return 0;
  }

  alignment->length = refSeq.length;
  alignment->locations = new ssize_t[refSeq.length];
  alignment->enumerations = new ssize_t[refSeq.length];
  alignment->hashLength = 1;
  refSeq.UnmaskNoCk();
  qrySeq.UnmaskNoCk();
  for (i = 0; i < alignment->length; i ++) {
    alignment->locations[i] = -1;
    alignment->enumerations[i] = 0;
  }

  DNASequence *longSeqPtr;
  DNASequence *shortSeqPtr;

  
  // Determine which sequence is longest 
  if (refSeq.length > qrySeq.length) {
    longSeqPtr = &refSeq;
    if (dir == 1) {
      shortSeqPtr = &qrySeq;
    }
    else {
      shortSeqPtr = &rcQrySeq;
    }
  }
  else {
    shortSeqPtr = &refSeq;
    if (dir == 1) {
      longSeqPtr = &qrySeq;
    }
    else {
      longSeqPtr = &rcQrySeq;
    }
  }
  
  // Align shorter sequence against part of longer sequence
  // the longer sequence that is aligned is the length of the shorter sequence
  // + the k-band.  The motivation behind this is that if the two
  // sequences differ by less than k indels, the optimal alignment
  // will be found here.
  forLocations = new ssize_t[shortSeqPtr->length + band];
  revLocations = new ssize_t[shortSeqPtr->length + band];
  forScores  = new double[shortSeqPtr->length + band];
  revScores  = new double[shortSeqPtr->length + band];

  // Initialize the alignment scores to unassigned (-1)
  for (i = 0; i < shortSeqPtr->length + band; i++) {
    forLocations[i] = revLocations[i] = -1;
    forScores[i] = revScores[i] = INFTY;
  }
  
  /**********************************************************/
  // Perform end alignment of the two sequences
  //
  BandedAlign(*longSeqPtr, *shortSeqPtr,
	      -4, 5, 5,  // assume the short sequence fits here, so 
                         // penalize for gaps.  todo: make biologically
                         // significant.
	      band,
	      forLocations, forScores, NULL);

  // Align in the reverse direction, starting from the ends of 
  // each string.
  DNASequence shortRev(shortSeqPtr->length);
  DNASequence longRev(shortSeqPtr->length+band);

  // Find the reverse of each sequence.  Note this is not 
  // the reverse complement, just the reverse (although the reverse
  // complement would do just same here with a little extra work).
  shortRev.StoreName("shortseq rev");
  shortRev._ascii= 1;
  longRev._ascii = 1;
  for (i = 0; i < shortSeqPtr->length; i++)
    shortRev.seq[i] = shortSeqPtr->seq[shortSeqPtr->length-1-i]; 
  
  for (i = 0; i < shortSeqPtr->length + band; i++) 
    longRev.seq[i] = longSeqPtr->seq[longSeqPtr->length - 1 - i];

  BandedAlign(longRev, shortRev,
	      -4, 5, 5,
	      band,
	      revLocations, revScores, NULL);

  
  // Now find the optimal combination of the two sequences
  double minScore = INFTY;
  ssize_t revdx;
  ssize_t minRev, minFor;
  minRev = minFor = -1;
  // Look for a minimum score that isn't a combination of
  // the two sequences (the optimal alignment simply fits 
  // one sequence into the end of the other ).
  for (i = 0; i < shortSeqPtr->length+band-2; i++) {
    if (minScore < forScores[i]) {
      minFor = i;
      minRev = -1;
    }
  }
  for (i = 0; i < shortSeqPtr->length+band-2; i++) {
    if (minScore < revScores[i]){ 
      minRev = i;
      minFor = -1;
    }
  }

  // Look for an optimal score that is a combination of the
  // two alignments with one gap in the long sequence, and 
  // possibly a gap in the shorter sequence.
  for (i = 0; i < shortSeqPtr->length+band-2-1; i++) {
    for (j = 0; j < shortSeqPtr->length+band-2-i-1; j++) {
      revdx = shortSeqPtr->length+band-1-j-1;
      if (forScores[i] + revScores[revdx] < minScore) {
	minScore = forScores[i] + revScores[revdx];
	minFor = i;
	minRev = j;
      }
    }
  }
  /*
   *******************************************
   Done aligning the ends of the two sequences.  Now combine the two alignments 
   into one.
   */
  
  if (qrySeq.length < refSeq.length) {
    // Now forLocations map from the short (query) sequence, and 
    // revLocations map from the reverse short sequence into the end of the 
    // reverse long sequence.  

    if (dir == 1) {
      for (i = 0; i <= minFor ; i++)
	alignment->locations[i] = forLocations[i];
      for (i = 0; i <= minRev; i++) {
	if (revLocations[i] != -1)
	  alignment->locations[longSeqPtr->length-1-i] = 
	    shortSeqPtr->length -1 - revLocations[i];
      }
    }
    else {
      // Locations are mapped to reverse complement, so unmap them
      for (i = 0; i <= minFor ; i++)
	alignment->locations[i] = shortSeqPtr->length - forLocations[i] - 1;
      for (i = 0; i <= minRev; i++) {
	if (revLocations[i] != -1)
	  /*	  explicit version of next statment 
	    alignment->locations[longSeqPtr->length-1-i] = 
	    shortSeqPtr->length-1- (shortSeqPtr->length -1 - revLocations[i]);*/
	  alignment->locations[longSeqPtr->length-1-i] = revLocations[i];
      }
    }
  }
  else {
    if (dir == 1) {
      // The qrySeq was the longer seq, and the short seq 
      for (i = 0; i < minFor; i++) {
	if (forLocations[i] != -1)
	  alignment->locations[forLocations[i]] = i;
    }
      // error here
      for (i = 0; i < minRev ; i++) {
	if (revLocations[i] != -1)
	alignment->locations[shortSeqPtr->length - 1 - revLocations[i]] = 
	  longSeqPtr->length - i - 1;
      }
    }
    else {
      // The qrySeq was the longer seq, and the short seq is reference
      for (i = 0; i < minFor; i++) {
	if (forLocations[i] != -1)
	  alignment->locations[forLocations[i]] = rcQrySeq.length-1 - i;
      }
      // error here
      for (i = 0; i < minRev ; i++) {
	if (revLocations[i] != -1)
	  /*	  explicit version of next statment:
	    alignment->locations[shortSeqPtr->length - 1 - revLocations[i]] = 
	    rcQrySeq.length - 1 - (longSeqPtr->length - i - 1);
	  */
	  alignment->locations[shortSeqPtr->length - 1 - revLocations[i]] = i;
      }
    }
  }
  // Now that the alignment is created, the sequence should be masked
  // and enumerations assigned according to the alignment.
  for (i = 0; i < refSeq.length; i++) {
    if (alignment->locations[i] == -1) {
      alignment->enumerations[i] = 0;
      refSeq.IsACTG(i) && refSeq.MarkPosNotCommon(i);
    }
    else {
      alignment->enumerations[i] = dir;
      refSeq.IsACTG(i) && refSeq.UnmaskPos(i);
    }
  }

  for (i = 0; i < refSeq.length; i++) {
    if ((!refSeq.IsACTG(i) && alignment->locations[i] != -1) || !refSeq.IsPosMasked(i)) {
      if (refSeq.IsACTG(i)) 
	assert(alignment->enumerations[i] != 0 && 
	       alignment->locations[i] != -1);
    }
    else {
      assert(alignment->locations[i] == -1);
      assert(alignment->enumerations[i] == 0);
    }
  }
  
  delete[] forLocations;
  delete[] revLocations;
  delete[] forScores;
  delete[] revScores;
  delete[] shortRev.seq;
  delete[] longRev.seq;
  return minScore;
}


void AffineAlign2Seq(DNASequence &refSeq, DNASequence &qrySeq, 
		     T_Alignment*&alignment,
		     T_Params params, std::string padding,
		     ssize_t dir) {
  ssize_t i; 
  // Find out what strand was aligned to get this fragment. 
  ssize_t *locations = new ssize_t[refSeq.length];
  ssize_t strand = 1;
  double score;
  if (dir == 1) {
    score = AffineAlign(refSeq, qrySeq, -1, 1, 2, 6, 0, locations, NULL);
    alignment->locations = locations;
  }
  else {
    DNASequence qryRCSeq;
    MakeRC(qrySeq, qryRCSeq);
    score = AffineAlign(refSeq, qryRCSeq, -1, 1, 2, 6, 0, locations, NULL);
    alignment->locations = locations;
    ssize_t i;
    // Translate the reverse aligned positions to forward
    RotateAlignment(alignment->locations, refSeq.length, qrySeq.length);
    strand = -1;
    delete []qryRCSeq.seq;
  } 
  // Store the enumerations
  alignment->length = refSeq.length;
  alignment->enumerations = new ssize_t[refSeq.length];
  alignment->hashLength = 1;
  refSeq.UnmaskNoCk(0, refSeq.length);
  qrySeq.UnmaskNoCk(0, qrySeq.length);
  for (i = 0; i < refSeq.length ; i++) {
    if (alignment->locations[i] != -1) 
      alignment->enumerations[i] = strand;
    else 
      alignment->enumerations[i] = 0; 
  } 
  // Mask unaligned characters
  for ( i = 0; i < refSeq.length; i++) 
    if (alignment->locations[i] == -1)
      refSeq.IsACTG(i) && refSeq.MarkPosNotCommon(i);
}

T_Alignment* AlignSubsequences(DNASequence &refSeq, DNASequence &qrySeq,
			       ssize_t startRefPos, ssize_t endRefPos,
			       ssize_t startQryPos, ssize_t endQryPos,
			       ssize_t dir,
			       T_Alignment *alignment,
			       T_Params params,
			       T_ParallelParams parInfo, std::string padding){
  T_Params recParams;
  ssize_t temp;
  DNASequence refSubSeq, qrySubSeq;
  T_Alignment *recAlignment = NULL;
  // Initialize the subsequences
  assert(startRefPos < endRefPos || (printf("aligning reference sequence in reverse direction") == 0));
  //  refSubSeq.length = endRefPos - startRefPos + params.hashLength;
  refSubSeq.length = endRefPos - startRefPos;
  refSubSeq.seq    = &refSeq.seq[startRefPos] ;
  refSubSeq.startPosition = 0;
  
  if (startQryPos > endQryPos) {
    temp = endQryPos;
    endQryPos = startQryPos;
    startQryPos = temp;
  }
  qrySubSeq.length = endQryPos - startQryPos;
  
  qrySubSeq.seq = &qrySeq.seq[startQryPos];
  qrySubSeq.startPosition = 0;

  // Set the recursive parameters
  recParams = params;
  recParams.hashLength = std::min(std::min(refSubSeq.length, qrySubSeq.length) - 1, 
				  params.hashLength - params.hashStep);
  recParams.wordLength = (ssize_t) std::min(11.0, std::max(2.0, (log2(std::max(refSubSeq.length, qrySubSeq.length)) / 1.2)));

  if ( refSubSeq.length > params.distThreshold && 
       qrySubSeq.length > params.distThreshold ) {
    // The alignment computed for this region will be overwritten
    // by the (hopefully) more fine grained recursive alignment.  
    // Reset all the positions in the current alignment that will be recomputed.
    ssize_t i;
    for (i = 0; i < refSubSeq.length - params.hashLength; i++) {
      assert((i + startRefPos) < alignment->length ||
      (printf("inconsistent ref pos %d at %d for alignment of length: %d\n", 
	      startRefPos, i, alignment->length) == 0));
      alignment->locations[i+startRefPos]    = -1;
      alignment->enumerations[i+startRefPos] = 0;
    }
    
    EnumerateRecursively(refSubSeq, qrySubSeq, recAlignment, dir,
			 recParams, parInfo, padding + " ");
    
    if (recAlignment != NULL) {
      recAlignment->refPos = startRefPos;
      recAlignment->qryPos = startQryPos;
      recAlignment->level  = alignment->level + 1;
    }

    return recAlignment;
  }
  else {
    return NULL;
  }
}


void EnumerateRecursively(DNASequence &refSeq, DNASequence &qrySeq, 
			  T_Alignment *&alignment,
			  ssize_t dir,
			  T_Params params, T_ParallelParams &parInfo, 
			  std::string padding) {
  // Parameters say we cannot further align.

  /*
    Given a reference sequence and query sequence (refSeq, qrySeq), 
    find an enumeration of nucleotides of the mapping of each
    nucleotide in refSeq to qrySeq.
  */
  if (params.hashLength <= params.minHashLen) {
    std::cout << "not aligning sequencess of lengths: " << refSeq.length 
	      << "and " << qrySeq.length << std::endl;
    return;
  }

  ssize_t i;
  ssize_t *qryEnumerations, *refLocations, *qryLocations;
  if ( ! params.loadCheckpoint) {
    alignment = new T_Alignment;
    alignment->hashLength = params.hashLength;
    
    // Check to see if alignment should be done with dynamic programming
    // via either affine alignment (when both sequences are short)
    // or end aligning, when one sequence is long and the other is short.
    if (qrySeq.length / double(refSeq.length) < params.ratioThreshold
				||
				refSeq.length / double(qrySeq.length) < params.ratioThreshold) { 
      EndAlign2Seq(refSeq, qrySeq, alignment, dir); 
      return; 
    }
    if ( (double(refSeq.length) * double(qrySeq.length) < params.productSize) {
      /*      std::cout << padding << "affine aligning " << refSeq.length 
	      << " " << qrySeq.length << std::endl;*/
	
      AffineAlign2Seq(refSeq, qrySeq, alignment, params, padding, dir);
      return;
    }
    // It was notnot appropriate to do dynamic programming alignment.  Recursively
    // align finding common tuples (anchors).
    ssize_t countEnumerated; // result of how many positions in ref are mapped to qry
    countEnumerated = EnumerateSequences(refSeq, qrySeq, 
					 params, parInfo, alignment);

    if (countEnumerated == 0) {
      delete alignment;
      alignment = NULL;
      return;
    }
  }
  
  // For now, assume all recursive alignment is done in serial.
  parInfo.numFragments = 1;
  parInfo.numChildren  = 0;
#ifdef _DEBUG_VERSION
  // validate aligned sequences.
  for (i = 0; i < alignment->length - alignment->hashLength + 1; i++) {
    assert(alignment->locations[i] >= -1);
    if (refSeq.IsACTG(i)) {
      if(refSeq.IsPosMasked(i)) {
	assert(alignment->locations[i] == -1);
	assert(alignment->enumerations[i] == 0);
      }
      else {
	assert(alignment->locations[i] != -1 &&
	       alignment->enumerations[i] != 0);
      }
      if (alignment->enumerations[i] != 0)
	assert(alignment->locations[i] != -1);
    }
  }
#endif
  T_Alignment *recAlignment = NULL;
  T_Alignment *curAlignment = NULL;

  //  if (params.hashLength - params.hashStep >= params.minHashLen) { 
  DNASequence refSubSeq, qrySubSeq;
  ssize_t numTuples;
  ssize_t numDiscardedStrips;
  // Advance to the next smaller hash size.
  
  //    if (params.alignStrips) {
  // Transform the enumeration into strips
  // The first step for this is to compress the enumerations array
  // from per nucleotide enumeration to enumeration of locations.
  CreateStripArrays(refSeq,  // description of first seq
		    params, 
		    alignment->enumerations, alignment->locations,  
		    numTuples,
		    qryEnumerations,          //    condensed enumerations
		    refLocations, qryLocations//    Condensed locations
		    );
  T_Strips *strips = new T_Strips;
  GetStrips(qryEnumerations, refLocations, refLocations,
	    qryLocations, qryLocations, NULL,   
	    NULL, numTuples, *strips, params.hashLength, params.alignStrips);
  
  // This information is now contained in the strips structure
  delete []qryEnumerations;
  delete []refLocations;
  delete []qryLocations;
  /*
   * Go from postions to enumerated strips.
   */
  //  EnumerateStrips(*strips);
  ssize_t stripStart;
  T_StripQQry   stripsSortedByQry;
  T_Strips::iterator stripIt, nextIt;

  // Recursive alignment is alignment between two already anchored positions. 
  // Assume that the ends of the sequenes align, and are in proper
  // orientation. 
  stripIt = strips->begin();
  if (stripIt->startRefPos != 0) {
    std::cout << padding << "added start " << std::endl;
    // Insert a new strip that marks the beginning of the strand
    Strip beg;
    beg.startQryEnum = 0;
    beg.endQryEnum = 0;
    beg.endRefEnum = 0;
    beg.startRefEnum = 0;
    beg.startRefPos = 0;
    beg.startQryPos = 0;
    beg.endRefPos = params.wordLength;
    beg.endQryPos = params.wordLength;
    strips->push_front(beg);
  }
  stripIt = strips->end();
  --stripIt;
  if (stripIt->endRefPos != refSeq.length-1) {
    std::cout << padding << "added end " << std::endl;
    Strip end;
    end.startQryEnum = stripIt->endQryEnum+(size_t) +1;//1*Sign(stripIt->startQryEnum);
    end.endQryEnum   = stripIt->endQryEnum+(size_t) +1; //1*Sign(stripIt->startQryEnum);
    end.startRefEnum = stripIt->endRefEnum+1;
    end.endRefEnum   = stripIt->endRefEnum+1;
    end.startRefPos  = refSeq.length - params.hashLength - 1;
    end.endRefPos    = refSeq.length-1;
    end.startQryPos  = qrySeq.length - params.hashLength - 1;
    end.endQryPos    = qrySeq.length - 1;
    strips->push_back(end);
  }
  // Re-assign the enumeration to reflect the beginnign and end.
  //  EnumerateStrips(*strips);
  ssize_t matchExt;
  stripIt = strips->begin();
  while (stripIt != strips->end()) {
    if (Sign(stripIt->startQryEnum) > 0) 
      matchExt = ExtendMatch(refSeq, qrySeq, 
			     stripIt->startRefPos, stripIt->startQryPos,
			     Sign(stripIt->startQryEnum), 0.03);
    else
      matchExt = ExtendMatch(refSeq, qrySeq, 
			     stripIt->startRefPos, stripIt->startQryPos + params.hashLength - 1,
			     Sign(stripIt->startQryEnum), 0.03);

    /*    std::cout << "extended strip of length: " << stripIt->startRefPos << "  " 
	      << Sign(stripIt->startQryEnum)*stripIt->size << " ext: " 
	      << matchExt << std::endl;*/
    stripIt->size = matchExt;
    ++stripIt;
  }

  std::cout << "before min threading: " << std::endl;
  PrintEnumerations(*strips, std::cout);  
  if (strips->size() < 500)
    MinThreadingPath(strips, strips->size());
  else
    MaximizeStrips(strips, params.lisWindow, params.lisThreshold);

  ReannotateAlignment(refSeq, strips, alignment);
  
  // Since filtering has been performed, some strips have been
  // removed.   The resulting order of strips should be determined
  // so that adjacent strips can be aligned recursively.

  std::cout << "after min threading: " << std::endl;  
  PrintEnumerations(*strips, std::cout);

  EnumerateStrips(*strips); // Find which strips are adjacent to eachother.
  //  AlignStrips(refSeq, qrySeq, alignment, curAlignment, strips, params, parInfo, "");
  AlignBetweenStrips(refSeq, qrySeq, alignment, curAlignment, strips, params, parInfo, "");

  delete strips;
}

void AlignStrips(DNASequence &refSeq,
		 DNASequence &qrySeq,
		 T_Alignment *alignment, 
		 T_Alignment *&curAlignment,
		 T_Strips *&strips,
		 T_Params params,
		 T_ParallelParams &parInfo,
		 std::string padding) {
  ssize_t strand;
  ssize_t startRefPos, endRefPos;
  ssize_t startQryPos, endQryPos;
  T_Alignment *recAlignment;
  T_Strips::iterator stripIt;
  recAlignment = NULL;
  params.AlignGaps();
  // Refine the alignment between preserved k-mers
  ssize_t i;
  ssize_t s = 0;
  for (stripIt = strips->begin(); 
       stripIt != strips->end(); 
       ++stripIt){
    ++s;
    strand = Sign(stripIt->startQryEnum);
    recAlignment = NULL;
    // Reset the alignment in this region.
    startRefPos = stripIt->startRefPos; // params.hashLength;
    endRefPos   = stripIt->endRefPos;
    startQryPos = stripIt->startQryPos;// + 1*strand; //params.hashLength*strand;
    endQryPos   = stripIt->endQryPos;// - 1*strand;
    std::cout << "aligning strip: " << stripIt->startQryEnum << " "
	      << startRefPos << " " << endRefPos << " "
	      << startQryPos << " " << endQryPos << std::endl;
    if (startRefPos < endRefPos) {
      for (i = startRefPos; i < endRefPos; i++) {
	refSeq.MarkPosNotCommon(i);
	alignment->locations[i] = -1;
	alignment->enumerations[i] = 0;
      }
    }
    recAlignment = AlignSubsequences(refSeq, qrySeq,
				     startRefPos, endRefPos,
				     startQryPos, endQryPos,
				     strand,
				     alignment, params, parInfo, padding);
    AddAlignment(alignment, curAlignment, recAlignment);
  }
}


void AlignBetweenStrips(DNASequence  &refSeq, 
			DNASequence  &qrySeq,
			T_Alignment *&alignment,
			T_Alignment *&curAlignment,
			T_Strips    *&strips,
			T_Params       params,
			T_ParallelParams &parInfo,
			std::string padding) {

  T_Strips::iterator stripIt, nextIt;
  T_Alignment *recAlignment;
  ssize_t strand;
  ssize_t startRefPos, endRefPos;
  ssize_t startQryPos, endQryPos;
  recAlignment = NULL;
  nextIt = strips->begin();
  ++nextIt;
  ssize_t i;
  for (stripIt = strips->begin(); 
       stripIt != strips->end(); 
       ++stripIt){
    if (nextIt != strips->end() && 
	nextIt->startQryEnum == stripIt->endQryEnum+1) {
      params.AlignGaps();
      strand = Sign(stripIt->endQryEnum);
      // Reset the alignment between strips
      startRefPos = stripIt->endRefPos; //params.hashLength;
      endRefPos   = nextIt->startRefPos;// - 1;
      startQryPos = stripIt->endQryPos; //params.hashLength*strand;
      endQryPos   = nextIt->startQryPos;// - 1*strand;
      std::cout << "Aligning " << startRefPos << " ... " << endRefPos 
		<< " vs " << startQryPos << " ... " << endQryPos << std::endl;
      if (startRefPos < endRefPos) { 
	for (i = startRefPos; i < endRefPos; i++) {
	  refSeq.MarkPosNotCommon(i);
	  alignment->locations[i] = -1;
	  alignment->enumerations[i] = 0;
	}
	recAlignment = AlignSubsequences(refSeq, qrySeq,
					 startRefPos, endRefPos,
					 startQryPos, endQryPos,
					 strand,
					 alignment,
					 params, parInfo, padding);
	AddAlignment(alignment, curAlignment, recAlignment);
      }
    }
    if (nextIt != strips->end())
      ++nextIt;
  }
}

void BuildScaffold(DNASequence &refSeq, DNASequence &qrySeq, 
		   T_Alignment *&alignment,
		   T_Strips *&strips,
		   T_Params params, T_ParallelParams &parInfo,
		   std::string padding) {

  ssize_t i;
  ssize_t *qryEnumerations, *refLocations, *qryLocations;
  alignment = new T_Alignment;
  alignment->hashLength = params.hashLength;
  
  // Find common tuples (may take a while for long sequences).
  ssize_t countEnumerated; // result of how many positions in ref are mapped to qry
  countEnumerated = EnumerateSequences(refSeq, qrySeq, 
				       params, parInfo, alignment);
  
  std::cout << "got " << countEnumerated << " tuples. " << std::endl;
  // If there weren't any common tuples (strange), exit here.
  if (countEnumerated == 0) {
    delete alignment;
    alignment = NULL;
    return;
  }

#ifdef _DEBUG_VERSION
  // validate aligned sequences.
  for (i = 0; i < alignment->length - alignment->hashLength + 1; i++) {
    assert(alignment->locations[i] >= -1);
    if (refSeq.IsACTG(i)) {
      if(refSeq.IsPosMasked(i)) {
	assert(alignment->locations[i] == -1);
	assert(alignment->enumerations[i] == 0);
      }
      else {
	assert(alignment->locations[i] != -1 &&
	       alignment->enumerations[i] != 0);
      }
      if (alignment->enumerations[i] != 0)
	assert(alignment->locations[i] != -1);
    }
  }
#endif

  T_Alignment *recAlignment = NULL;
  T_Alignment *curAlignment = NULL;

  DNASequence refSubSeq, qrySubSeq;
  ssize_t numTuples;
  ssize_t numDiscardedStrips;
  // Advance to the next smaller hash size.
  
  // Transform the enumeration into strips
  // The first step for this is to compress the enumerations array
  // from per nucleotide enumeration to enumeration of locations.
  CreateStripArrays(refSeq,  // description of first seq
		    params, 
		    alignment->enumerations, alignment->locations,  
		    numTuples,
		    qryEnumerations,          //    condensed enumerations
		    refLocations, qryLocations//    Condensed locations
		    );
  strips = new T_Strips;
  GetStrips(qryEnumerations, refLocations, refLocations,
	    qryLocations, qryLocations, NULL,   
	    NULL, numTuples, *strips, params.hashLength, params.alignStrips);
  
  PrintEnumeration( refSeq,
		    params.hashLength,
		    alignment->enumerations,
		    alignment->locations, std::cout);


  
  // This information is now contained in the strips structure
  delete []qryEnumerations;
  delete []refLocations;
  delete []qryLocations;
  /*
   * Go from postions to enumerated strips.
   */

  // Clean up likely spurious strips.  There are ways for cleaning
  // up the strips when aligning strips and aligning gaps.
  // 
  // Aligning strips is done to get a large-grained alignment.  No
  // assumption on ordering and orientation of strips is made. Adjacent
  // common k-mers are merged into strips, and small, out-of-order
  // strips are removed to create larger strips.
  // 
  // Aligning gaps is done to get a finer detailed alignment.  It
  // is assumed that the two sequences are in order, although the
  // set of markers (common k-mers) found when enumerating them
  // may not be (due to low sequence identity, and spurious
  // matches, etc.).  The set of markers is cleaned by finding the
  // longest increasng subset of markers.

  MaximizeStrips(strips, params.lisWindow, params.lisThreshold);
  EnumerateStrips(*strips);
  ReannotateAlignment(refSeq, strips, alignment);
}


void CopyAlignment(DNASequence &seq,
		   T_Alignment *alignment, 
		   ssize_t parentRefOffset, ssize_t parentQryOffset, 
		   ssize_t *locations, ssize_t *enumerations, ssize_t lenLocations,
		   ssize_t *qryLocations, ssize_t lenQryLocations) {
  ssize_t i;
  ssize_t pos, qryPos;
  
  ssize_t positionsCopied = 0;
  ssize_t positionsUnmasked = 0;
  for (i = 0; i < alignment->length - alignment->hashLength + 1; i++) {

    pos = i + alignment->refPos + parentRefOffset;
    assert("alignment->locations overrunning boundary " &&
	   + parentRefOffset < lenLocations);

    if (alignment->enumerations[i] != 0)
      assert(alignment->locations[i] != -1);

    assert(alignment->locations[i] + alignment->qryPos + parentQryOffset < lenQryLocations);
    
    // Look to see if the nucleotide at pos i has been assigned a location in the query genome
    if (alignment->locations[i] != -1) {
      // It is possible that this position on the query sequence has
      // already been mapped to the reference sequence.  In this case,
      // do not map the position in the reference sequence to the
      // already mapped query position.
      // qryLocations contains a list of locations in the query sequence that are mapped to 
      // the reference sequence (a reverseof the alignment->locations array).
      qryPos = alignment->locations[i] + alignment->qryPos + parentQryOffset;

      // Try to give power to the outer alignment.  This position has already 
      // been assigned in the query sequence by a recursive alignment, but a parent
      // aligment also aligned it's location.  
      if (qryLocations[qryPos] == -1 || qryLocations[qryPos] == pos) {
	/*	locations[qryLocations[qryPos]] = -1;
	  enumerations[qryLocations[qryPos]] = 0;
	  qryLocations[qryPos] = 0;
	  seq.MarkPosNotUnique(qryLocations[qryPos]);
	  }*/
	// Position in the query sequence is not yet mapped.
	// Map it here.
	if (!seq.IsPosMasked(pos) ) {
	  // Map the positions in the two sequences.
	  /**********************************************	
	  Experimental code 11/30. 
	    It seems that when anchors are defined, the entire length of the 
	    aligned sequence should be included in the alignment, not just the point
	    where the alignment started. Do that here.
	  */
	  ssize_t a;
	  ssize_t s;
	  s = Sign(alignment->enumerations[i]);
	  for (a = 0; a < alignment->hashLength; a++) {
	    if (s == 1) {
	      if (qryLocations[qryPos+a] == -1 || qryLocations[qryPos+a] == pos+a) {
		locations[pos + a]       = qryPos + a;
		enumerations[pos + a]    = s;
		qryLocations[qryPos + a] = pos + a;
		seq.UnmaskPosNoCk(pos + a); 
	      }
	      else {
		locations[pos + a] = -1; 
		enumerations[pos + a] = 0; 	  
		seq.MarkPosNotUnique(pos + a); 
	      }
	    }
	    else {
	      if (s == -1) {
		if (qryLocations[qryPos+ alignment->hashLength - 1 -a] == -1 || 
		    qryLocations[qryPos+ alignment->hashLength - 1 -a] == pos+a) {
		  locations[pos + a]       = qryPos + alignment->hashLength - 1 - a;
		  enumerations[pos + a]    = s;
		  qryLocations[qryPos + alignment->hashLength - 1 - a] = pos + a; 
		  seq.UnmaskPosNoCk(pos + a); 
		}
		else {
		  locations[pos + a] = -1; 
		  enumerations[pos + a] = 0;
		  seq.MarkPosNotUnique(pos + a); 
		}
	      }
	      else {
		assert(s == 1 || s == -1 || printf("invalid  enumeration %d at %d\n ",s, pos+a) == 0);
	      }
	    }
	    // It was possible that this position was masked by a
	    // recursive alignment.  There is no reason this position
	    // should be removed in the alignment, and so unmask it.  
	    // Masked positions are counted as marked
	  }
	}
      }
      else {
	// Position is already mapped.  This position in  the
	// reference sequence was marked as assigned (unmasked),
	// because in one of the recursive alignment steps did align
	// it.  Since that alignment position is not being used, mark
	// that position as masked.
	if (qryLocations[qryPos] != pos) {
	  locations[pos] = -1;
	  enumerations[pos] = 0;	  
	  seq.MarkPosNotUnique(pos);
	}
      }
    }
    if (i > 0 && pos > 0 
	&& alignment->locations[i-1] != -1 
	&& alignment->locations[i] != -1) {
      assert(locations[pos-1] == -1 ||
	     locations[pos-1] != locations[pos]);
    }
    assert(locations[pos] >= -1);
  }
}

void AssignAlignmentLocations(DNASequence &seq,
			      T_Alignment *alignment, 
			      ssize_t parentRefOffset, 
 			      ssize_t parentQryOffset,
			      ssize_t *locations,
			      ssize_t *enumerations,
			      ssize_t lenLocations, 
			      ssize_t *qryLocations,
			      ssize_t qryLen,
			      ssize_t level) {

  // Alignment specifies an alignment between two sequences by storing
  // in the array alignment->locations all of the locations of the
  // second sequence relative to the first.  This is enriched by the
  // recursive call to assign alignments.

  ssize_t i;
  assert("AssignAlignment overran boundary" 
	 && (parentRefOffset + alignment->length - 1) < lenLocations);
  while (alignment != NULL) {
    // Branch into child if there is one.
    assert(alignment != alignment->next);

    if (alignment->child != NULL)
      AssignAlignmentLocations(seq,
			       alignment->child,
			       parentRefOffset + alignment->refPos,
			       parentQryOffset + alignment->qryPos,
			       locations, enumerations, lenLocations,
			       qryLocations, qryLen, level++);

    CopyAlignment(seq, alignment, parentRefOffset, parentQryOffset, 
		  locations, enumerations, lenLocations,
		  qryLocations, qryLen);    

    alignment = alignment->next;
  }
}

void FreeAlignments(T_Alignment *alignment) {
  T_Alignment *curAlignment;
  while (alignment != NULL) {
    if (alignment->child != NULL) 
      FreeAlignments(alignment->child);
    curAlignment = alignment;
    alignment = alignment->next;
    delete []curAlignment->enumerations;
    delete []curAlignment->locations;
    delete curAlignment;
  }
}	   

void LocationsToEnumeration(ssize_t *locations, ssize_t *enumerations, ssize_t locLength, 
			    ssize_t *enumeration, ssize_t enumLength ) {
  // The locations array map positions in one sequence to
  // their orthologous positions in another sequence.  Each position in
  // the other sequence is assigned an enumeration.  The enumeration
  // is ordered from 1 .. m-u, where m is the length of the query
  // sequence, and u is the number of unmapped positions.
  ssize_t pos;
  ssize_t e;
  pos = 0;
  // Mark the positions that are to be enumerated
  for (pos = 0; pos < locLength; pos++) {
    if (locations[pos] != -1) {
      assert("position mapped outside qry string " && 
	     locations[pos] < enumLength);
      /*      if (pos < (locLength-1))
	      assert(locations[pos] != locations[pos+1]);*/
      //      assert(enumeration[locations[pos]] == 0);
      enumeration[locations[pos]] = enumerations[pos];
    }
  }
  // Enumerate the positions
  e = 1;
  for (pos = 0; pos < enumLength; pos++) {
    if (enumeration[pos] != 0) {
      enumeration[pos] = e * Sign(enumeration[pos]);
      e++;
    }
  }
}

void StoreEnumerationInAlignment(ssize_t *locations,
				 ssize_t *enumerations, //T_Alignment &alignment,
				 ssize_t length,
				 ssize_t *qryEnumerations,
				 ssize_t enumLength) {
  // qryEnumerations stores a list of the enumeration at each position
  // of the qry sequence.  Masked nucleotides are not enumerated.
  
  ssize_t pos;
  for (pos = 0; pos < length; pos++) {
    if (locations[pos] != -1) {
      enumerations[pos] =  qryEnumerations[locations[pos]]; 
    }
  }
}



void ReannotateAlignment(DNASequence &seq,
			 T_Strips *strips,
			 T_Alignment *alignment) {
  
  // Given the longest increasing subset within a strip, all 
  // other annotated positions are spurious and need to be removed from the 
  // stored alignment.  

  // Reset the alignment
  ssize_t i;
  for (i = 0; i < alignment->length; i++) {
    alignment->locations[i] = -1;
    alignment->enumerations[i] = 0;
    if (!seq.IsPosMasked(i))
      seq.MarkPosNotCommon(i);
  }

  // Add all tuples to the alignment
  T_Strips::iterator it, end;
  end = strips->end();
  ssize_t startRefPos, startQryPos;
  for (it = strips->begin(); it != end; ++it) {
    alignment->locations[it->startRefPos] = it->startQryPos;
    alignment->enumerations[it->startRefPos] = it->startQryEnum;
    seq.UnmaskPos(it->startRefPos);
  }

}


void ReadParamFile(std::string paramFileName, 
		   T_Params &parameters, 
		   T_ParallelParams &parallelParams) {
  std::ifstream in;
  in.open(paramFileName.c_str());
  if (!in.good()) {
    std::cout << "could not open param file, bailing out. " << std::endl;
    exit(1);
  }
  
  std::string word;
  while (!in.eof() && in.good()) {
    in >> word;
    // skip this line if it is a comment
    if (word.size() > 0 && word[0] == '#') {
      in.ignore(1000000, '\n');
      continue;
    }
    std::cout << word  << "  ";
    if (word == "wordLength") {
      in >> parameters.wordLength; std::cout << parameters.wordLength << std::endl; }
    else if (word == "neighborDist"){
      in >> parameters.neighborDist; std::cout << parameters.neighborDist << std::endl; }
    else if (word == "lisWindow"){
      in >> parameters.lisWindow; std::cout << parameters.lisWindow << std::endl; }
    else if (word == "lisThreshold"){
      in >> parameters.lisThreshold; std::cout << parameters.lisThreshold << std::endl; }
    else if (word == "distThreshold"){
      in >> parameters.distThreshold; std::cout << parameters.distThreshold << std::endl; }
    else if (word == "mergeThreshold"){
      in >> parameters.mergeThreshold; std::cout << parameters.mergeThreshold << std::endl; }
    else if (word == "productSize"){
      in >> parameters.productSize; std::cout << parameters.productSize << std::endl; }
    else if (word == "hashLength"){
      in >> parameters.hashLength; std::cout << parameters.hashLength << std::endl; }
    else if (word == "minHashLength"){
      in >> parameters.minHashLen; std::cout << parameters.minHashLen << std::endl; }
    else if (word == "hashStep"){
      in >> parameters.hashStep; std::cout << parameters.hashStep << std::endl; }
    else if (word == "outputTuple"){
      in >> parameters.outputTuple; std::cout << parameters.outputTuple << std::endl; }
    else if (word == "outputGlue") {
      in >> parameters.outputGlue; std::cout << parameters.outputGlue << std::endl; }
    else if (word == "maskRepeats"){
      in >> parameters.maskRepeats; std::cout << parameters.maskRepeats << std::endl; }
    else if (word == "numFragments"){
      in >> parallelParams.numFragments; std::cout << parallelParams.numFragments << std::endl; }
    else if (word == "numChildren"){
      in >> parallelParams.numChildren; std::cout << parallelParams.numChildren << std::endl; }
    else if (word == "minGapDist"){
      in >> parameters.minGapDist; std::cout << parameters.minGapDist << std::endl; }
    else if (word == "maxGapDist"){
      in >> parameters.maxGapDist; std::cout << parameters.maxGapDist << std::endl; }
    else if (word == "ratioThreshold"){
      in >> parameters.ratioThreshold; std::cout << parameters.ratioThreshold << std::endl; }
    else if (word == "alignStrips") {
      parameters.alignStrips = 1;
      parameters.alignGaps   = 0;
    }
    else if (word == "alignGaps") {
      parameters.alignGaps = 1;
      parameters.alignStrips = 0;
    }
    else {
      std::cout << "Unrecognized option '" << word << "' in parameter file" << std::endl;
      exit(1);
    }
    //    in.getline();
  }
  in.close();
}

ssize_t factorial(ssize_t k) {
  ssize_t i, res;
  res = 1;
  for (i = 1; i < k; i++) {
    res = res * i;
  }
  return res;
}

ssize_t ExtendMatch(DNASequence &refSeq, DNASequence &qrySeq, 
		ssize_t refPos, ssize_t qryPos, 
		ssize_t dir, double perror) {
  // Extend a match until the frequency of errors is high enough that
  // the are occuring with frequency much higher than perror.

  ssize_t r, q, t;
  r = refPos;
  q = qryPos;
  t = 0;
  ssize_t mma, mmb; // mismatch positions
  ssize_t nerrors;
  double walkprob;
  nerrors = 0;
  ssize_t extend;
  extend = 1;
  mma = -1;
  mmb = -1;
  ssize_t eq;
  while ( r < refSeq.length && q < qrySeq.length && extend) {
    if (dir > 0) 
      eq = (DNASequence::UnmaskedValue(refSeq.seq[r]) != 
	    DNASequence::UnmaskedValue(qrySeq.seq[q])) ? 1 : 0;
    else
      eq  = (DNASequence::UnmaskedValue(refSeq.seq[r]) != 
	     NucRC(DNASequence::UnmaskedValue(qrySeq.seq[q]))) ? 1 : 0;

    if (eq) {
      nerrors++;
      // Deterine if the number of errors seen so far is too many.
      ssize_t k = nerrors - 1;
      walkprob = std::pow(perror*t, k)*exp(-perror*t)/factorial(k);
      if (walkprob < 0.05) {
	// Finished with this extension.
	extend = 0;
      }
    }
    r += dir;
    q += dir;
    t++;
  }
  return t;
}
