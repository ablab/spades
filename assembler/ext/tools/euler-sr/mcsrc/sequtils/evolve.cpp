/***************************************************************************
 * Title:          evolve.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "SeqUtils.h"
#include "compatibility.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>

void PrintUsage();
void InitEnv(int argc, char *argv[], 
	     std::string &refSeqFile,
	     double &divergence,
	     double &fractionIndels,
	     double &rearrangementFraction,
	     ssize_t &numInversions,
	     ssize_t &numTranslocations,
	     ssize_t &randSeed,
	     std::string &outfile);

ssize_t GetRandom(ssize_t max) {
  long rval = random();
  double ratio = double(rval)/RANDOM_MAX;
  return (ssize_t) floor(ratio*max);
}


double GetRandomFraction() {
  long rval = random();
  return double(rval)/RANDOM_MAX;
}  

std::string RearrangeGenome(DNASequence &genome, ssize_t numInversions, ssize_t numTranslocations, double rearrangementFraction) {

  ssize_t i;
  ssize_t type;
  ssize_t length, sourceStart, sourceEnd;
  ssize_t destStart, destEnd;
  std::stringstream sstr;
  ssize_t p;
  for (i = 0; i < numInversions; i++) {
    std::cout << "inverting " << std::endl;
    length   = (ssize_t) ceil(rearrangementFraction * genome.length * GetRandomFraction());
    sourceStart = GetRandom(genome.length - length);
    sourceEnd = sourceStart + length - 1;
    sstr << "("<<sourceStart << ".-." << sourceEnd << ")";
    destStart = sourceEnd;
    destEnd = sourceStart;

    // Reverse the sequence
    unsigned char tmp;
    for (p = 0; p < ceil(length/2.0); p++) {
      tmp = genome.seq[sourceStart + p];
      genome.seq[sourceStart + p] = comp_bin[genome.seq[sourceEnd - p]];
      genome.seq[sourceEnd - p] = comp_bin[tmp];
    }
    if (numTranslocations > 0) {
      type = GetRandom(2);
      std::cout << "type: " << type << std::endl;
      if (type == 1) {
	std::cout << "translocating inversion " << std::endl;
	sstr << "("<<sourceStart << "("<< length << ") ->" << destStart << ")";
	numTranslocations--;
	destStart = GetRandom(genome.length - length);
	destEnd = destStart + length - 1;

	DNASequence tempSeq(length);
	for (p = 0; p < length; p++) 
	  tempSeq.seq[p] = genome.seq[sourceStart + p];
	
	if (destStart > sourceStart) {
	  // Fragment will be inserted after original location
	  for (p = 0; p < (destStart - sourceStart+1); p++) 
	    genome.seq[sourceStart + p] = genome.seq[sourceStart + p + length];
	}
	else {
	  // Fragment will be inserted before original location, make room before 
	  // the fragment
	  for (p = sourceStart - destStart; p >= 0; p--)
	    genome.seq[destStart + p + length] = genome.seq[destStart + p];
	}
	// Copy the translocated sequence
	for (p = 0; p < length; p++)
	  genome.seq[destStart + p] = tempSeq[p];
      }
    }
  }
  // Translocate sequence
  for (i = 0; i < numTranslocations; i++) {
    std::cout << "translocating " << std::endl;
    length   = (ssize_t) ceil(rearrangementFraction * genome.length * GetRandomFraction());
    sourceStart = GetRandom(genome.length - length);
    sourceEnd = sourceStart + length - 1;
    destStart = GetRandom(genome.length - length);
    destEnd = destStart + length - 1;
    sstr << "("<<sourceStart << "("<< length << ") ->" << destStart << ")";
    std::cout << sstr.str() << std::endl;
    DNASequence tempSeq(length);
    for (p = 0; p < length; p++) 
      tempSeq.seq[p] = genome.seq[sourceStart + p];
    
    if (destStart > sourceStart) {
      // Fragment will be inserted after original location
      std::cout << "to after source " << std::endl;
      for (p = 0; p < (destStart - sourceStart+1); p++) {
	genome.seq[sourceStart + p] = genome.seq[sourceStart + p + length];
      }
    }
    else {
	  // Fragment will be inserted before original location, make room before 
	  // the fragment
      std::cout << "moving to before source " << std::endl;
      for (p = (sourceStart - destStart); p >= 0; p--)
	genome.seq[destStart + p + length] = genome.seq[destStart + p];
    }
    // Copy the translocated sequence
    for (p = 0; p < length; p++) {
      genome.seq[destStart + p] = tempSeq[p];
    }
  }
  return sstr.str();
}

int main(int argc, char* argv[]) {
  double divergence, fractionIndels, rearrangementFraction;
  ssize_t numInversions, numTranslocations;
  ssize_t randSeed;
  std::string outfile, refSeqFile;
  divergence = 0.05;
  fractionIndels = 0.02;
  numInversions = 0;
  numTranslocations = 0;
  randSeed = 0;
  rearrangementFraction = 0.25;
  InitEnv(argc, argv, refSeqFile,
	  divergence, fractionIndels,rearrangementFraction,
	  numInversions, numTranslocations, 
	  randSeed, outfile);

  // Initialize random number generation
  time_t t;
  if (randSeed == 0) 
    randSeed = (ssize_t) time(&t);
  srandom((_UINT_) randSeed);

  // Read the reference sequence
  std::ifstream in;
  in.open(refSeqFile.c_str());
  if (! in.good() ) {
    std::cout << "Could not open " << refSeqFile << std::endl;
    exit(1);
  }
  SeqReader reader(&in);
  DNASequence refSeq;
  reader.GetSeq(refSeq);
  in.close();
  // Go about making modifications to the sequence.

  ssize_t numEvents, numMutations, numIndels;

  numEvents = (ssize_t) floor(refSeq.length * divergence);
  // muations done * 3/4 = actual rate, because of 1/4 odds of mutating
  // back to the same nuc.
  // so actual rate = 4/3 * mutations done
  
  numMutations = (ssize_t) floor((1 - fractionIndels) * numEvents * 4/3); 
  numIndels    = (ssize_t) floor(fractionIndels * numEvents);

  // Create the copied sequence
  DNASequence workingSeq(refSeq.length);
  workingSeq._ascii = 0;
  ssize_t *indels = new ssize_t[refSeq.length];
  memcpy(workingSeq.seq, refSeq.seq, refSeq.length);
  memset(indels, 0, sizeof(ssize_t) * refSeq.length);


  ssize_t i;
  // Evolve the sequence
  ssize_t pos;
  ssize_t nuc;
  for (i = 0; i < numMutations; i++) {
    pos = GetRandom(refSeq.length);
    nuc = GetRandom(4);
    workingSeq.seq[pos] = nuc;
  }

  // Create indels
  ssize_t type;
  //UNUSED// ssize_t lengthChange = 0;
  ssize_t numDels, numIns;
  numDels = 0; numIns = 0;
  for (i = 0; i < numIndels; i++) {
    pos = GetRandom(refSeq.length);
    // Distribute evenly between indels
    type = GetRandom(2);
    if (type == 0) {
      indels[pos]--;
      numDels++;
    }
    else {
      indels[pos]++;
      numIns++;
    }
  }
  

  ssize_t w = 0;
  ssize_t j;
  ssize_t newLength = 0;

  // Calculate the length
  for (w = 0; w < workingSeq.length;) {
    if (indels[w] < 0) {
      // indel here.  Don't add anything to evolvedSeq, and move to next nucleotide
      w+= -indels[w];
    }
    if (indels[w] >= 0) {
      for (j = 0; j < indels[w]; j++) {
	newLength++;
      }
      newLength++;
      w++;
    }
  }

  DNASequence evolvedSeq(newLength);
  evolvedSeq._ascii = 0;
  // Calculate the length
  ssize_t e = 0;
  for (w = 0; w < workingSeq.length;) {
    if (indels[w] < 0) {
      // indel here.  Don't add anything to evolvedSeq, and move to next nucleotide
      w+= -indels[w];
    }
    if (indels[w] >= 0) {
      for (j = 0; j < indels[w]; j++) {
	evolvedSeq.seq[e] = workingSeq[w];
	e++;
      }
      evolvedSeq.seq[e] = workingSeq[w];
      e++;
      w++;
    }
  }

  std::ofstream outputfile;
  std::ostream *output;
  if (outfile == "")
    output = &std::cout;
  else {
    outputfile.open(outfile.c_str());
    if (!outputfile.good()) {
      std::cout << "Error opening output file " << outfile << std::endl;
      exit(0);
    }
    output = &outputfile;
  }
  std::string rstr;
  rstr = RearrangeGenome(evolvedSeq, numInversions, numTranslocations, rearrangementFraction);
  std::stringstream titleStream;
  titleStream << refSeq.namestr << ".evolved " << divergence << " ("<<fractionIndels << ") " << rstr;
  evolvedSeq.StoreName((char*)titleStream.str().c_str());
  evolvedSeq.PrintSeq(*output);

  if (outfile != "")
    outputfile.close();
    
  return 0;
}



void InitEnv(int argc, char *argv[], 
	     std::string &refSeqFile,
	     double &divergence,
	     double &fractionIndels,
	     double &rearrangementFraction,
	     ssize_t &numInversions,
	     ssize_t &numTranslocations,
	     ssize_t &randSeed,
	     std::string &outfile) {
  ssize_t copt;
  while ( (copt=getopt(argc, argv, "i:d:f:o:v:t:r:F:")) != EOF) {
    switch(copt) {
    case 'i':
      refSeqFile = optarg;
      continue;
    case 'o':
      outfile = optarg;
      continue;
    case 'd':
      divergence = atof(optarg);
      continue;
    case 'f':
      fractionIndels = atof(optarg);
      continue;
    case 'r':
      randSeed  = atoi(optarg);
      continue;
    case 'v':
      numInversions = atoi(optarg);
      continue;
    case 't':
      numTranslocations = atoi(optarg);
      continue;
    case 'F':
      rearrangementFraction = atof(optarg);
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  if (refSeqFile == "") {
    PrintUsage();
    exit(1);
  }
}

void PrintUsage() {
  std::cout << "Create a sequence that has evolved from another one.  Both small scale (indels & mutations) " << std::endl;
  std::cout << "and large scale (inversions) happen in the alignment. "<< std::endl;
  std::cout << "Eventually, insertions of canonical repeats and deletions will happen. " << std::endl;
  std::cout << "Also, there should be a map file that has the proper alignment of nucleotides (output in the strip format). " << std::endl;
  std::cout << "usage: evolve -i inputseq -o outputseq -d divergence -f fractionIndels -v numInversions -t numTranslocations " << std::endl;
  std::cout << "   -i inputSeq, the reference sequence " << std::endl;
  std::cout << "   -o outputseq, the evolved sequence " << std::endl;
  std::cout << "   -d divergence, sequence divergence, as a percentage " << std::endl;
  std::cout << "   -f fractionIndels,  fraction of all modifications made in sequence divergence that are indels " << std::endl;
  std::cout << "   -v numInversions, number of inversions.  I'm not sure what size distribution to use yet. " << std::endl;
  std::cout << "   -t numTranslocations, number of translocations.  ditto for the size dist on this. " << std::endl;
  std::cout << "   -F reversalFraction, amount of genome to consider in reversals. " << std::endl;
}

