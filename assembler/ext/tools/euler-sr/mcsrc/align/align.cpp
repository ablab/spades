/***************************************************************************
 * Title:          align.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <stdlib.h>
#include <string>
#include <unistd.h>

#include "DNASequence.h"
#include "SeqReader.h"
//#include "EnumUtils.h"
#include "AlignmentPrinter.h"
#include "Alignment.h"
#include "align/alignutils.h"

void PrintUsage();
void InitEnv(int argc, char *argv[], 
	     std::string &inputFileA,
	     std::string &inputFileB,
	     ssize_t &mismatch,
	     ssize_t &gapOpen,
	     ssize_t &gapExtend,
	     ssize_t &nonAffineGap,
	     ssize_t &match,
	     ssize_t &kband,
	     std::string &scoreMatFile,
	     std::string &scoreMatStr,
						 std::string &alignType,
						 std::string &outputType);

int main(int argc, char* argv[]) {
  
  DNASequence *seq1;
  DNASequence *seq2;
  ssize_t kband;
  std::string alignType;
  ssize_t mismatch, gapOpen, gapExtend, nonAffineGap, match;
  match        = -100;
  mismatch     = 200;
  nonAffineGap = 800;
  gapOpen      = 400;
  gapExtend    = 30;
  kband        = 2;
  
  std::string scoreMatFileName, scoreMatString;
  IntMatrix scoremat;
  std::string inA, inB;
  alignType = "affine";
	std::string outputType = "full";
  inA = "";
  inB = "";

  scoreMatFileName = "";
  scoreMatString   = "";
  InitEnv(argc, argv, 
					inA, inB, 
					mismatch, gapOpen, gapExtend, nonAffineGap, match,
					kband,  scoreMatFileName, scoreMatString, alignType, outputType);

  if (inA == "" || inB == "") {
    PrintUsage();
    exit(1);
  }

  if (scoreMatFileName != "") {
    ReadScoreMatFile(scoreMatFileName, scoremat);
  }
  else 
    InitScoreMat(scoremat, match, mismatch);

  std::ifstream ifa, ifb;

  ifa.open(inA.c_str());
  if (! ifa.good() ) { std::cout << "could not open " << inA << std::endl; exit(0);}
  
  ifb.open(inB.c_str());
  if (! ifb.good() ) { std::cout << "could not open " << inB << std::endl; exit(0);}
  
  SeqReader readerA(&ifa);
  SeqReader readerB(&ifb);

  double score = 0.0;
  ssize_t length, refAlignStart, qryAlignStart;
  ssize_t *qryLocations, *refLocations, *enumerations, *locations;
  double *optScores;
    
  while(1) {
    if (! (readerA.GetSeq(seq1, SeqReader::noConvert)))
      break;
    if (! (readerB.GetSeq(seq2, SeqReader::noConvert)))
      break;

    if (seq1 == NULL || seq2 == NULL) {
      std::cout << "Error getting sequences from " << argv[1] << std::endl;
    }
    locations    = new ssize_t[seq1->length + 1];
    optScores    = new double[seq1->length + 1];
    qryLocations = new ssize_t[seq1->length + 1];
    refLocations = new ssize_t[seq1->length + 1];
    enumerations = new ssize_t[seq1->length + 1];
    IntMatrix pathMat;
    IntMatrix scoreMat;
    IntMatrix mismatchMat;
    InitScoreMat(mismatchMat, 0, 1);
    ssize_t i;
    for (i = 0; i < seq1->length; i++) {
      qryLocations[i] = refLocations[i] = -1;
      enumerations[i] = 0;
    }

    for (i = 0; i < seq1->length; i++) {
      locations[i] = -1;
      optScores[i] = INF;
    }
    if (alignType == "affine")
      score = AffineAlign(*seq1, *seq2, 
			  match, mismatch, nonAffineGap, gapOpen, gapExtend, 
			  locations, scoremat);
    else if (alignType == "band")
      score = BandedAlign(*seq1, *seq2, 
													match, mismatch, nonAffineGap, kband, 
													locations, scoreMat, pathMat, mismatchMat);
    else if (alignType == "fit") {
			IntMatrix scoreCache;
			IntMatrix pathCache;
			CreateMatrix(scoreCache, seq1->length+1, seq2->length+1);
			CreateMatrix(pathCache,  seq1->length+1, seq2->length+1);
      score = FitAlign(*seq1, *seq2, 
											 match, mismatch, nonAffineGap, 
											 locations, scoremat, scoreCache, pathCache);
		}
    else if (alignType == "overlap")
      score = OverlapAlign(*seq1, *seq2, 
			   match, mismatch, nonAffineGap,
			   locations, scoremat);
    else if (alignType == "local")
      score = LocalAlign(*seq1, *seq2, 
			 match, mismatch, nonAffineGap, 
			 locations, scoremat);
    else if (alignType == "affloc") 
      score = AffineLocalAlign(*seq1, *seq2, 
			       match, mismatch, nonAffineGap, gapOpen, gapExtend, 
			       locations, scoremat);
    else if (alignType == "nw") 
      score  = Align(*seq1, *seq2, 
		     match, mismatch, nonAffineGap, 
		     locations, scoremat);
    else if (alignType == "graph") {
      score = ScoreBandedAffineAlign(*seq1, *seq2, 
				     1500, scoremat, gapOpen, gapExtend, 3000, 
				     locations, length, 
				     refAlignStart, qryAlignStart, 0, 0);

  
      ssize_t e;
      e = 0;
      for (i = 0; i < length; i++) {
	if (locations[i] != -1) {
	  enumerations[e] = e+1;
	  refLocations[e] = i;
	  qryLocations[e] = locations[i];
	  //      std::cout << i << " " << locations[i] << std::endl;
	  e++;
	}
      }

			if (outputType == "full") {
				PrintAlignment(*seq1, *seq2, 0, 0, locations, length + 1, std::cout);
			}
			else if (outputType == "summary") {
				ssize_t nMatch, nMisMatch, nGapQry, nGapRef;
				nMatch = nMisMatch = nGapQry = nGapRef = 0;
				ComputeAlignmentStatistics(locations, length + 1,
																	 (char*) seq1->seq, (char*) seq2->seq,
																	 nMatch, nMisMatch, nGapQry, nGapRef);
				std::cout << nMatch<< " " << nMisMatch << " " << nGapQry << " " << nGapRef << std::endl;
			}
			//      PrintAlignment(&std::cout, *seq1, *seq2, enumerations, refLocations, qryLocations, e);
      return 0;
    }  

    ssize_t e;
    e = 0;
    for (i = 0; i < seq1->length; i++) {
      if (locations[i] != -1) {
	enumerations[e] = e+1;
	refLocations[e] = i;
	qryLocations[e] = locations[i];
	e++;
      }
    }
			if (outputType == "full") {
				std::cout << seq1->namestr << " " << seq1->length << std::endl;
				std::cout << seq2->namestr << " " << seq2->length << std::endl;
				std::cout << score << std::endl;
		
				PrintAlignment(*seq1, *seq2, 0, 0, locations, seq1->length, std::cout);
			}
			else if (outputType == "summary") {
				ssize_t nMatch, nMisMatch, nGapQry, nGapRef;
				nMatch = nMisMatch = nGapQry = nGapRef = 0;
				ComputeAlignmentStatistics(locations, length + 1,
																	 (char*) seq1->seq, (char*) seq2->seq,
																	 nMatch, nMisMatch, nGapQry, nGapRef);
				std::cout << nMatch<< " " << nMisMatch << " " << nGapQry << " " << nGapRef << std::endl;
			}

    /*    PrintAlignment(&std::cout, *seq1, *seq2, 
		   enumerations, refLocations, qryLocations, e);
    */
    delete [] enumerations;
    delete [] refLocations;
    delete [] qryLocations;
    delete [] locations;
    delete [] optScores;
    seq1->Reset();
    seq2->Reset();
    delete seq1;
    delete seq2;
    seq1 = NULL;
    seq2 = NULL;
  }
  if (seq1 != NULL) {
    seq1->Reset();
    delete seq1;
  }
  if (seq2 != NULL) {
    seq2->Reset();
    delete seq2;
  }
  ifa.close(); ifb.close();
}

void InitEnv(int argc, char *argv[], 
	     std::string &inputFileA,
	     std::string &inputFileB,
	     ssize_t &mismatch, 
	     ssize_t &gapOpen,
	     ssize_t &gapExtend,
	     ssize_t &nonAffineGap,
	     ssize_t &match,
	     ssize_t &kband,
	     std::string &scoreMatFile,
	     std::string &scoreMatStr,
						 std::string &alignType,
						 std::string &outputType) {
  ssize_t copt;
  while ( (copt=getopt(argc, argv, "s:S:m:M:o:g:e:ht:k:O:")) != EOF) {
    switch(copt) {
    case 's':
      scoreMatFile = optarg;
      continue;
    case 'S':
      scoreMatStr = optarg;
      continue;
    case 'm':
      match = atoi(optarg);
      continue;
    case 'M':
      mismatch = atoi(optarg);
      continue;
    case 'o':
      gapOpen = atoi(optarg);
      continue;
    case 'e':
      gapExtend = atoi(optarg);
      continue;
    case 'g':
      nonAffineGap = atoi(optarg);
      continue;
    case 'k':
      kband = atoi(optarg);
      continue;
    case 't':
      alignType = optarg;
      continue;
		case 'O':
			outputType = optarg;
			continue;
    case 'h':
      PrintUsage();
      exit(0);
    default:
      std::cout << "error with option: " << copt << " arg: " << optarg << std::endl;
      PrintUsage();
      exit(1);
    }
  }
  ssize_t i;
  i = optind;
  if (i < argc) {
    inputFileA = argv[i];
  }
  i++;
  if (i < argc) {
    inputFileB = argv[i];
  }
}

void PrintUsage() {
  std::cout << "Align two sequences from a file\n" << std::endl;
  std::cout << "Find MINImum alignment score of two sequences " << std::endl;
  std::cout << "usage:  align [-m match -M mismatch -o gapOpen -e gapExtend -g nonAffineGap] refseq qryseq" << std::endl;
  std::cout << "Alignment corresponds to minimum penalty path (lowest scoring). " << std::endl;
  std::cout << "   -t alignType , where aligntype can be : "<< std::endl;
  std::cout << "                  nw     - Needleman Wunsch (unoptimized) " << std::endl;
  std::cout << "                  affine - for affine alignment " << std::endl;
  std::cout << "                  band   - for banded alignment " << std::endl;
  std::cout << "                  fit    - for fitting alignment where full query sequence is fit into " << std::endl
	    << "                           reference sequence " << std::endl;
  std::cout << "                  local  - local S.W. alignment " << std::endl;
  std::cout << "                  affloc - local-affine alignment " << std::endl;
  std::cout << "                  graph - dynamic graph based alignment.  " << std::endl;
  std::cout << "                          explores high scoring alignment space of two sequences." << std::endl;
  std::cout << "   -m match (-1.0) nucleotide match penalty " << std::endl;
  std::cout << "   -M mismatch (1.0) nucleotide mismatch penalty " << std::endl;
  std::cout << "   -o gapOpen (6.0) penalty for opening a gap. " << std::endl;
  std::cout << "   -e gapExtend (0.0) penalty for extending a gap. " << std::endl;
  std::cout << "   -g gap (1.0) penalty for nonaffine gap. " << std::endl;
  std::cout << "   -s scoreMat scoring matrix file " << std::endl;
  std::cout << "   -S \"scorematstring \" a string with 16 numbers in order a c g t in rows & cols " << std::endl;
}
