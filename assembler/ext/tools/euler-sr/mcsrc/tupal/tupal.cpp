/***************************************************************************
 * Title:          tupal.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "TupalLib.h"

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

ssize_t counter;
void printusage();
void initenv(int argc, char *argv[], 
	     T_Params &params, 
	     T_ParallelParams &parInfo,
	     std::vector<std::string> & inpfiles, 
	     std::string &outfile);



std::ostream *ostr;
int main(int argc, char* argv[]) {
  ssize_t mytid;                  /* my task id */
  ssize_t n, nproc, numt, who, msgtype, nhost, narch;

  // always need an i and a j
  ssize_t i, j, s;
  ssize_t g, f;
  counter = 0;
  std::vector<std::string> inpfiles;
  std::string outfile;
  std::ifstream in;
  std::vector<ssize_t> len_seq_vect;
  DNASequence *fragment;
  ostr = &std::cout;
  struct pvmhostinfo *hostp;
    
  ssize_t spawnResult;

  T_Params mainParams; // comes with default parameters
  T_ParallelParams parInfo;
  
  /*
    Initial values
    hashLength   = 20;   // number of nucleotides to hash on
    wordLength   = 14;   // hash table size = 4^word length
    hashDecrease = 2;

    // Search options
    dist         = 1;    // Edit distance, remove all tuples within 
                         // dist of a tuple.
    doIndel      = 0;    // Don't search for insertions and deletions
                         // of tuples (takes time).
  */

  // Initially retain all tuples
  initenv(argc,argv, mainParams, parInfo,
	  inpfiles, outfile);

  if (outfile == "") {
    std::cout << "no outfile specified " << std::endl;
    printusage();
    if (parInfo.numChildren > 1) 
      pvm_exit();
    exit(0);
  }

  if (parInfo.numChildren > 1) {
    parInfo.mytid = pvm_mytid();
    if (mytid < 0) {
      std::cout << "pvm error " << std::endl;
      pvm_perror(argv[0]);
      return 1;
    }
    char hostname[200];
    gethostname(hostname,200);
    parInfo.childIds = new ssize_t[parInfo.numChildren];
    spawnResult = pvm_spawn("task",
			    (char**)NULL, PvmTaskDefault,    
			    NULL, parInfo.numChildren, parInfo.childIds);  
    if (spawnResult != parInfo.numChildren) {
      std::cout << "error spawning tasks " << std::endl;
      return 1;
    }
  }

  // Step 0.  Read in the genomes.

  std::vector<std::string>::iterator inpit;
  T_Sequences *forSequencePtr;
  T_GenomeVect genomes;
    
  DNASequence *seq;
  T_Alignment *alignment;
  DNASequence *refSeq;
  DNASequence *qrySeq;
  DNASequence xptRefSeq;
  DNASequence xptQrySeq;
  ssize_t genomeLength = 0;
  if (mainParams.loadCheckpoint) {
    alignment = new T_Alignment;
    ReadCheckpoint(*alignment, mainParams, xptRefSeq, xptQrySeq, mainParams.checkpointFile);
    refSeq = &xptRefSeq;
    qrySeq = &xptQrySeq;
    mainParams.loadCheckpoint = 1;
    mainParams.doCheckpoint = 0;
  }
  else {
    for (inpit = inpfiles.begin(); inpit != inpfiles.end(); ++inpit) {
      // Read in each genome, possibly multiple contigs.
    std::cout << "reading " << *inpit << std::endl;
    in.open(inpit->c_str());
    if (! in.good() ) {
      std::cout << "could not open " << inpit->c_str() << std::endl;
      exit(0);
    }

    forSequencePtr =  new T_Sequences;
    
    SeqReader reader(&in);
    if (mainParams.maskRepeats)
      reader.MaskRepeats();

    // Read in all sequences in the file. Currently only the first
    // sequences are aligned.  They may be concatenated in the future,
    // but that's not for sure.
    while (reader.GetSeq(seq)) { 
      forSequencePtr->push_back(seq);
    }
    // Store this genome.
    genomes.push_back(forSequencePtr);
    in.close();
    in.clear();
    }
    std::cout <<"Done reading, start aligning " << std::endl;
    // Arrays that hold the enmueration and location of each nucleotide, 
    // and have a position masked where there is no enumeration for the
    // nucleotide.
    ssize_t dist = 1;

    if (genomes.size() < 2) {
      std::cout << "you must specify 2 sequences to align " << std::endl;
      exit(0);
    }
    refSeq = genomes[0][0][0];
    qrySeq = genomes[1][0][0];

  }
  T_Strips *scaffold;
  BuildScaffold(*refSeq, *qrySeq, alignment, scaffold,
		mainParams, parInfo, "");
  std::cout << "scaffold: " << std::endl;
  PrintEnumerations(*scaffold, std::cout);
  mainParams.hashLength -= mainParams.hashStep;
  // Now given the scaffold, perform recursive alignment
  T_Alignment *curAlignment;
  //  AlignStrips(*refSeq, *qrySeq, alignment, curAlignment, scaffold, mainParams, parInfo, "");
  AlignBetweenStrips(*refSeq, *qrySeq, alignment, curAlignment, scaffold, mainParams, parInfo, "");

  ssize_t *qrySeqArray;
  qrySeqArray = new ssize_t[qrySeq->length];
  // qrySeqArray is first used to hold locations in the 
  // query sequence that are assigned positions in the reference 
  // sequence.  When there are multiple positions in the reference
  // sequene that map to the query sequence, this array will resolve
  // conflicts in order of the first position mapped.
  for (i = 0; i < qrySeq->length; i++) 
    qrySeqArray[i] = -1;
  // 0 = no offset
  
  ssize_t *locations, *enumerations;
  locations  = new ssize_t[alignment->length];
  enumerations = new ssize_t[alignment->length];
  for (i = 0; i < alignment->length; i++) {
    locations[i] = -1;
    enumerations[i] = 0;
  }
  if (alignment != NULL) {
    AssignAlignmentLocations(*refSeq,
			     alignment, 0, 0, 
			     locations, //			     alignment->locations, 
			     enumerations, //			     alignment->enumerations,
			     alignment->length, 
			     qrySeqArray, qrySeq->length, 
			     0);
    /*    PrintLocations(alignment->locations, alignment->length,
	  std::cout);*/

    // Now qrySeqArray is used for holding enumerations, reset 
    // to an un-enumerated value, 0.
    for (i = 0; i < qrySeq->length; i++) 
      qrySeqArray[i] = 0;

    LocationsToEnumeration(locations, 
			   enumerations,
			   //alignment->locations, 
			   //			   alignment->enumerations, 
			   alignment->length, 
			   qrySeqArray, qrySeq->length);

    StoreEnumerationInAlignment(locations, //*alignment,
				enumerations,
				alignment->length,
				qrySeqArray,
				qrySeq->length);

    delete []qrySeqArray;

    std::ofstream out;
    out.open(outfile.c_str());
    if (!out.good()) {
      std::cout << "could not open " << outfile << std::endl;
      exit(0);
    }
    
    if (mainParams.outputGlue) {
      PrintGlue(alignment->locations, alignment->length, out);
    }
    else {
      PrintEnumeration(*refSeq, mainParams.hashLength, 
		       enumerations, //		       alignment->enumerations, 
		       locations, //alignment->locations,
		       out);
    }
    out.close();
  }

  // a new nodification.
  // Step 4.  Print the result.

  
  // Shut down children and exit pvm
  if (parInfo.numChildren > 1) {
    for (i = 0; i < parInfo.numChildren; i++)
      CommunicationLib::InitiateTask(parInfo.childIds[i], doExit);
    pvm_exit();
  }
  for (g = 0; g < genomes.size(); g++) {
    forSequencePtr = genomes[g];
    for (i = 0; i < forSequencePtr->size(); i++) {
      if ((*forSequencePtr)[i]->seq != NULL)
	delete [] (*forSequencePtr)[i]->seq;
      if ((*forSequencePtr)[i]->name != NULL)
	delete [] (*forSequencePtr)[i]->name;
      delete (*forSequencePtr)[i];
    }
    delete forSequencePtr; 
  }

  delete []enumerations;
  delete []locations;
  FreeAlignments(alignment);
}

void initenv(int argc, char *argv[], T_Params &params, T_ParallelParams &parInfo,
	     std::vector<std::string> & inpfiles, 
	     std::string &outfile){
  ssize_t copt;
  std::string inpfile;
  while ( (copt=getopt(argc, argv, "i:o:w:n:p:c:f:Itms:g:G:d:a:l:uP:S:W:T:M:r:R:k:K:")) != EOF) {
    switch(copt) {
    case 'a':
      if (strcmp(optarg, "s") == 0) {
	params.alignStrips = 1;
	params.alignGaps   = 0;
      }
      else if (strcmp(optarg, "g") == 0) {
	params.alignGaps   = 1;
	params.alignStrips = 0;
      }
      else {
	printusage();
	exit(0);
      }
      continue;
    case 'i':
      inpfile = optarg;
      inpfiles.push_back(inpfile);
      continue;
    case 'I':
      params.doIndel = 1;
      continue;
    case 'k':
      params.checkpointFile = optarg;
      params.doCheckpoint = 1;
      continue;
    case 'K':
      params.checkpointFile = optarg;
      params.loadCheckpoint = 1;
      continue;
    case 'o':
      outfile = optarg;
      continue;
    case 'R':
      params.ratioThreshold = atof(optarg);
      continue;
    case 'W':
      sscanf(optarg, "%d", &params.lisWindow);
      continue;
    case 'T':
      sscanf(optarg, "%d", &params.lisThreshold);
      continue;
    case 'd':
      sscanf(optarg, "%d", &params.distThreshold);
      continue;
    case 'M':
      sscanf(optarg, "%d", &params.mergeThreshold);
      continue;
    case 'w':
      sscanf(optarg,"%d", &params.wordLength);
      continue;
    case 'P':
      params.productSize = atoi(optarg);
      continue;
    case 'p':
      sscanf(optarg,"%d", &params.hashLength);
      continue;
    case 's':
      sscanf(optarg, "%d", &params.hashStep);
    case 't':
      params.outputTuple = 1;
      continue;
    case 'u':
      params.outputGlue  = 1;
      continue;
    case 'm':
      params.maskRepeats = 1;
      continue;
    case 'l':
      sscanf(optarg, "%d", &params.minHashLen);
      continue;
    case 'n':
      sscanf(optarg, "%d", &params.neighborDist);
      continue;
    case 'f':
      sscanf(optarg, "%d", &parInfo.numFragments);
      continue;
    case 'c':
      sscanf(optarg, "%d", &parInfo.numChildren);
      continue;
    case 'g':
      sscanf(optarg, "%d", &params.minGapDist);
      continue;
    case 'G':
      sscanf(optarg, "%d", &params.maxGapDist);
      continue;
    case 'r':
      ReadParamFile(optarg, params, parInfo);
      continue;
    default:
      printusage();
      exit(1);
    }
  }
  if (inpfiles.size() == 0) {
    printusage();
    exit(1);
  }
  if (params.hashLength < params.minHashLen)
    params.minHashLen = params.hashLength;
}

void printusage() {
  std::cout << "usage:  tupal  -i seq1 [-i seq2 -i ...] -o outfile [wnpcfItmsgGdaluPSWTMr] -p probelength -w wordlength " << std::endl;
  std::cout << "-i input sequence1 ... : sequences to find unique probes in." << std::endl;
  std::cout << "-o output              : output file." << std::endl;
  std::cout << "-t                     : print corresponding tuple " << std::endl;
  std::cout << "-r paramFile           : read params from file, rather than tediously typing them every time. " << std::endl;
  std::cout << "------------------------ k-mer alignment options  ------------------- " << std::endl;
  std::cout << "-p probelength         : length of the probes to search for (20)." << std::endl;
  std::cout << "-w wordlength          : number of nucleotides in hash table (14). " << std::endl;
  std::cout << "-n neighborThresh      : distance to search for a neighbor (0) " << std::endl;
  std::cout << "-I                     : search insertions and deletions (false)."<<std::endl;
  std::cout << "-m                     : mask repeats. " << std::endl;
  std::cout << "------------------------ parallelization options ------------------- " << std::endl;
  std::cout << "-f                     : number of fragments to break genome into " << std::endl;
  std::cout << "-c                     : number of child processes to "
	    << "start "  << std::endl;
  std::cout << "------------------------ recursive alignment options --------------- " << std::endl;
  std::cout << "-g minGapDist          : minimium distance of gap to recursively re-align (3) " 
	    << std::endl;
  std::cout << "-G maxGapDist          : maximum distance of gap to recursively re-align (infinity) " 
	    << std::endl;
  std::cout << "-a [s|g]               : align either gaps (g), or "
    "strips (s)" << std::endl;
  std::cout << "-s step                : amount to decrease hash"
    " length by at each recursion (0 implies no recursion)" <<
    std::endl;
  std::cout << "-l minHashLength       : minimum hash length to use in\n"
    " recursion (default probelength, no recursion) " << std::endl;
  std::cout << "-W lisWindow          : size of window to compute LIS for" << std::endl;
  std::cout << "-T lisThreshold       : threshold to consider adjacent integers for longest increasing subset"
	    << std::endl;
  std::cout << "-M mergeThreshold     : discard strips of sizes 1 to mergeThreshold when aligning " << std::endl;
}

