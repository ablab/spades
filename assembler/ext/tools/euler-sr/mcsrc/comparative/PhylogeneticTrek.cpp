/***************************************************************************
 * Title:          PhylogeneticTrek.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <iterator>
#include <map>
#include <set>

#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "lav/LAVFile.h"
#include "lav/LAVReader.h"
#include "lav/LAVUtils.h"
#include "FragmentUtils.h"


ssize_t  ReadSequenceIds(std::string &seqIdFileName, SequenceIds &seqIds);
void StoreLociIds(Loci &loci, SequenceIds &lociIds);
ssize_t  ReadUnknownSequences(SequenceIds &unknownIds, 
			std::map<std::string, DNASequence*> &sequences,
			std::string &seqDir);
void ReadAlignment(std::string ref, std::string qry, 
		   std::string &alignmentDirectory,
		   LAVFile &alignment);


int main(int argc, char* argv[]) {

  std::string locusFileName;
  std::string seqIdFileName;
  std::string locusOutName;
  std::string sequenceDir, alignDir;
  std::string verboseOut;
  verboseOut = "";
  std::ofstream vbOut;
  sequenceDir = ".";
  if (argc < 4) {
    std::cout << "usage: " << argv[0] << " locusFileName seqIdFileName locusOutName" 
	      << std::endl;
    std::cout << "options:  [-seqdir seqdir] [-aligndir aligndir] [-verbose outFile] [-del]" 
	      << std::endl;
    std::cout << "           -seqdir and -aligndir specify the base directory for sequences and " << std::endl;
    std::cout << "           alignments. " << std::endl;
    std::cout << "           -verbose prints progress to outfile (must be specified) " << std::endl;
    std::cout << "           -del checks to see if the query has been deleted " << std::endl;
    exit(1);
  }
  int argi = 1;
  locusFileName = argv[argi++];
  seqIdFileName = argv[argi++];
  locusOutName  = argv[argi++];
  while (argi < argc) {
    if (strcmp(argv[argi], "-seqdir") == 0) {
      ++argi;
      sequenceDir = argv[argi];
    }
    if (strcmp(argv[argi], "-aligndir") == 0) {
      ++argi;
      alignDir = argv[argi];
    }
    if (strcmp(argv[argi], "-verbose") == 0) {
      ++argi;
      verboseOut = argv[argi];
      openck(verboseOut, vbOut, std::ios::out);
    }
    ++argi;
  }

  // Read in all necessary data
  Loci loci;
  if (ReadLoci(locusFileName, loci)) {
    std::cout << "error reading components " << std::endl;
    exit(1);
  }

  SequenceIds sequenceIds, knownIds, unknownIds;
  ReadSequenceIds(seqIdFileName, sequenceIds);
  StoreLociIds(loci, knownIds);
  std::set_difference(sequenceIds.begin(), sequenceIds.end(),
		      knownIds.begin(), knownIds.end(),
		      std::insert_iterator<SequenceIds>(unknownIds, 
							unknownIds.begin()));
  
  ReadLociSequences(loci, sequenceDir);


  // Cache the unknown sequences to reduce io
  std::map<std::string, DNASequence*> seqMap;
  ReadUnknownSequences(unknownIds, seqMap, sequenceDir);  
  SequenceIds::iterator seqIdIt;
  if (verboseOut != "") {
    vbOut << "read " << seqMap.size() + loci.size() 
	  << " sequences " << std::endl;
  }
  // Prepare the name of the coordinates of frames
  
  std::string qryPosName, queryName, sbjctName;
  MakeTempName(queryName, "qry.fasta");
  MakeTempName(sbjctName, "sbjct.fasta");
  ssize_t curSize;
  ssize_t curLocus = 0;
  ssize_t outer, inner;
  outer = 0; inner = 0;
  while (curLocus < loci.size()) {
    ssize_t l;
    curSize = loci.size();
    // search all loci that haven't been searched yet.
    if (verboseOut != "" ){
      vbOut << "on outer iter: " << outer << std::endl;
      ++outer;
    }
    inner = 0;
    for (l = curLocus; l < loci.size(); l++ ) {

      seqIdIt = unknownIds.begin();
      inner = 0;
      while ( seqIdIt != unknownIds.end() ){
	if (verboseOut != "") {
	  vbOut << "starting " << loci[l]->seqid << " " << *seqIdIt << std::endl;
	  vbOut << "on inner iter: " << inner << std::endl;
	}
	++inner;
	// Fetch the alignment (won't cache these... small file size
	// but much processing will be easy on fwgrid).c
	LAVFile alignment;
	//	std::cout << "reading blocks "; std::cout.flush();
	ReadAlignment(loci[l]->seqid, *seqIdIt, alignDir, alignment);
	//	std::cout << " done " << std::endl;
	std::vector<LAVBlock*> refOrderBlocks; 
	BlockReferenceOrder refOrder;
	//	std::cout << "sorting blocks "; std::cout.flush();
	StoreBlockArray(alignment, refOrderBlocks, refOrder);
	//	std::cout << "done " << std::endl;
	
	// Locate the region to search
	ssize_t sbjctStart, sbjctEnd;
	ssize_t foundExtension = 0;
        if (verboseOut != "") {
	  vbOut << "searching for " << loci[l]->seqid << " " << loci[l]->start << "\t" << loci[l]->end
	        << " in " << *seqIdIt << std::endl;
        }
	ssize_t regionFound;
	//	std::cout << "searching for surround region "; std::cout.flush();
	regionFound = FindSurroundingRegion(loci[l]->start, loci[l]->end,
					    refOrderBlocks,
					    sbjctStart, sbjctEnd,
					    seqMap[*seqIdIt]->length);
	//	std::cout << "done " << std::endl;
	if ( regionFound and 
	     sbjctStart < sbjctEnd and
	     (loci[l]->end - loci[l]->start)*10 > sbjctEnd - sbjctStart) {
	  if (verboseOut != "") {
	    vbOut << "found corresponding region: " << sbjctStart << " " << sbjctEnd << std::endl;
	    vbOut << "gap sizes: " << loci[curLocus] << loci[l]->end - loci[l]->start << " " 
		  << sbjctEnd - sbjctStart << std::endl;
          }
	  if (sbjctStart >= sbjctEnd) {
             ++seqIdIt;
	     continue;
          }
          DNASequence qrySeq, subSeq;
          qrySeq = loci[l]->locusSeq;
          qrySeq.HardMask();

	  // adjust the subject boundaries so that we search the nearby region
	  sbjctStart = std::max(0, sbjctStart - 5000);
	  sbjctEnd   = std::min(seqMap[*seqIdIt]->length - 1, sbjctEnd + 5000);

          subSeq.Copy(*seqMap[*seqIdIt], sbjctStart, sbjctEnd+1);
          subSeq.HardMask();
          if (verboseOut != "") {
//            vbOut << "checking seq: " << std::endl;
 //           qrySeq.PrintSeq(vbOut);  vbOut << std::endl;
          }
	  // No problems in locating the region
	  // Output the subject sequence to a file
	  //	  std::cout << "printing temporary sequences "; std::cout.flush();
	  PrintTempSeq(subSeq, 0, subSeq.length - 1, sbjctName);
	  PrintTempSeq(qrySeq, 0, qrySeq.length - 1, queryName);
	  //	  std::cout << " done " << std::endl;
	  Pos position;
	  ssize_t blastResult;
	  //	  std::cout << "running blast "; std::cout.flush();
	  blastResult = RunBlast(queryName, sbjctName, position);
	  //	  std::cout << std::endl;
	  if (blastResult) {
            if (verboseOut != "" ) {
	    vbOut << "sequence " << loci[l]->seqid << " found extension in " 
		  << *seqIdIt << " at "
		  << position.rBegin << " " << position.rEnd << " " 
		  << position.qBegin << " " << position.qEnd << " " 
		  << position.eValue << std::endl;
            }
            if (position.rBegin == -1 or position.qBegin == -1 or
             	position.rEnd == -1 or position.qEnd == -1) { 
	      continue;
            }
  	    Locus *locus = new Locus;
	    locus->seqid = *seqIdIt;
            if (position.qBegin < position.qEnd) {
	      /*
		if (position.qBegin - position.rBegin > 0) 
	        locus->start = sbjctStart + position.qBegin - position.rBegin;
		else
	        locus->start = sbjctStart;
		
		if (locus->start + loci[l]->locusSeq.length - 1 < sbjctEnd)
	        locus->end  = locus->start + loci[l]->locusSeq.length - 1;
		else
		locus->end  = sbjctEnd-1;
	      */
	      // more constrained attempt
	      locus->start = sbjctStart + position.qBegin;
	      locus->end   = sbjctStart   + position.qEnd;
            }
            else { 
	      /*
		if (position.qBegin - position.rEnd > 0)
                locus->start = sbjctStart + position.qBegin - position.rEnd;
		else
                locus->start = sbjctStart;
		
		if (locus->start + loci[l]->locusSeq.length -1 < sbjctEnd)
                locus->end  = locus->start + loci[l]->locusSeq.length - 1;
		else 
                locus->end  = sbjctEnd - 1;
	      */
	      locus->start = sbjctStart + position.qEnd;
	      locus->end   = sbjctStart + position.qBegin;
            }
	    locus->locusSeq.Copy(*seqMap[*seqIdIt], locus->start, locus->end+1);
            if (verboseOut != "" ) {
               vbOut << "new locus : " << locus->start << " " << locus->end << std::endl;
            }
	    loci.push_back(locus);
	    SequenceIds::iterator curIt = seqIdIt; 
	    ++curIt;
	    unknownIds.erase(seqIdIt);
	    seqIdIt = curIt;
	    foundExtension = 1;
	  }
	  // Clean up
	  system(std::string("rm " + sbjctName).c_str());
	  system(std::string("rm " + queryName).c_str());
	}
	if (!foundExtension) 
	  ++seqIdIt;
      }
    }
    // At this point the number of sites 
    // that have been searched are 
    curLocus = curSize;
  }

  PrintLoci(loci, locusOutName);
  if (verboseOut != "")
    vbOut.close();
  return 0;
}

void ReadAlignment(std::string ref, std::string qry, 
		   std::string &alignmentDirectory,
		   LAVFile &alignment) {
  std::stringstream alignNameStrm;
  alignNameStrm << alignmentDirectory << "/" 
		<< ref << "."  << qry 
		<< ".net.lav";
	
  // Read in the alignment, and store the blocks in
  // a searchable format
  LAVReader::ReadLAVFile(alignNameStrm.str(), alignment);
};
void StoreLociIds(Loci &loci, SequenceIds &lociIds) {
  ssize_t l;
  for (l = 0; l < loci.size(); l++)
    lociIds.insert(loci[l]->seqid);
}


ssize_t ReadUnknownSequences(SequenceIds &unknownIds, 
			 std::map<std::string, DNASequence*> &sequences,
			 std::string &seqDir) {
  SequenceIds::iterator seqIdIt;
  for (seqIdIt = unknownIds.begin(); 
       seqIdIt != unknownIds.end();
       ++seqIdIt) {
     
    std::stringstream nameStrm;
    nameStrm << seqDir << "/" << *seqIdIt << ".fasta";
    //    std::cout << "reading " << nameStrm.str() << std::endl;
    DNASequence *seq = new DNASequence;
    SeqReader::GetSeq(nameStrm.str(), *seq, SeqReader::noConvert);
    sequences[*seqIdIt] = seq;
  }
}




    
ssize_t ReadSequenceIds(std::string &seqIdFileName, SequenceIds &seqIds) {
  std::ifstream in;
  openck(seqIdFileName, in, std::ios::in);
  std::string id;
  while (in >> id) {
    seqIds.insert(id);
  }
}

