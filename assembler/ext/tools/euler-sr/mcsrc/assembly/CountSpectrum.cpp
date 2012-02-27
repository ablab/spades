/***************************************************************************
 * Title:          CountSpectrum.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "compatibility.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "DeBruijnGraph.h"
#include "HashedSpectrum.h"
#include "ListSpectrum.h"
#include "MultTuple.h"
#include <vector>
#include <iostream>
#include <fcntl.h>
#include "IntegralTupleStatic.h"


void PrintUsage() {
	std::cout << "usage: countSpectrum reads kmer_out " << std::endl;
	std::cout << "  [-tupleSize tupleSize ] (20)  The size of tuple for which to count frequency." 
						<< std::endl;
	std::cout << "  [-noPrintCount|-printPos]     Print tuples as a list of indices into the read list." 
						<< std::endl;
  std::cout << "                                This can be a more space-efficient way of representing tuples." 
						<< std::endl;
	std::cout << "  [-trimFront T] Trim first T characters from each read before" 
						<< std::endl << "                  counting the spectrum" << std::endl;
	std::cout << "  [-trimEnd   T] Trim T characters from the end of each read"
						<< std::endl << "                  before ocunting the spectrum. " << std::endl;
	std::cout << "    [-sort]   Write tuples in sorted order." << std::endl;   
	std::cout << "  [-lockFile file]  Wait on locked file 'file' to read input." << std::endl;
	std::cout << "  This is convenient when running on a grid that has an NFS bottle neck." 
						<< std::endl;
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
		PrintUsage();
    exit(0);
  }
  std::string readsFile = argv[1];
  std::string tupleOutFile = argv[2];
  int argi = 3;
  SimpleSequenceList seqList;
  HashValueFunctor hashFunction;
  int tupleSize = 20;
  ssize_t printCount = 1;
  ssize_t printPosition = 0;
  ssize_t DoSort = 0;
	ssize_t trimFront = 0;
	ssize_t trimEnd   = 0;
	ssize_t sortTuples= 0;
	ssize_t doLock = 0;
	std::string lockFileName = "";
  hashFunction.hashLength = 11;

  while (argi < argc) {
    if (strcmp(argv[argi], "-tupleSize") == 0 or
				strcmp(argv[argi], "-vertexSize") == 0) {
      ++argi;
      tupleSize = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-noPrintCount") == 0) {
      printCount = 0;
    }
    else if (strcmp(argv[argi], "-printPos") == 0) {
      printPosition =1;
    }
    else if (strcmp(argv[argi], "-sort") == 0) {
      DoSort = 1;
    }
		else if (strcmp(argv[argi], "-trimFront") == 0){ 
			trimFront = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-trimEnd") == 0) {
			trimEnd   = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-sort") == 0) {
			sortTuples = 1;
		}
		else if (strcmp(argv[argi], "-lockFile") == 0) {
			doLock = 1;
			lockFileName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-hashLength") == 0) {
			// TODO: Not documted in PrintUsage.
			hashFunction.hashLength = atoi(argv[argi++]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(1);
		}
    ++argi;
  }

	std::string reportFileName = FormReportName(readsFile);
	ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);

  if (tupleSize < hashFunction.hashLength) 
    hashFunction.hashLength = tupleSize;


  CountedReadPos::hashLength = tupleSize;
  CountedReadPos::sequences  = &seqList;
	hashFunction.sequences     = &seqList;


  std::ifstream seqIn;

  DNASequence read, rc;
  //UNUSED// ssize_t curRead = 0;
  SimpleSequence simpleSeq;
  //UNUSED// ssize_t storeSeq;
  //UNUSED// ssize_t seqNumber = 0;
	HashedSpectrum spectrum(hashFunction);
	_INT_ lockFileDes;
	//UNUSED// struct flock fldes;
	if (doLock) {
		WaitLock(lockFileName, lockFileDes);
	}

  ReadSimpleSequences(readsFile, seqList, report);
	
	if (doLock) {
		ReleaseLock(lockFileDes);
	}
  AppendReverseComplements(seqList);
	spectrum.hashFunction.sequences = &seqList;
	spectrum.SetHashSize(tupleSize);
	spectrum.StoreSpectrum(seqList, trimFront, trimEnd );

  // now print the hash
  std::ofstream tupleOut;
  openck(tupleOutFile, tupleOut, std::ios::out, report);
  ReadPosHashTable::iterator it, end;
  
  ssize_t i;
  CountedReadPos readPos;
  ssize_t nReadPos = 0;
  // Step 1, count the number of read positions.
	ReadPosHashTable::HashPage *page;
	ReadPosHashTable::Data *data, *pageEnd;
  for (i = 0; i < hashFunction.MaxHashValue(); i++ ) {
		page = spectrum.hashTable.table[i].head;
		while (page != NULL) {
			pageEnd = &(*page).page[page->size];
			for (data = &(*page).page[0]; data != pageEnd; ++data) {
				++nReadPos;
			}
			page  = page->next;
		}
	}
	/*    end = spectrum.hashTable.table[i].End();
    it  = spectrum.hashTable.table[i].Begin();
    while (it < end) {
      ++nReadPos;
      ++it;
    }
		}*/

	if (sortTuples) {
		ListSpectrum<MultTuple> tupleList;
		tupleList.Resize(nReadPos);
		ssize_t index =0;
		//		MultTuple::tupleSize = tupleSize;
		MultTuple::SetTupleSize(tupleSize);
		for (i = 0; i < hashFunction.MaxHashValue(); i++ ) {
			page = spectrum.hashTable.table[i].head;
			while (page != NULL) {
				pageEnd = &(*page).page[page->size];
				for (data = &(*page).page[0]; data != pageEnd; ++data) {

					/*				
										end = spectrum.hashTable.table[i].End();
										it  = spectrum.hashTable.table[i].Begin();
										while (it < end) {
					*/
					tupleList[index].copy((char*) &seqList[data->read].seq[data->pos], tupleSize);
					++index;
				}
				page = page->next;
			}
		}
		std::sort(tupleList.tupleList.begin(), tupleList.tupleList.end());
		tupleOut.close();
		tupleOut.clear();
		if (doLock) {
			WaitLock(lockFileName, lockFileDes);
		}
		tupleList.Write(tupleOutFile);
	}
	else {
		tupleOut << nReadPos << std::endl;
		if (doLock) {
			WaitLock(lockFileName, lockFileDes);
		}
		for (i = 0; i < hashFunction.MaxHashValue(); i++ ) {
			page = spectrum.hashTable.table[i].head;
			while (page != NULL) {
				pageEnd = &(*page).page[page->size];
				for (data = &(*page).page[0]; data != pageEnd; ++data) {

					//					end = spectrum.hashTable.table[i].End();
					//			it  = spectrum.hashTable.table[i].Begin();
					//			while (it < end) {
					if (printPosition) {
						tupleOut << data->read << " " << data->pos << std::endl;
					}
					else {
						PrintTuple(seqList, *data, tupleSize, tupleOut); 
						if (printCount) {
							tupleOut << " " << data->count;
						}
						tupleOut << std::endl;
					}
				}
				page = page->next;
			}
		}
	}
	if (doLock){ 
		ReleaseLock(lockFileDes);
	}

	EndReport(report);
	report.close();
}
