/***************************************************************************
 * Title:          ReadSimulator.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "SimpleStats.h"
#include "utils.h"

#include <string>
#include <iterator>
#include <stdlib.h>
#include <algorithm>
using namespace std;
void PrintUsage() {
  std::cout << "usage: readsimulator reference_seq output_file [options] " << std::endl;
  std::cout << "     -rl read_length (50)\t the length of a read. " << std::endl;
  std::cout << "     -rs read_stddev (0)\t standard deviation of read length sizes " << std::endl;
  std::cout << "     -ml mate_length (0)\t the distance spanning ends of reads in a mate-pair " << std::endl;
  std::cout << "                                  note: clone_length = 2*read_length + mate_length " << std::endl;
  std::cout << "     -ms mate_stddev (0)\t the deviation of mate lengths " << std::endl;
  std::cout << "     -coverage coverage (40)\t simulate coverageX reads. #reads = coverage*length(reference_seq)/read_length " << std::endl;
  std::cout << "     -stratify  spacing nPasses(0)\t rather than randomly simulating reads, " << std:: endl
						<< "                                  sample a read at 'spacing' between samples " << std::endl;
  std::cout << "                                  a value of 0 means random sampling.  Nonzero valus override coverage. " << std::endl;
	std::cout << "                        nPasses is the number of times to go through the genome and sample." << std::endl;
  std::cout << "     -mutate mutationRate (0)\t each nucleotide mutates with this probability " << std::endl;
  std::cout << "     -indel  indelRate (0)\t randomly create an indel with this frequency " << std::endl;
  std::cout << "     -startWell well (0)\t start counting the well at well (used for db reads) " << std::endl;
	std::cout << "     -plate P\t\tUse 'p' as plate to identify clone type." << std::endl;
  std::cout << "     -pairFile file    Print the mate pair information to 'file'" << std::endl;
  std::cout << "     -printGap gapFile Prints the sequence between mate ends. This " << std::endl 
						<< "                         is useful when benchmarking mate-end assembly " << std::endl;
  std::cout << "     -coordsFile file  Print the start and end coordinates of every read to 'file'" << std::endl;
  std::cout << "     -exactReads file  Print unmutated reads to 'file'" << std::endl;
  std::cout << "     -cloneLib cloneLen stddev pct  : Generate 'pct' percentage of the reads " << std::endl
						<< "                       using a clone length of cloneLen with deviation 'stddev' " << std::endl;
  std::cout << "     -ruleFile  file   Print clone naming conventions to 'file'." << std::endl;
	std::cout << "     -errorProfile file Read error probabilities from 'file'." << std::endl;
	std::cout << "     -catasProfile file Read catastrophe probabilities from 'file'" << std::endl;
	std::cout << "     -catasFraction fract Make a read catastrophic with probability 'fract'" << std::endl;
	std::cout << "     -readCoords file  Read read coordinates from 'file'." << std::endl;
}

ssize_t CopySequence(DNASequence &sourceSeq, ssize_t sourceStart, ssize_t sourceEnd,
									double mutationRate, double indelRate, DNASequence &newSeq, DNASequence &origSeq);

void AddPosDependentErrors(DNASequence &sourceSeq, DNASequence &newSeq,
													 std::vector<double> &errorProb, std::vector<ssize_t> &nerrors);

void MakeCatastrophe(DNASequence &sourceSeq, DNASequence &newSeq,
										 std::vector<double> &catasProb);

int main(int argc, char* argv[]) {
  ssize_t readLength;
  double readLengthStddev;
  ssize_t mateLength;
  ssize_t useCloneLib, useErrorProb, useCatasProb;
  double mateLengthStddev;
  ssize_t startWell;
	ssize_t plate;
  double mutationRate;
  double indelRate;
  ssize_t stratify;
	ssize_t nPasses;
  double coverage;
	double catasFraction = 0;
  std::string sourceSeqName;
  std::string readsFileName;
  std::string gapFile;
  std::string exactFile;
  std::string ruleFile;
	std::string errorProbFileName, catasProbFileName;
	errorProbFileName = "";
	catasProbFileName = "";
  readLength = 50;
  mateLength = 0;
  coverage   = 40;
  readLengthStddev = 0.0;
  mateLengthStddev = 0.0;
  mutationRate = -1.0;
  indelRate    = -1.0;
  startWell    = 0;
	plate        = 1;
  stratify     = 0;
	nPasses      = 0;
  useCloneLib  = 0;
	useErrorProb = 0;
	useCatasProb = 0;
  gapFile   = "";
  std::vector<ssize_t> cloneLengths;
  std::vector<ssize_t> cloneStddevs;
  std::vector<ssize_t> clonePcts;
	std::vector<double> errorProb, catasProb;
	//UNUSED// double catasFrac;
	
  if (argc < 3) {
    PrintUsage();
    exit(0);
  }
  SeedRandom();
  std::string pairFile = "";
  std::string coordsFile = "";
	std::string readCoordsFile = "";
	ssize_t doReadCoords  = 0;
  char cloneKeys[12][2] = {{'a', 'b'}, {'c', 'd'}, {'e', 'f'}, {'g', 'h'}, 
													 {'i', 'j'}, {'k', 'l'}, {'m', 'n'}, {'o', 'p'},
													 {'q', 'r'}, {'u', 'v'}, {'w', 'x'}, {'y', 'z'}};
  sourceSeqName = argv[1];
  readsFileName = argv[2];
  int argi = 3;
  ssize_t libML, libMSD, libPCT;
	
  while (argi < argc) {
    if (strcmp(argv[argi], "-rl") == 0) {
      argi++;
      readLength = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-rs") == 0) {
      argi++;
      readLengthStddev = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-ml") == 0) {
      argi++;
      mateLength = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-cloneLib") == 0){ 
      useCloneLib = 1;
      libML = atoi(argv[++argi]);
      libMSD = atoi(argv[++argi]);
      libPCT = atoi(argv[++argi]);
      if (cloneLengths.size() == 12) {
				std::cout << "ERROR: no more 13 different clone types are supported here." << std::endl;
				exit(1);
      }
      cloneLengths.push_back(libML);
      cloneStddevs.push_back(libMSD);
      clonePcts.push_back(libPCT);
    }
    else if (strcmp(argv[argi], "-ms") == 0) {
      argi++;
      mateLengthStddev = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-stratify") == 0) {
      stratify = atoi(argv[++argi]);
			nPasses  = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-coverage") == 0) {
      argi++;
      coverage = atof(argv[argi]);
    }
    else if (strcmp(argv[argi], "-mutate") == 0) {
      argi++;
      mutationRate = atof(argv[argi]);
    }
    else if (strcmp(argv[argi], "-indel") == 0) {
      argi++;
      indelRate = atof(argv[argi]);
    }
    else if (strcmp(argv[argi], "-pairFile") ==0) {
      pairFile = argv[++argi];
    }
    else if (strcmp(argv[argi], "-printGap") == 0) {
      gapFile = argv[++argi];
    }
    else if (strcmp(argv[argi], "-coordsFile") == 0) {
      coordsFile = argv[++argi];
    }
    else if (strcmp(argv[argi], "-exactReads") == 0) {
      exactFile = argv[++argi];
    }
    else if (strcmp(argv[argi], "-ruleFile") == 0) {
      ruleFile = argv[++argi];
    }
		else if (strcmp(argv[argi], "-errorProfile") == 0) {
			errorProbFileName = argv[++argi];
			useErrorProb = 1;
		}
		else if (strcmp(argv[argi], "-catasProfile") == 0) {
			catasProbFileName = argv[++argi];
			useCatasProb = 1;
		}
		else if (strcmp(argv[argi], "-catasFraction") == 0){ 
			catasFraction = atof(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-readCoords") == 0) {
			readCoordsFile = argv[++argi];
			doReadCoords = 1;
		}
		else if (strcmp(argv[argi], "-startWell") == 0) {
			startWell = atoi(argv[++argi]);
		} 
		else if (strcmp(argv[argi], "-plate") == 0) {
			plate = atoi(argv[++argi]);
		} 
		else {
      PrintUsage();
      std::cout << "bad option: " << argv[argi] << std::endl;
      exit(0);
    }
    argi++;
  }

  DNASequence sourceSeq;
  SeqReader::GetSeq(sourceSeqName, sourceSeq);

  ssize_t nReads;
  if (stratify == 0) {
    nReads= (ssize_t) (sourceSeq.length * coverage / readLength);
  }
  else {
    nReads = ((sourceSeq.length - readLength) / stratify) * nPasses;
  }
  std::cout << "using: " << nReads << " reads." << std::endl;
	std::ifstream errorProbFile;
	if (useErrorProb) {
		// read in the error probabilities in a slick way
		openck(errorProbFileName, errorProbFile, std::ios::in);
		std::copy(std::istream_iterator<double>(errorProbFile), 
							std::istream_iterator<double>(), 
							std::back_insert_iterator<std::vector<double> >(errorProb));
	}
	if (useCatasProb) {
		std::ifstream catasProbFile;
		openck(catasProbFileName, catasProbFile, std::ios::in);
		std::copy(std::istream_iterator<double>(catasProbFile), 
							std::istream_iterator<double>(), 
							std::back_insert_iterator<std::vector<double> >(catasProb));
		ssize_t c;
		// Represent as cdf for easy uniform-> parametric conversion
		for (c = 1; c < catasProb.size(); c++) {
			catasProb[c] = catasProb[c] + catasProb[c-1];
		}
		if (catasProb[catasProb.size()-1] == 0) {
			std::cout << "ERROR! 0 probability for something that should be close to 1!"<<std::endl;
			exit(0);
		}
		// normalize to what is 1.
		for (c = 0; c < catasProb.size(); c++ ){
			catasProb[c] /= catasProb[catasProb.size()-1];
		}
	}
	std::ifstream readCoordsIn;
	if (doReadCoords) {
		openck(readCoordsFile, readCoordsIn, std::ios::in);
	}
		
  std::ofstream readsOut;
  openck(readsFileName, readsOut, std::ios::out);

	std::ofstream shortReadsOut;
	openck(readsFileName + ".short", shortReadsOut, std::ios::out);

	if (cloneLengths.size() > 0 and pairFile == "") 
		pairFile = readsFileName + ".mates";

  std::ofstream pairOut, gapOut, exactOut;
  if (pairFile != "") 
    openck(pairFile, pairOut, std::ios::out);

  std::vector<char> cloneStart;
  if (gapFile != "" ) {
    openck(gapFile, gapOut, std::ios::out);
    cloneStart.resize(sourceSeq.length);
    std::fill(cloneStart.begin(), cloneStart.end(), 0);
  }

  if (exactFile != "") {
    openck(exactFile, exactOut, std::ios::out);
  }

  ssize_t n;
  //UNUSED// char readName[150];
  //UNUSED// ssize_t readStart;
  ssize_t adjustReadLength, adjustReadMateLength, adjustMateLength;
  ssize_t minSampleStart, maxSampleStart;
  ssize_t readDir = 0;
  DNASequence read, readRC, exactRead;
  read._ascii = 0; readRC._ascii = 0;
  //UNUSED// ssize_t readPos, seqPos;
  ssize_t sourceStart = -1, sourceEnd = -1; // TODO: check initialization values
  unsigned char* tmpRead;
  ssize_t maxReadLength = readLength * 3;
  //UNUSED// ssize_t mutatedReadLength = 0;
  tmpRead = new unsigned char[maxReadLength];
  //UNUSED// ssize_t fullReadLength;
  adjustReadMateLength = 0;
  ssize_t mateStart, mateEnd;
  ssize_t nSkipped = 0;
  std::ofstream readCoords, gapCoords;
  DNASequence exactRC;
	read._ascii = 1;
	readRC._ascii = 1;
	ssize_t nCatastrophe = 0;
	std::vector<ssize_t> nerrors;
	nerrors.resize(readLength);
	for (n = 0;n < readLength; n++ ) nerrors[n] = 0;
	
  // open all the files here, in case there are any erros we don't want to have
  // to wait until the very end.
  if (coordsFile != "") {
    openck(coordsFile, readCoords, std::ios::out);
    if (pairFile != "") {
      openck(coordsFile + ".gap", gapCoords, std::ios::out);
    }
  }
  std::ofstream ruleOut;

  if (stratify) {
    sourceStart = -stratify;
  }
  ssize_t curClone = 0;
  ssize_t numCurClone = 0;

  std::vector<ssize_t> cloneCounts;
  ssize_t i;
  for (i = 0; i < cloneLengths.size(); i++) {
    cloneCounts.push_back((ssize_t) floor(nReads/2 * (clonePcts[i]/100.0)));
  }
  if (useCloneLib) {
    mateLength = cloneLengths[curClone];
    mateLengthStddev = cloneStddevs[curClone];
  }
  else {
    mateLength = 0;
    mateLengthStddev = 0;
  }
  // Calculate the cumulative percetages so that
  // when creating a stratified genome, the
  // clones may be spread across the specified clone library
  std::vector<ssize_t> cumPcts;
  ssize_t cumPct = 0;
  for(i = 0; i < clonePcts.size(); i++) {
    cumPct += clonePcts[i];
    cumPcts.push_back(cumPct);
  }

  n=0;
	ssize_t curReadIndex = 0;
	ssize_t curPass = 0;
	ssize_t totalErrors = 0;
	ssize_t totalNucleotides = 0;
  while(doReadCoords or n < nReads) {
    numCurClone++;
    n++; // output at least one sequecne each loop.  two for mate-pairs
    // If sampling at specific intervals, and at the end of the genome
    // stop sampling.
    //    std::cout << "stratify: " << stratify << std::endl;
		
		if (doReadCoords) {
			if (!(readCoordsIn >> readDir >> sourceStart))
				// done reading sequences
				break;
		}
		

    if (stratify and sourceStart > sourceSeq.length - (readLength )) {
	++curPass;
        if (curPass >= nPasses)  {
          break;
        }
        else {
           sourceStart = -stratify;
        }
     } 

    // If we are simluating a clone library, check to see 
    // if the quota for the current clone has been met.  If so
    // move to the next clone size, or start generating unpaired reads.
    if (stratify ) {
      sourceStart += stratify;
      ssize_t cloneBin = Random(100);
      //		  std::cout << "clonebin: " << cloneBin;
      for (i = 0; i < cumPcts.size(); i++ ) {
				//		    std::cout << " " << cumPcts[i];
				if (cloneBin < cumPcts[i]) {
					break;
				}
      }
      //		  std::cout << " found bin: " << i << std::endl;
      if (i <cumPcts.size()) {
				curClone = i;
				useCloneLib      = 1;
				mateLengthStddev = cloneStddevs[curClone];
				mateLength       = cloneLengths[curClone];
      }
      else {
				useCloneLib = 0;
				mateLengthStddev = 0;
				mateLength = 0;
      }
    }

    if (!stratify and useCloneLib) {
      if (numCurClone > cloneCounts[curClone]) {
				if (curClone == cloneCounts.size() - 1) {
					// on the last clone in the library, use unpaired reads here out
					useCloneLib = 0;
					mateLengthStddev = 0;
					mateLength = 0;
				}
				else {
					// Start generating clones from the next size insert
					curClone++;
					mateLengthStddev = cloneStddevs[curClone];
					mateLength       = cloneLengths[curClone];
					numCurClone = 0;
				}
      }
    }

    // Nudge the read and mate lengths a bit.
    adjustReadLength = (ssize_t) (readLength + NormalRandom(0, readLengthStddev));
    // the read's mate
    if (useCloneLib) {
      adjustReadMateLength = (ssize_t) (readLength + NormalRandom(0, readLengthStddev));
      adjustMateLength     = (ssize_t) (cloneLengths[curClone] + NormalRandom(0, mateLengthStddev));
    }
    else {
      adjustReadMateLength = 0;
      adjustMateLength = 0;
    }
		
    if (stratify and sourceStart + adjustReadLength + adjustReadMateLength + adjustMateLength >= sourceSeq.length) {
      // The mate would sample past the end of the sequence, don't use it.
      useCloneLib = 0;
      adjustReadMateLength = 0;
      adjustMateLength = 0;
      mateLength = 0;
    }

		// Simulate the read in a random direction.
		if (!doReadCoords) {
			if (stratify > 0 or Random(2) > 0 ) {
				// make read in forward direction
				readDir = 0;
			}
			else {
				// read in reverse
				readDir = 1;
			}
		}

    // Pick a spot in the genome to get a read

		if (readDir == 0) {
			minSampleStart = 0;
			maxSampleStart = (ssize_t) (sourceSeq.length - (adjustReadLength + adjustReadMateLength + adjustMateLength));
		}
		else {
			minSampleStart = adjustReadLength + adjustReadMateLength + adjustMateLength;
			maxSampleStart = sourceSeq.length - adjustReadLength;
		}

    if (stratify == 0 and doReadCoords == 0) {
      sourceStart = minSampleStart + Random(maxSampleStart - minSampleStart);
      sourceEnd   = sourceStart + adjustReadLength;
    }
    else if (stratify) {
      sourceEnd   = sourceStart + adjustReadLength;
    }
		if (doReadCoords) {
			if (sourceStart + readLength > sourceSeq.length) 
				sourceEnd= sourceSeq.length;
			else
				sourceEnd = sourceStart + readLength;
		}
    read.Reset(maxReadLength);
		//		std::cout << "source start: " << sourceStart << " se: " << sourceEnd;
		read.Reset(0);
		exactRead.Reset(0);

		if (sourceEnd >= sourceSeq.length) {
			sourceStart -= (sourceSeq.length - sourceEnd);
			sourceEnd -= (sourceSeq.length - sourceEnd);
		}
		if (sourceStart <= 0) {
			sourceEnd += -sourceStart;
			sourceStart = 0;
		}
    totalErrors += 
			CopySequence(sourceSeq, sourceStart, sourceEnd,  mutationRate, indelRate, read, exactRead);
		
		totalNucleotides += read.length;
    if (readDir == 1) {
			readRC.Reset(0);
			exactRC.Reset(0);
      MakeRC(read, readRC);
      read.seq = 0;
      read = readRC;
			
      MakeRC(exactRead, exactRC);
      exactRead = exactRC;
    }

		if (useErrorProb) {
			AddPosDependentErrors(exactRead, read, errorProb, nerrors);
		}
		if (useCatasProb) {
			
			if (Uniform() < catasFraction) {
				MakeCatastrophe(exactRead, read, catasProb);
				nCatastrophe++;
			}
		}

    std::stringstream titlestrm;
    if (useCloneLib) 
      titlestrm << "sim_" << plate << "_" << curReadIndex + startWell << "." << cloneKeys[curClone][0];
    else 
      titlestrm << "sim_" << plate << "_" << curReadIndex + startWell << ".s";

    read.namestr = titlestrm.str();
		read.AddPosKey(sourceStart);
		read.AddStrandKey(readDir);
		
    read.PrintSeq(readsOut,200);
    readsOut << std::endl;

		read.length = 35;
		read.PrintSeq(shortReadsOut,200);
		shortReadsOut << std::endl;

    if (coordsFile !="") {
      readCoords << sourceStart << " " << sourceEnd << std::endl;
    }
    if (exactFile != "" ){
      exactRead._ascii = 0;
      exactRead.namestr = read.namestr;
      exactRead.PrintSeq(exactOut,200);
      exactOut << std::endl;
    }
    // There is a mate for this read
    if (mateLength > 0 ) {
			if (stratify == 0)
				n++; // count
			if (readDir == 0) {
				mateStart = sourceStart + adjustReadLength + adjustMateLength;
				mateEnd   = mateStart   + adjustReadMateLength;
			}
			else {
				mateStart = sourceStart - (adjustMateLength + adjustReadMateLength);

				mateEnd   = mateStart + adjustReadMateLength;
			}
			
			// Shift up if necessary
			if (mateStart < 0) {
				mateEnd  += -mateStart;
				mateStart = 0;
			}

			// shift down if necessary.
			if (mateEnd >= sourceSeq.length) {
				mateStart -= (sourceSeq.length - mateEnd);
				mateEnd   -= (sourceSeq.length - mateEnd);
			}

			read.Reset(0);
			exactRead.Reset(0);
      CopySequence(sourceSeq, mateStart, mateEnd, mutationRate, indelRate, read, exactRead);
			
      if (readDir == 0) {
				readRC.Reset(0);
				MakeRC(read, readRC);
				read.seq = 0; // unpoint this sequence
				read = readRC;
				exactRC.Reset(0);
				MakeRC(exactRead, exactRC);
				exactRead = exactRC;
				exactRead._ascii = 0;
      }
      titlestrm.str("");
      titlestrm << "sim_" << plate << "_" << curReadIndex + startWell << "." << cloneKeys[curClone][1];
      read.namestr = titlestrm.str();
      read.PrintSeq(readsOut,200);
      readsOut << std::endl;
			
			read.length = 35;
      read.PrintSeq(shortReadsOut,200);
      shortReadsOut << std::endl;
			

      if (coordsFile != "") {
				readCoords << mateStart <<" " << mateEnd << std::endl;
      }
      if (pairFile != "") {
				pairOut << "sim_" << plate << "_" << (n-2) + startWell << " " << (n-2) + startWell 
								<< " " << (n-1) + startWell << " " << curClone << std::endl;
      }

      if (exactFile != "" ){ 
				exactRead.namestr = read.namestr;
				exactRead.PrintSeq(exactOut,200);
				exactOut << std::endl;
      }
      if (gapFile != "" ) {
				// In order to not print too many full gapped sequences,
				// only print those that are not overlapped by other mates
				ssize_t searchSpan = readLength - 10;
				//UNUSED+// ssize_t p;
				ssize_t i ;
				ssize_t searchStart = (sourceStart + sourceEnd)/2 - searchSpan;
				if (searchStart < 0)
					searchStart = 0;
				ssize_t searchEnd   = (sourceStart + sourceEnd)/2 + searchSpan;
				if (searchEnd >= sourceSeq.length)
					searchEnd = sourceSeq.length-1;

				ssize_t overlap = 0;
				for (i = searchStart; i < searchEnd; i++ ) {
					if (cloneStart[i]) {
						overlap = 1;
						break;
					}
				}
				if (overlap == 0) {
					cloneStart[sourceStart] = 1;
				
					//UNUSED// ssize_t gapLength = mateEnd - sourceStart;
					read.Reset(0);
					exactRead.Reset(0);
					CopySequence(sourceSeq, sourceStart, mateEnd, 0, 0, read, exactRead);
					
					*read.titlestream << "gap_" << plate << "_" << n + startWell;
					read.PrintSeq(gapOut,200);
					gapOut << std::endl;
					if (coordsFile != "") {
						gapCoords << sourceStart << " " << mateEnd << std::endl;
					}
				}
				else {
					nSkipped++;
				}
      }
    }
		curReadIndex++;
  }
  if (gapFile != "") {
    std::cout << "skipped writing " << nSkipped << std::endl;
  }
	cout << totalErrors << " errors out of " << totalNucleotides << " bases." << endl;
}

ssize_t CopySequence(DNASequence &sourceSeq, ssize_t sourceStart, ssize_t sourceEnd,
									double mutationRate, double indelRate,
									DNASequence &newSeq, DNASequence &origSeq) {
  //UNUSED// ssize_t length = sourceEnd - sourceStart;
  ssize_t readPos;
  ssize_t seqPos;
  readPos = 0;
  seqPos = sourceStart;
  //UNUSED// double randVal;
	//  newSeq.Reset(length);
	//  std::vector<unsigned char> newSeqVect;
  origSeq.Copy(sourceSeq, sourceStart, sourceEnd);
	newSeq.Copy(sourceSeq, sourceStart, sourceEnd);
	ssize_t i, cur;
	ssize_t numErrors = 0;
	for (i = cur = 0; i < newSeq.length and cur < newSeq.length; i++) {
    if (Uniform() < mutationRate) {
			newSeq.seq[cur] = MutateNuc(origSeq.seq[i]);
      cur++;
			++numErrors;
    }
		else if (Uniform() < indelRate) {
			++numErrors;
			if (Uniform() < 0.5) {
				// insert
				newSeq.seq[cur] = RandomNuc();
				cur++;
			}
			else {
				// delete
				i++;
			}
		}
		else {
			newSeq.seq[cur] = origSeq[i];
			cur++;
		}
  }
	return numErrors;
}

void AddPosDependentErrors(DNASequence &sourceSeq, DNASequence &newSeq,
													 std::vector<double> &errorProb, std::vector<ssize_t> &nerrors) {
	ssize_t i;
	sourceSeq.Copy(newSeq);
	assert(errorProb.size() >= sourceSeq.length);
	for (i = 0;i < sourceSeq.length; i++) {
		if (Uniform() < errorProb[i]) {
			newSeq.seq[i] = tolower(MutateNuc(sourceSeq.seq[i]));
			nerrors[i]++;
		}
	}
}

void MakeCatastrophe(DNASequence &sourceSeq, DNASequence &newSeq,
										 std::vector<double> &catasProb) {
	//UNUSED// ssize_t catasPos;
	double catasRand = Uniform();
	ssize_t i;
	sourceSeq.Copy(newSeq);
	for (i = 0; i < catasProb.size()-1; i++ ){
		if (catasRand >= catasProb[i] && catasRand <= catasProb[i+1]) 
			break;
	}
	for (; i < sourceSeq.length; i++) {
		newSeq.seq[i] = tolower(MutateNuc(newSeq.seq[i]));
	}
}
