/***************************************************************************
 * Title:          ContigMapper.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/25/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"
#include "MapContigs.h"
#include "SortedTupleList.h"
#include "utils.h"
#include "mctypes.h"
#include <fstream>
#include "ContigMap.h"
#include "IntegralTupleStatic.h"


void PrintUsage() {
	std::cout << "contigMapper refSeq contigFile" << std::endl;
	std::cout << "   [-tupleSize t]  (20) Find exact matches of size t" << std::endl;
	std::cout << "   [-readTable table_name]      Read a table of reference positions instead of compute" << std::endl;
	std::cout << "   [-writeTable table_name]     Write a table of reference positions" << std::endl;
	std::cout << "   [-writeContigMap mapname] Write the coordinates of contigs to a file " << std::endl;
	std::cout << "   [-writeUnique]  (false) When writing coordinates, only print uniquely mapped contigs." << std::endl;
	std::cout << "   [-minBlockCount c] (10) Discard mapped blocks with less than c words" << std::endl;
}

int main(int argc, char* argv[]) {

	std::string refSeqFile, contigFile;

	std::ifstream refIn, contigIn;
	
	if (argc <= 2) {
		PrintUsage();
		exit(1);
	}

	refSeqFile = argv[1];
	contigFile = argv[2];

	int argi;
	argi = 3;
	int tupleSize = 20;
	std::string table = "";
	std::string tableOutName = "";
	std::string contigMapName = "";
	ssize_t readTable = 0;
	ssize_t writeTable = 0;
	ssize_t writeUnique = 0;
	ssize_t writeContigMap = 0;
	ssize_t minBlockCount = 10;

	while (argi < argc) {
		if (strcmp(argv[argi], "-tupleSize") == 0) {
			tupleSize = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-readTable") == 0) {
			table = argv[++argi];
			readTable = 1;
		}
		else if (strcmp(argv[argi], "-writeTable") == 0){ 
			tableOutName = argv[++argi];
			writeTable = 1;
		}
		else if (strcmp(argv[argi], "-writeUnique") == 0) {
			writeUnique = 1;
		}
		else if (strcmp(argv[argi], "-writeContigMap") == 0) {
			contigMapName = argv[++argi];
			writeContigMap = 1;
		}
		else if (strcmp(argv[argi], "-minBlockCount") == 0) {
			minBlockCount = atoi(argv[++argi]);
		}
		else {
			std::cout << "bad option: " << argv[argi] << std::endl;
			PrintUsage();
			exit(1);
		}
		argi++;
	}

	DNASequence refSeq, refSeqRC;
	openck(refSeqFile, refIn, std::ios::in);
	if (!SeqReader::GetSeq(refIn, refSeq, SeqReader::noConvert)) {
		std::cout << "could not read" << refSeqFile << std::endl;
		exit(1);
	}
	MakeRC(refSeq, refSeqRC);
	
	SimpleSequenceList ref;
	ref.resize(2);
	ref[0].seq = refSeq.seq;
	ref[0].length = refSeq.length;
	ref[1].seq = refSeqRC.seq;
	ref[1].length = refSeqRC.length;
	
	ssize_t p;
	std::vector<CountedReadPos> refPositions;

	if (readTable == 0) {
		std::cout << "making sorted list. " << std::endl;
		MakeSortedTupleList(ref, tupleSize, refPositions);
		
		std::cout << refPositions.size() << " positions." << std::endl;
		RemoveDuplicatedTuples(ref, tupleSize, refPositions);
		std::cout << " after removing duplicates: " << refPositions.size() << std::endl;

		if (writeTable == 1) {
			std::ofstream tableOut;
			openck(tableOutName, tableOut, std::ios::out);
			tableOut << refPositions.size() << std::endl;
			for (p = 0; p < refPositions.size(); p++) {
				tableOut << refPositions[p] << std::endl;
			}
		}
	}
	else {
		std::cout << "reading table" << std::endl;
		std::ifstream tableIn;
		openck(table, tableIn);
		ssize_t numPos;
		tableIn >> numPos;
		refPositions.resize(numPos);
		for (p = 0 ;p < numPos; p++) {
			tableIn >> refPositions[p];
		}
		tableIn.close();
	}
	DNASequenceList contigs;
	ReadDNASequences(contigFile, contigs);

	//UNUSED+// ssize_t b;
	ssize_t c ;
	typedef std::list<Block> BlockList;
	BlockList blocks;
	BlockList::iterator blockIt;
	ContigMap map;
	if (writeContigMap) {
		std::ofstream contigMap;
		openck(contigMapName, contigMap, std::ios::out);
		for (c = 0; c < contigs.size(); c++) { 
			blocks.clear();
			MapContig(contigs[c], ref, refPositions, tupleSize, minBlockCount, blocks);
			//		std::cout << "mapped contig of size " << contigs[c].length << std::endl;
			if (!writeUnique or blocks.size() == 1) {
				//				contigMap << c << " " << blocks.size() << " ";
				for (blockIt = blocks.begin(); blockIt != blocks.end(); ++blockIt) {
					map.refPos = (*blockIt).refPos;
					map.refEnd = (*blockIt).refPos + (*blockIt).length;
					map.qryPos = (*blockIt).qryPos;
					map.qryEnd = (*blockIt).qryPos + (*blockIt).length;
					map.length = (*blockIt).length;
					//					contigMap << map;
					if ((*blockIt).strand == 0) {
						contigMap << (*blockIt).refPos << " " << (*blockIt).length + (*blockIt).refPos << " " 
											<< (*blockIt).length<< std::endl;
					}
					else {
						contigMap << refSeq.length - ((*blockIt).refPos + (*blockIt).length) << " "
											<< refSeq.length - ((*blockIt).refPos) << " "
											<< (*blockIt).length << std::endl;
					}
					/*					contigMap << " " << (*blockIt).refPos << " " 
										<< (*blockIt).refPos + (*blockIt).length << " "
										<< (*blockIt).qryPos << " " << (*blockIt).length;*/
				}
				/*				if (blocks.size() == 0) {
					contigMap << " " << contigs[c].length;
				}
				contigMap << std::endl;
				*/
			}
			else if (writeUnique and blocks.size() != 1) {
				//				contigMap << c  << " 0 " << contigs[c].length << std::endl;
			}
		}
		contigMap.close();
	}

	/*

	for (c = 0; c < contigs.size(); c++) { 
		std::cout << "Alignment of " << contigs[c].namestr << " " << contigs[c].length << std::endl;

		blocks.clear();
		MapContig(contigs[c], ref, refPositions, tupleSize, minBlockCount, blocks);
		//		std::cout << "mapped contig of size " << contigs[c].length << std::endl;
		b = 0;
		for (blockIt = blocks.begin(); blockIt != blocks.end(); ++blockIt) {
			std::cout << b << " " << (*blockIt).refPos << " ... " 
								<< (*blockIt).refPos + (*blockIt).length << " "
								<< (*blockIt).qryPos 
								<< " " << (*blockIt).length << std::endl;
		}
		b++;
	}
	*/
}
