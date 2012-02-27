/***************************************************************************
 * Title:          IntegralCountSpectrum.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "VectorHashedSpectrum.h"
#include "utils.h"
#include "IntegralTupleStatic.h"
#include "hash/VectorHashTable.h"


//typedef	TupleHash<CountedIntegralTuple,11,32> CountedSpectrumHash;
//typedef TupleHash<IntegralTuple,11,32> SpectrumHash;

typedef	VectorTupleHash<CountedIntegralTuple,11, 100000> CountedSpectrumHash;
typedef VectorTupleHash<IntegralTuple,11, 65536> SpectrumHash;


using namespace std;

void WriteSpectrum(SpectrumHash *spectrum, std::ostream &spectOut, ssize_t printBinary, ssize_t directed) {
	ssize_t i;
	// straight duplication of code. need to learn more about templates to
	// get rid of this.
	//	SpectrumHash::HashTable::iterator end, it;
	ssize_t nTuples = 0;
	/*
		SpectrumHash::HashTable::Page *page;
	SpectrumHash::HashTable::Data *data, *pageEnd;
	*/
	for (i = 0; i < spectrum->hashTable.size; i++ ) {
		//		page = spectrum->hashTable.table[i].head;
		//UNUSED// ssize_t j;
		if (spectrum->hashTable.table[i] != NULL) {
			nTuples += spectrum->hashTable.table[i]->size();
		}
		/*				
		while (page != NULL) {
			pageEnd = &(*page).page[page->size];
			for (data = &(*page).page[0]; data != pageEnd; ++data) {
				++nTuples;
			}
			page = page->next;
		}
	*/
	}

	if (!directed)
		nTuples *=2;
	ssize_t tupleSize_SSZT = IntegralTuple::tupleSize;
	if (printBinary) {
		spectOut.write((const char*) &tupleSize_SSZT, sizeof(ssize_t));
		spectOut.write((const char*) &nTuples, sizeof(ssize_t));
	}
	else {
		spectOut << tupleSize_SSZT << std::endl;
		spectOut << nTuples << std::endl;
	}
	IntegralTuple tuple;
	std::string tupStr;

	for (i = 0; i < spectrum->hashTable.size; i++ ) {
	/*
		page = spectrum->hashTable.table.table[i].head;
		while (page != NULL) {
			pageEnd = &(*page).page[page->size];
			for (data = &(*page).page[0]; data != pageEnd; ++data) {
	*/
		if (spectrum->hashTable.table[i] != NULL) {
			ssize_t j;
			for (j = 0; j < spectrum->hashTable.table[i]->size(); j++) {
				if (printBinary) {
					spectOut.write((const char*) &(*spectrum->hashTable.table[i])[j],
												 sizeof(IntegralTuple));
					if (directed == 0) {
						IntegralTuple tupleRC;
						(*spectrum->hashTable.table[i])[j].MakeRC(tupleRC);
						spectOut.write((const char*) &tupleRC, sizeof(IntegralTuple));
					}
				}
				else {
					(*spectrum->hashTable.table[i])[j].ToString(tupStr);
					spectOut << tupStr << std::endl;
					if (directed == 0) {
						IntegralTuple tupleRC;
						(*spectrum->hashTable.table[i])[j].MakeRC(tupleRC);
						tupleRC.ToString(tupStr);
						spectOut << tupStr << std::endl;
					}
				}
				++nTuples;
			} // end for.
			//			page = page->next;
		}
	}
}


void WriteCountedSpectrum(CountedSpectrumHash *countedSpectrum, std::ostream &spectOut, ssize_t printBinary, ssize_t directed, ssize_t minMult) {
	ssize_t i;
	ssize_t nTuples = 0;
	//	CountedSpectrumHash::HashTable::iterator end, it;
	/*	CountedSpectrumHash::HashTable::Page *page;
	CountedSpectrumHash::HashTable::Data *data, *pageEnd;
	*/
	for (i = 0; i < countedSpectrum->hashTable.size; i++ ) {
		
		/*
		//end = countedSpectrum->hashTable.table.table[i].End();
		//		it  = countedSpectrum->hashTable.table.table[i].Begin();
		// 		page = countedSpectrum->hashTable.table.table[i].head;
		while (page != NULL) {
			pageEnd = &(*page).page[page->size];
			for (data = &(*page).page[0]; data != pageEnd; ++data) {
				if (data->count >= minMult)
					++nTuples;
			}
			page = page->next;
		}
		*/
		if (countedSpectrum->hashTable.table[i] != NULL) 
			nTuples += countedSpectrum->hashTable.table[i]->size();
		
	}
	if (directed == 0)
		nTuples*=2;
		
	ssize_t tupleSize_SSZT = IntegralTuple::tupleSize;
	if (printBinary) {
		std::cout << "wrote " << nTuples << " in binary to a file." << std::endl;
		spectOut.write((const char*) &tupleSize_SSZT, sizeof(ssize_t));
		spectOut.write((const char*) &nTuples, sizeof(ssize_t));
	}
	else {
		spectOut << tupleSize_SSZT << std::endl;
		spectOut << nTuples << std::endl;
	}
	IntegralTuple tuple;
	std::string tupStr;
	nTuples = 0;
	for (i = 0; i < countedSpectrum->hashTable.size; i++ ) {
		/*		page = countedSpectrum->hashTable.table.table[i].head;
		while (page != NULL) {
			pageEnd = &(*page).page[page->size];
			for (data = &(*page).page[0]; data != pageEnd; ++data) {
				++nTuples;
				if (data->count < minMult)
					continue;
				if (printBinary) {
					spectOut.write((char*) &*data, sizeof(CountedIntegralTuple));
					if (!directed) {
						CountedIntegralTuple tupleRC;
						data->MakeRC(tupleRC);
						tupleRC.count = data->count;
						spectOut.write((char*) &tupleRC, sizeof(CountedIntegralTuple));
					}
				}
				else {
					data->ToString(tupStr);
					spectOut << tupStr << " " << data->count << endl;
					if (!directed) {
						CountedIntegralTuple tupleRC;
						data->MakeRC(tupleRC);
						//					cout << "making rc " << tupleRC.tuple << endl;  
						tupleRC.ToString(tupStr);
						spectOut << tupStr << " " << data->count << endl;
					}
				}
			}
			page = page->next;
		*/
		if (countedSpectrum->hashTable.table[i] != NULL) {
			ssize_t j;
			for (j = 0; j < countedSpectrum->hashTable.table[i]->size(); j++) {
				if (printBinary) {
					spectOut.write((char*) &((*countedSpectrum->hashTable.table[i])[j]), 
												 sizeof(CountedIntegralTuple));
					if (!directed) {
						CountedIntegralTuple tupleRC;
						(*countedSpectrum->hashTable.table[i])[j].MakeRC(tupleRC);
						tupleRC.count = (*countedSpectrum->hashTable.table[i])[j].count;
						spectOut.write((char*) &tupleRC, sizeof(CountedIntegralTuple));
					}
				}
				else {
					(*countedSpectrum->hashTable.table[i])[j].ToString(tupStr);
					spectOut << tupStr << " " << (*countedSpectrum->hashTable.table[i])[j].count << endl;
					if (!directed) {
						CountedIntegralTuple tupleRC;
						(*countedSpectrum->hashTable.table[i])[j].MakeRC(tupleRC);
						//					cout << "making rc " << tupleRC.tuple << endl;  
						tupleRC.ToString(tupStr);
						spectOut << tupStr << " " << (*countedSpectrum->hashTable.table[i])[j].count << endl;
					}
				}
			}
		}
	}	
}

void PrintUsage() {
	cout << "usage: integralCountSpectrum reads tupleSize spectrum [-printCount] [-binary]" 
			 << endl;
	cout << "     -skipGapped   Do not count sequences with 'N' or '.' in them." 
			 << endl << endl
			 << "     -trimEnd n    Skip the trailing 'n' bases in every read."
			 << endl << endl
			 << "     -trimFront n  Skip the first 'n' bases in every read." 
			 << endl << endl
			 << "     -directed     Hash only the forward direction of every read " 
			 << endl
			 << "                    (the default is to hash forward and reverse complement)" 
			 << endl << endl
			 << "     -checkpoint file N  Flush the hash table to a file every 'N' added bases." 
			 << endl << endl
			 << "     -minMult M    Only print tuples that appear more than M times." 
			 << endl;
}

int main(int argc, char* argv[]) {
	
	string readsFileName, spectrumFileName;
	string cpBaseName;
	int tupleSize;
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	int argi = 1;
	readsFileName    = argv[argi++];
	tupleSize        = atoi(argv[argi++]);
	spectrumFileName = argv[argi++];

	ssize_t printCount  = 0;
	ssize_t printBinary = 0;
	ssize_t trimEnd     = 0;
	ssize_t trimFront   = 0;
	ssize_t skipGapped  = 0;
	ssize_t cpCount     = 0;
	ssize_t cpIndex     = 0;
	ssize_t directed    = 0;
	ssize_t minMult     = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-printCount") == 0) {
			printCount = 1;
		}
		else if (strcmp(argv[argi], "-binary") == 0) {
			printBinary = 1;
		}
		else if (strcmp(argv[argi], "-trimEnd") == 0) {
			trimEnd  = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-trimFront") == 0) {
			trimFront = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-skipGapped") == 0) {
			skipGapped = 1;
		}
		else if (strcmp(argv[argi], "-checkpoint") == 0) {
			cpBaseName = argv[++argi];
			cpCount    = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-directed") == 0) {
			directed = 1;
		}
		else if (strcmp(argv[argi], "-minMult") == 0) {
			minMult = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(1);
		}
		argi++;
	}
	std::ifstream readsIn;
	std::ofstream spectOut;
	std::ofstream report;
	string reportFileName = FormReportName(readsFileName);
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);

	openck(readsFileName, readsIn, std::ios::in, report);
	if (printBinary) 
		openck(spectrumFileName, spectOut, std::ios::out | std::ios::binary, report);
	else 
		openck(spectrumFileName, spectOut, std::ios::out, report);
	IntegralTuple::SetTupleSize(tupleSize);
	//	IntegralTuple::tupleSize = tupleSize;

	CountedSpectrumHash *countedSpectrum = (CountedSpectrumHash *) NULL;
	SpectrumHash *spectrum;
	// Creation of a spectrum uses a lot of memory, so only
	// make the required one.
	if (printCount) {
		countedSpectrum = new CountedSpectrumHash;
		spectrum = (SpectrumHash*) countedSpectrum;
	}
	else {
		spectrum = new SpectrumHash;
	}

	DNASequence read;
	DNASequence simpleRead;
	DNASequence readRC;
	ssize_t r = 0; 
	while(SeqReader::GetRead(readsIn, read)) { //, SeqReader::noConvert)) {
		simpleRead.seq = read.seq;
		simpleRead.length = read.length;
		if (skipGapped) {
			ssize_t seqIsGapped = 0;
			ssize_t curPos;
			for (curPos = 0; curPos < read.length; curPos++ ){ 
				if (unmasked_nuc_index[read.seq[curPos]] >= 4) {
					seqIsGapped = 1;
					break;
				}
			}
			if (seqIsGapped) {
				continue;
			}
		}
				
		if (printCount) {
			if (directed == 0)
				countedSpectrum->HashSequenceUndirected(simpleRead, trimFront, trimEnd);
			else
				countedSpectrum->HashSequence(simpleRead, trimFront, trimEnd);

			if (cpCount and countedSpectrum->hashTable.count >= cpCount) {
				stringstream outNameStream;
				outNameStream.str("");
				outNameStream << cpBaseName << "." << cpIndex << ".cp";
				ofstream cpOut;
				openck(outNameStream.str(), cpOut, std::ios::out, report);
				countedSpectrum->Flush();
				WriteCountedSpectrum(countedSpectrum, spectOut, printBinary, directed, 0);
				//				spectrum->hashTable.table.Free();
			}
		}
		else {
			if (directed == 0) 
				spectrum->HashSequenceUndirected(simpleRead, trimFront, trimEnd);
			else 
				spectrum->HashSequence(simpleRead, trimFront, trimEnd);

			if (cpCount and spectrum->hashTable.count >= cpCount) {
				stringstream outNameStream;
				outNameStream.str("");
				outNameStream << cpBaseName << "." << cpIndex << ".cp";
				ofstream cpOut;
				openck(outNameStream.str(), cpOut, std::ios::out, report);
				spectrum->FlushDirected();
				WriteSpectrum(spectrum, spectOut, printBinary, directed);
				//				spectrum->hashTable.table.Free();
			}
		}

		// Free up some memory
		read.Reset();
		PrintStatus(++r);
	}
	//	spectrum->hashTable.table.Summarize();
	cout << "read: " << r << " reads." << endl;

  //UNUSED// ssize_t i;
  // Step 1, count the number of read positions.
	
	if (printCount) {
		countedSpectrum->FlushDirected();
		countedSpectrum->hashTable.Summarize();
		WriteCountedSpectrum(countedSpectrum, spectOut, printBinary, directed, minMult);
	}

	else {
		spectrum->FlushDirected();
		WriteSpectrum(spectrum, spectOut, printBinary, directed);
	}

	std::cout << std::endl << "printed: " << spectrum->hashTable.count << " tuples." << std::endl;
	EndReport(report);
	report.close();
	return 0;
}

