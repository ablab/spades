/***************************************************************************
 * Title:          CombineIntegralTupleLists.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntegralTupleStatic.h"
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

ssize_t FilesRemain(std::vector<ifstream*> &files) {
	ssize_t i;
	
	for (i = 0; i < files.size(); i++ ){
		if (*files[i])
			return 1;
	}
	return 0;
}

ssize_t AdvanceFile(ifstream &file, CountedIntegralTuple &ref, CountedIntegralTuple &cur) {
	if (cur == ref) {
		if (!file)
			return 0;
		if (file) {
			file.read((char*) &cur, sizeof(CountedIntegralTuple));
			return 1;
		}
	}
	return 1;
}

ssize_t AdvanceFile(ifstream &file, ssize_t &tupleIndex, ssize_t numTuples, CountedIntegralTuple &cur) {
	if (tupleIndex >= numTuples)
		return 0;
	file.read((char*) &cur, sizeof(CountedIntegralTuple));
	tupleIndex++;
	if (file) return 1;
	else return 0;
}

void GetMin(vector<CountedIntegralTuple> &tupleList,
						CountedIntegralTuple &min) {
	ssize_t i;
	min.count = -1;
	for (i =0 ;i < tupleList.size(); i++ ){
		if (min.count == -1 or 
				min > tupleList[i] // compare tuples
#if 0
				min.tuple > tupleList[i].tuple
#endif
				) {
			min = tupleList[i];
		}
	}
}


void CountCurTupleMultiplicity(vector<CountedIntegralTuple> &tupleList,
															 CountedIntegralTuple &curTuple) {
	ssize_t i;
	for (i = 0; i < tupleList.size(); i++ ){ 
		if (tupleList[i] == curTuple) {
			curTuple.count += tupleList[i].count;
		}
	}
}



int main(int argc, char* argv[]) {

	string combinedFileName;
	ssize_t minMult;
	vector<string> inputFileNames;

	int tupleSize = 0;

	combinedFileName = argv[1];
	minMult = atoi(argv[2]);
	int argi = 3;
	while(argi < argc) {
		if (strlen(argv[argi]) > 0 and argv[argi][0] == '-') {
			// parse an option
			cout << "No options are supported." << endl;
			exit(0);
		}
		else {
			inputFileNames.push_back(argv[argi]);
		}
		++argi;
	}

	std::vector<ifstream*> inFiles;
	ssize_t i;
	ifstream *fptr;
	ssize_t numInFiles = inputFileNames.size();
	cout << "searching: " << numInFiles << endl;
	inFiles.resize(numInFiles);
	for (i = 0; i < numInFiles; i++ ){
		fptr = new ifstream;
		inFiles[i] = fptr;
	}

	ofstream combinedOut;
	openck(combinedFileName, combinedOut, std::ios::out | std::ios::binary);
	
	std::vector<CountedIntegralTuple> tuples;
	std::vector<CountedIntegralTuple> curTuples;
	std::vector<ssize_t> inCurTuples;
	std::vector<ssize_t> inNumTuples;
	ssize_t nFiles = inputFileNames.size();
	curTuples.resize(nFiles);
	inCurTuples.resize(nFiles); std::fill(inCurTuples.begin(), inCurTuples.end(), 0);
	inNumTuples.resize(nFiles); std::fill(inNumTuples.begin(), inNumTuples.end(), 0);
	ssize_t stage;
	ssize_t nTuples = 0;
	for (stage = 0; stage < 2; stage++) {
		if (stage == 1) {
			ssize_t tupleSize_SSZT = tupleSize;
			combinedOut.write((const char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
			combinedOut.write((const char*) &nTuples, sizeof(ssize_t));
			nTuples = 0;
		}

		for (i = 0; i < numInFiles; i++) {
			// open the input files in stage 0
			// reset then reopen the input files in stage 1.
			if (stage == 0) {
				openck(inputFileNames[i], *inFiles[i], std::ios::in | std::ios::binary);

				ssize_t tupleSize_SSZT;
				(*inFiles[i]).read((char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
				if (tupleSize != tupleSize_SSZT) {
					if (tupleSize == 0) {
						tupleSize = tupleSize_SSZT;
						CountedIntegralTuple::SetTupleSize(tupleSize);
					} else {
						std::cerr << "Warning: tupleSize=" << tupleSize_SSZT 
											<< " in file[i]=\"" << inputFileNames[i]
											<< "\" but was " << tupleSize << " in previous files." << std::endl;
					}
				}

				(*inFiles[i]).read((char*) &inNumTuples[i], sizeof(ssize_t));
			}
			else {
				(*inFiles[i]).clear();
				(*inFiles[i]).seekg(0);
				std::fill(inCurTuples.begin(), inCurTuples.end(), 0);
				for (i = 0; i < numInFiles; i++) {
					inFiles[i]->read((char*) &inNumTuples[i], sizeof(ssize_t));
				}
			}
		}
		//
		// Prep the files with the first tuple in each file.
		//

		ssize_t fileAdvanced = 0;
		ssize_t adv;
		for (i = 0; i < numInFiles; i++) {
			adv = AdvanceFile(*inFiles[i], inCurTuples[i], inNumTuples[i], curTuples[i]);
			fileAdvanced |= adv;
		}
		if (fileAdvanced == 0) {
			// No files were advanced, done.
			combinedOut.close();
			return 0;
		}
		CountedIntegralTuple minTuple;
		while (fileAdvanced) {
			//
			// Find the next tuple to output

			// 
			// First time around, there is no min.
			//
			minTuple.count = -1;
			GetMin(curTuples, minTuple);

			// See how often it occurs in the input files.
			CountCurTupleMultiplicity(curTuples, minTuple);
			//			cout << "got min: " << minTuple.tuple << " " << minTuple.count << endl;
			if (minTuple.GetMult() >= minMult) {
				nTuples++;
				if (nTuples % 100000 ==0 ){
					cout << "counted: " << nTuples << endl;
				}
				if (stage == 0) { 
				}
				else {
					combinedOut.write((const char*) &minTuple, sizeof(CountedIntegralTuple));
				}
			}
			fileAdvanced = 0;
			
			for (i =0 ; i < numInFiles; i++) { 
				if (curTuples[i].tuple == minTuple.tuple) {
					adv = AdvanceFile(*inFiles[i], inCurTuples[i], inNumTuples[i], curTuples[i]);
					fileAdvanced |= adv;
				}
			}
			if (fileAdvanced == 0) {
				if (stage == 1) {
					// 
					// All done.
					//
					combinedOut.close();
					return 0;
				}
				else {
					cout << "counted: " << nTuples << " frequent tuples." << endl;
				}
			}
		}
	}
}
