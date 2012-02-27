/***************************************************************************
 * Title:          SortIntegralTupleList.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include "IntegralTupleStatic.h"

void PrintUsage() {
	std::cout << "usage: sortIntegralTupleList tupleList.spect" << std::endl;
	std::cout << "     -printCount   Both read in and print out the count of each tuple" << std::endl;
	std::cout << "     -minMult M    Discards tuples that occur less than M times." << std::endl
						<< std::endl;
	exit(0);
}

int main(int argc, char *argv[]) {
	
	std::string tupleListName;
	int argi = 1;
	if (argc < 2) {
		PrintUsage();
	}
	tupleListName = argv[argi++];
	ssize_t readCount = 0;
	ssize_t minMult   = 0;
	while(argi < argc) {
		if (strcmp(argv[argi], "-printCount") == 0) {
			readCount = 1;
		}
		else if (strcmp(argv[argi], "-minMult") == 0) {
			minMult = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
		}

		argi++;
	}

	std::string reportFileName = FormReportName(tupleListName);
	ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);

	
	//std::ifstream listIn;
	std::ofstream listOut;
	//openck(tupleListName, listIn, std::ios::in | std::ios::binary, report);
	ssize_t nTuples;
	//	listIn.read((char*) &nTuples, sizeof(int));

	if (readCount) {
		CountedIntegralTuple* tupleList;
		/*		CountedIntegralTuple *tupleList = new CountedIntegralTuple[nTuples];
		listIn.read((char*) tupleList, 
		sizeof(CountedIntegralTuple)*nTuples);

		listIn.close();
		*/
		ReadBinaryTupleList(tupleListName, &tupleList, nTuples, minMult, report);


		std::sort(tupleList, tupleList + nTuples);
		openck(tupleListName, listOut, std::ios::out | std::ios::binary, report);

		ssize_t tupleSize_SSZT = CountedIntegralTuple::tupleSize;
		listOut.write((const char*) &tupleSize_SSZT, sizeof(ssize_t));
		listOut.write((const char*) &nTuples, sizeof(ssize_t));
		listOut.write((const char*) tupleList, sizeof(CountedIntegralTuple)*nTuples);
	}
	else {
		IntegralTuple *tupleList;
		ReadBinaryTupleList(tupleListName, &tupleList, nTuples, 0, report);
		ssize_t tupleSize_SSZT = IntegralTuple::tupleSize;
		std::sort(tupleList, tupleList + nTuples);
		openck(tupleListName, listOut, std::ios::out | std::ios::binary, report);
		listOut.write((const char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
		listOut.write((const char*) &nTuples, sizeof(ssize_t));
		listOut.write((const char*) tupleList, sizeof(IntegralTuple)*nTuples);
	}

	EndReport(report);
	report.close();

	return 0;
}
