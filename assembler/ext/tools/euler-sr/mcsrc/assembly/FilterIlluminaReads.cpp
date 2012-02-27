/***************************************************************************
 * Title:          FilterIlluminaReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "utils.h"

using namespace std;

void PrintUsage() {
	cout << "usage: filterIlluminaReads in out [-max[ACTG] pct] [-minTrim trim] [-v]" << endl;
	cout << "    -maxC C (0.85)  Trim the ends of reads if they are more than C%"<<endl;
	cout << "    -maxA A (0.85)  Trim the ends of reads if they are more than A%"<<endl;
	cout << "    -maxG G (0.85)  Trim the ends of reads if they are more than G%"<<endl;
	cout << "    -maxT T (0.85)  Trim the ends of reads if they are more than T%"<<endl;
	cout << "    -maxN N       Trim a read if it has one nucleotide more than N% in the suffix" <<endl;
	cout << "    -minTrim T (25) only trim if the run is longer than T" << endl;
	cout << "    -v            Be verbose with trimming" << endl;
}

int main(int argc, char* argv[]) {

	string inFile, outFile;

	ifstream in;
	ofstream out;
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	
	inFile =  argv[1];
	outFile = argv[2];
	if (inFile == outFile) {
		cout << "ERROR: the input file and output file cannot be the same." << endl;
		exit(1);
	}
	int argi = 3;
	double maxC = 0.85;
	double maxA = 0.85;
	double maxG = 0.85;
	double maxT = 0.85;
	ssize_t verbose = 0;
	ssize_t minTrim = 25;
	while (argi < argc) {
		if (strcmp(argv[argi], "-maxC") == 0) {
			maxC = atof(argv[++argi]);
			if (maxC > 1.0) { maxC = maxC / 100.0; }
		}
		else if (strcmp(argv[argi], "-maxA") == 0) {
			maxA = atof(argv[++argi]);
			if (maxA > 1.0) { maxA = maxA / 100.0; }
		}
		else if (strcmp(argv[argi], "-maxG") == 0) {
			maxG = atof(argv[++argi]);
			if (maxG > 1.0) { maxG = maxG / 100.0; }
		}
		else if (strcmp(argv[argi], "-maxT") == 0) {
			maxT = atof(argv[++argi]);
			if (maxT > 1.0) { maxT = maxT / 100.0; }
		}
		else if (strcmp(argv[argi], "-maxN") == 0) {
			double trimPct = atof(argv[++argi]);
			if (trimPct > 1.0) { trimPct = trimPct / 100.0; }
			maxT = maxA = maxC = maxG = trimPct;
		}		
		else if (strcmp(argv[argi], "-minTrim") == 0){ 
			minTrim = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-v") == 0) {
			verbose = 1;
		}
		++argi;
	}
	openck(inFile, in, std::ios::in);
	openck(outFile, out, std::ios::out);

	DNASequence read;
	while(SeqReader::GetSeq(in, read, SeqReader::noConvert)) {
		ssize_t nA, nC, nT, nG;
		nT = nG = nA = nC = 0;
		ssize_t p = read.length - 1;
		ssize_t trim = 0;

		// Don't try and trim reads that are too short.
		if (read.length < minTrim) {
			read.PrintlnSeq(out);
			continue;
		}
		// Compute the pct A/C in the suffix.

		ssize_t end = read.length - minTrim;
		for (p = read.length - 1; p >= end; p--) {
			switch(toupper(read.seq[p])) {
			case 'A':
				nA++;
				break;
			case 'C':
				nC++;
				break;
			case 'T':
				nT++;
				break;
			case 'G':
				nG++;
				break;
			}
		}

		if ((double)nA / minTrim > maxA) {
			// This triggers a trim, continue moving back while the read is 
			trim = minTrim;
			while (p >= 0 and toupper(read.seq[p]) == 'A') {
				p--;
				trim++;
			}
			if (verbose) {
				cout << "trimming A " << trim << endl;
				read.PrintlnSeq(cout);
			}
			read.length -= trim;
		}
		else if ((double)nC / minTrim > maxC) {
			trim = minTrim;
			while (p >= 0 and toupper(read.seq[p]) == 'C') {
				p--;
				trim++;
			}
			if (verbose) {
				cout << "trimming C " << trim << endl;
				read.PrintlnSeq(cout);
			}
			read.length -= trim;
		}
		else if ((double)nT / minTrim > maxT) {
			trim = minTrim;
			while (p >= 0 and toupper(read.seq[p]) == 'T') {
				p--;
				trim++;
			}
			if (verbose) {
				cout << "trimming T " << trim << endl;
				read.PrintlnSeq(cout);
			}
			read.length -= trim;
		}
		else if ((double)nG / minTrim > maxG) {
			trim = minTrim;
			while (p >= 0 and toupper(read.seq[p]) == 'G') {
				p--;
				trim++;
			}
			if (verbose) {
				cout << "trimming G " << trim << endl;
				read.PrintlnSeq(cout);
			}
			read.length -= trim;
		}
		if (read.length > 0)
			read.PrintlnSeq(out);
	}
	return 0;
}	
