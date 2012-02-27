/***************************************************************************
 * Title:          CountPaired.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include "DNASequence.h"
#include "SeqReader.h"
#include "compatibility.h"
#include "../IntegralTupleStatic.h"


using namespace std;


class LongUnpairedTuple {
public:
	ssize_t pos;
	ssize_t masked;
	static char *genome;
	static ssize_t k;

	bool operator<(const LongUnpairedTuple &rhs) const {
		return (strncmp(&genome[pos], &genome[rhs.pos], k) < 0);
	}
	bool operator==(const LongUnpairedTuple &rhs) const {
		return (strncmp(&genome[pos], &genome[rhs.pos], k) == 0);
	}
	friend ostream& operator<<(std::ostream &out, const LongUnpairedTuple &pt) {
		string astr, bstr;
		astr.insert(0, &pt.genome[pt.pos], pt.k);
		out << pt.pos << " " << pt.k << " " << " " << astr << " " << pt.masked;
		return out;
	}
	LongUnpairedTuple& operator=(const LongUnpairedTuple &rhs) {
		if (this != &rhs) {
			pos = rhs.pos;
		}
		return *this;
	}
	LongUnpairedTuple() {
		pos = -1;
		masked = 0;
	}
};
	
class LongPairedTuple {
public:
	ssize_t pos;
	ssize_t masked;
	static char *genome;
	static ssize_t k;
	static ssize_t d;
	bool operator<(const LongPairedTuple &rhs) const {
		ssize_t resa, resb;
		//		cout << "comparing " << *this << " with: " << rhs << endl;
		if ((resa = strncmp(&genome[pos], &genome[rhs.pos], k)) == 0) {
			resb = strncmp(&genome[pos+k+d], &genome[rhs.pos+k+d], k);
			return resb < 0;
		}
		else {
			return resa < 0;
		}
	}
	bool operator==(const LongPairedTuple &rhs) const {
		//
		return (strncmp(&genome[pos], &genome[rhs.pos], k) == 0 and
						strncmp(&genome[pos + k + d], &genome[rhs.pos+k+d], k) == 0);
	}
	friend ostream& operator<<(std::ostream &out, const LongPairedTuple &pt) {
		string astr, bstr;
		astr.insert(0, &pt.genome[pt.pos], pt.k);
		bstr.insert(0, &pt.genome[pt.pos + pt.k + pt.d], pt.k);
		out << pt.pos << " " << pt.k << " " << pt.d << " " << astr << " " << bstr << " " << pt.masked;
		return out;
	}
	LongPairedTuple& operator=(const LongPairedTuple &rhs) {
		if (this != &rhs) {
			pos = rhs.pos;
		}
		return *this;
	}
	LongPairedTuple() {
		pos = -1;
		masked = 0;
	}
};	

ssize_t LongPairedTuple::k = 0;
ssize_t LongPairedTuple::d = 0;
char* LongPairedTuple::genome = 0;
ssize_t LongUnpairedTuple::k = 0;
char* LongUnpairedTuple::genome = 0;


class PairedTuple {
public:
	IntegralTuple a, b;
	ssize_t pos;
	bool operator<(const PairedTuple &t) const {

		if (a < t.a) return 1; // compare tuples
		if (a == t.a) return b < t.b; // compare tuples
		return 0;


#if 0 // OLD CODE
		if (a.tuple < t.a.tuple) return 1;
		if (a.tuple == t.a.tuple) return b.tuple < t.b.tuple;
		return 0;
#endif
	}
	PairedTuple& operator=(const PairedTuple &t) {
		if (this != &t) {
			a = t.a;
			b = t.b;
			pos = t.pos;
		}
		return *this;
	}

	bool operator==(const PairedTuple &t) const {
		return a == t.a and b == t.b;

#if 0 // OLD CODE
		return a.tuple == t.a.tuple and b.tuple == t.b.tuple;
#endif
	}

	PairedTuple(IntegralTuple _a, IntegralTuple _b, ssize_t _pos) {
		a = _a;
		b = _b;
		pos = _pos;
	}
	friend std::ostream &operator<<(std::ostream &out, const PairedTuple &t)  {
		out << t.a.tuple << " " << t.b.tuple << " " << t.pos;
		return out;
	}
};

void PrintUsage() {
	cout << "Usage: countPaired seqFile tupleSize span" << endl;
	cout << " -histogram h  Print the number of times repeat counts occur to h" << endl
			 << " -contigLengths file Print the lengths of the non-repetitive contigs to file" << endl
			 << " -printNumDuplicated 'file' Print the number of duplicated k,d-mers to 'file'" << endl
			 << " -printUniqueDuplicated 'file' Print the number of unique duplicated K,d-mers to file" << endl
			 << " -printNumContigs 'file' Print hte number of contigs t o'file" <<endl
			 << " -printN50  'file'  Print the N50 to 'file'" << endl
			 << " -printMaskedGenome 'file' Print the repeat masked genome to file" << endl;
	
}

int main(int argc, char* argv[]) {
	
	if (argc < 4) {
		PrintUsage();
		exit(1);
	}
	string seqFileName = argv[1];
	int tupleSize = atoi(argv[2]);
	ssize_t span = atoi(argv[3]);

	int argi = 4;
	string contigLengthsFileName = "";
	ssize_t printContigLengths = 0;
	string histogramFileName;
	ssize_t printHistogram = 0;
	string numDuplicatedFileName = "";
	string uniqueDuplicatedFileName = "";
	ssize_t printUniqueDuplicated = 0;
	ssize_t printNumDuplicated = 0;
	ssize_t printNumContigs = 0;
	string numContigsFileName;
	ssize_t printN50 = 0;
	string n50FileName;
	ssize_t printMaskedGenome =0 ;
	string seqOutFileName;
	while (argi < argc) {
		if (strcmp(argv[argi], "-contigLengths") == 0) {
			contigLengthsFileName = argv[++argi];
			printContigLengths = 1;
		}
		else if (strcmp(argv[argi], "-histogram") == 0) {
			histogramFileName = argv[++argi];
			printHistogram =1;
		}
		else if (strcmp(argv[argi], "-printNumDuplicated") == 0) {
			numDuplicatedFileName = argv[++argi];
			printNumDuplicated =1;
		}
		else if (strcmp(argv[argi], "-printUniqueDuplicated") == 0) {
			uniqueDuplicatedFileName = argv[++argi];
			printUniqueDuplicated = 1;
		}
		else if (strcmp(argv[argi], "-printNumContigs") == 0) {
			printNumContigs = 1;
			numContigsFileName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-printN50") == 0) {
			printN50 = 1;
			n50FileName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-printMaskedGenome") == 0) {
			printMaskedGenome = 1;
			seqOutFileName = argv[++argi];
		}
		else {
			PrintUsage();
			cout << "bad option: " << argv[argi] << endl;
			exit(1);
		}
		
		argi++;
	}

	DNASequence seq;
	SeqReader::GetSeq(seqFileName, seq);

	LongPairedTuple::k = tupleSize;
	LongPairedTuple::d = span;
	LongPairedTuple::genome = (char*) seq.seq;

	LongUnpairedTuple::k = tupleSize;
	LongUnpairedTuple::genome = (char*) seq.seq;
	IntegralTuple::SetTupleSize(tupleSize);
	//	IntegralTuple::tupleSize = tupleSize;
	
	ssize_t i;
	ssize_t end = seq.length - span - tupleSize*2+1;

	IntegralTuple a, b;
	std::vector<LongPairedTuple> pairedTupleList;
	std::vector<LongUnpairedTuple> unpairedTupleList;
	if (span > 0) {
		pairedTupleList.resize(end);
		for (i =0 ; i < end; i++ ){ 
			pairedTupleList[i].pos = i;
		}
		cout << "done assigning positions, sorting." << endl;
		std::sort(pairedTupleList.begin(), pairedTupleList.end());
		cout << "done sorting. " << endl;
	}
	else {
		unpairedTupleList.resize(end);
		for (i =0 ; i < end; i++ ){ 
			unpairedTupleList[i].pos = i;
		}
		cout << "done assigning positions, sorting." << endl;
		std::sort(unpairedTupleList.begin(), unpairedTupleList.end());
		cout << "done sorting. " << endl;
	}
	// Count matched pairs, and mask the genome.
	ssize_t lastI;
	ssize_t cur;
	lastI = pairedTupleList.size() - 1;
	ssize_t numMasked = 0;

	map<ssize_t,ssize_t> tupleMultHist;
	ssize_t numDuplicated = 0;
	ssize_t numUniqueDuplicated = 0;

	if (span > 0) {
		for (i = 0; i < end; i++) {
			cur = i;
			ssize_t mult = 0;
			while (i < end and 
						 pairedTupleList[i+1] == pairedTupleList[cur]) {
				pairedTupleList[cur].masked = 1;
				pairedTupleList[i+1].masked = 1;
				i++;
				numMasked++;
				mult++;
			}
			if (mult > 0) {
				if(tupleMultHist.find(mult) != tupleMultHist.end()) {
					tupleMultHist[mult]++;
				}
				else {
					tupleMultHist[mult] = 1;
				}
				++numUniqueDuplicated;
				numDuplicated += mult;
			}
		}
		for (i = 0; i < end; i++ ) {
			//		cout << pairedTupleList[i] << endl;
			if (pairedTupleList[i].masked) {
				seq.seq[pairedTupleList[i].pos] = 'N';
				seq.seq[pairedTupleList[i].pos + LongPairedTuple::k + LongPairedTuple::d] = 'N';
			}
		}
	}
	else {
		for (i = 0; i < end; i++) {
			cur = i;
			ssize_t mult = 0;
			while (i < end and 
						 unpairedTupleList[i+1] == unpairedTupleList[cur]) {
				unpairedTupleList[cur].masked = 1;
				unpairedTupleList[i+1].masked = 1;
				i++;
				numMasked++;
				mult++;
			}
			if (mult > 0) {
				if(tupleMultHist.find(mult) != tupleMultHist.end()) {
					tupleMultHist[mult]++;
				}
				else {
					tupleMultHist[mult] = 1;
				}
				++numUniqueDuplicated;
				numDuplicated += mult;
			}
		}
		for (i = 0; i < end; i++ ) {
			//		cout << unpairedTupleList[i] << endl;
			if (unpairedTupleList[i].masked) {
				seq.seq[unpairedTupleList[i].pos] = 'N';
			}
		}
	}

	std::vector<ssize_t> contigLengths;

	if (printContigLengths or printNumContigs or printN50) {
		ssize_t pos, contigStart, contigEnd;
		pos = 0;
		// start contig on an unmasked nuc
		
		if  (seq.length > 0) {
			while (pos < seq.length and seq.seq[pos] == 'N')
				pos++;
			contigStart = pos;
			// Move forward 
			while (pos < seq.length - 1) {
				if (seq.seq[pos] == 'N' and seq.seq[pos+1] != 'N')
					contigStart = pos+1;
				else if (seq[pos] != 'N' and seq.seq[pos+1] == 'N') {
					contigEnd = pos+1;
					contigLengths.push_back(contigEnd - contigStart);
				}
				pos++;
			}
			if (seq.seq[pos] != 'N') {
				contigLengths.push_back(pos - contigStart);
			}
		}
	}
	if (printContigLengths) {
		std::ofstream contigsOut;
		openck(contigLengthsFileName, contigsOut, std::ios::out);
		ssize_t c;
		for (c = 0; c < contigLengths.size(); c++ ){ 
			contigsOut << contigLengths[c] << endl;
		}
	}
	
	if (printMaskedGenome) {
		std::ofstream maskedGenomeOut;
		openck(seqOutFileName, maskedGenomeOut, std::ios::out);
		seq.namestr += " (masked)";
		seq.PrintlnSeq(maskedGenomeOut);
		maskedGenomeOut.close();
	}

	if (printN50) {
		std::sort(contigLengths.begin(), contigLengths.end());
		ssize_t c;
		_SSZT_ totalLength = 0;
		for (c = 0; c < contigLengths.size(); c++ ){ 

			totalLength += contigLengths[c];
			//			cout << contigLengths[c] << " " << totalLength << endl;
		}
		ssize_t accumLength = 0;
		for (c = 0; c < contigLengths.size() and (accumLength *2 < totalLength); c++) { 
			accumLength += contigLengths[c];
		}
		std::ofstream n50Out;
		openck(n50FileName, n50Out, std::ios::out);
		if (c != 0)
			n50Out << contigLengths[c-1] << endl;
		else 
			n50Out << "N/A" << endl;
	}
		
	if (printHistogram) {
		std::ofstream histogramOut;
		openck(histogramFileName, histogramOut, std::ios::out);
		std::map<ssize_t,ssize_t>::iterator histIt, histEnd;
		for( histIt = tupleMultHist.begin(), histEnd = tupleMultHist.end();
				 histIt != histEnd; ++histIt) {
			histogramOut << histIt->first << " " << histIt->second << endl;
		}
	}
	if (printNumDuplicated) {
		std::ofstream numDuplicatedOut;
		openck(numDuplicatedFileName, numDuplicatedOut, std::ios::out);
		numDuplicatedOut << numDuplicated << endl;
	}
	if (printUniqueDuplicated) {
		std::ofstream uniqueDuplicatedOut;
		openck(uniqueDuplicatedFileName, uniqueDuplicatedOut, std::ios::out);
		uniqueDuplicatedOut << numUniqueDuplicated << endl;
	}
	if (printNumContigs) {
		std::ofstream numContigsOut;
		openck(numContigsFileName, numContigsOut, std::ios::out);
		numContigsOut << contigLengths.size() << endl;
	}
}
