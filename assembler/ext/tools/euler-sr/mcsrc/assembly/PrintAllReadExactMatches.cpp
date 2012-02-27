/***************************************************************************
 * Title:          PrintAllReadExactMatches.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"
#include "SeqUtils.h"
#include "IntegralTupleStatic.h"

void PrintSeq(char *seq) {
	char *s = new char[51];
	strncpy(s, seq, 50);
	s[50] = '\0';
	std::cout << s;
}

void ParseName(std::string name, ssize_t &pos, ssize_t &strand) {
	ssize_t start = 1;
	while(start < name.size() && name[start] != '_')
		start++;
	
	// pass the _
	start++;
	if (name.substr(start, 3) == "FOR")
		strand = 0;
	else
		strand = 1;
	
	start += 4;
	assert(start < name.size());
	assert(name[start] != '_');
	std::string posStr = name.substr(start);
	pos = atoi(posStr.c_str());
}

class SeqPos {
public:
	char *str;
	ssize_t  pos;
	ssize_t count;
	static ssize_t len;
	static ssize_t printSeq;
	bool operator< (const SeqPos& p) const {
		if (printSeq) {
			PrintSeq(this->str);
			std::cout << " ";
			PrintSeq(p.str);
			std::cout <<std::endl;
		}
		return (strncmp(this->str, p.str,len) < 0);
	}
	bool operator<= (const SeqPos &p) const {
		ssize_t cmpVal = strncmp(str, p.str, len);
		return (cmpVal <= 0);
	}
	bool operator> (const SeqPos& p) const {
		return (strncmp(str, p.str, len) > 0);
	}
	bool operator>= (const SeqPos &p) const {
		ssize_t cmpVal =strncmp(str, p.str, len);
		return (cmpVal >= 0);
	}

	bool operator== (const SeqPos& p) const {
		return (strncmp(str, p.str, len) == 0);
	}
};

class CompareSeqPos {
public:
	ssize_t operator()(const SeqPos &p1, const SeqPos &p2) {
		return p1 < p2;
	}
};
ssize_t LookupRead(std::vector<SeqPos>&genome, DNASequence &read);

ssize_t SeqPos::len = 0;
ssize_t SeqPos::printSeq = 0;
int main(int argc, char* argv[]) {

	std::string genomeName, readsFileName;
	ssize_t readLength;
	if (argc < 3) {
		std::cout << "usage: printAllReadExactMatches genome reads readLen" << std::endl;
		exit(1);
	}
	genomeName = argv[1];
	readsFileName  = argv[2];
	readLength = atoi(argv[3]);
	SeqPos::len = readLength;
	std::vector<SeqPos> positions;
	DNASequence genome, genomeRC;
	SeqReader::GetSeq(genomeName, genome, SeqReader::noConvert);
	MakeRC(genome, genomeRC);
	positions.resize(genome.length * 2);
	ssize_t p;
	for (p = 0; p < genome.length; p++) {
		positions[p].str = (char*) &genome.seq[p];
		positions[p].pos = p;
		positions[genome.length + p].str = (char*) &genomeRC.seq[p];
		positions[genome.length + p].pos = genome.length + p;
	}
	std::sort(positions.begin(), positions.end(), CompareSeqPos());
	p = 1;
	ssize_t count = 1;
	while (p < genome.length*2) {
		while (positions[p] == positions[p-1] && p < genome.length*2) {
			p++;
			count++;
		}
		positions[p-1].count = count;
		count = 1;
		p++;
	}
	
	std::ifstream readsIn;
	openck(readsFileName, readsIn, std::ios::in);
	DNASequence read;
	//UNUSED// ssize_t index;
	SeqPos readSeq;
	//UNUSED// ssize_t start;
	std::vector<SeqPos>::iterator indexIt, curIt;
	std::stringstream namestrm;
	ssize_t strand, pos;
	ssize_t readIndex = -1;
	std::string forStr = "FOR";
	std::string revStr = "REV";
	std::string strandStr;
	while (SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {
		++readIndex;
		//		index = LookupRead(positions, read);
		readSeq.str = (char*) read.seq;
		readSeq.pos = 0;
		
		indexIt = std::lower_bound(positions.begin(), positions.end(), readSeq, CompareSeqPos());

		//		std::cout << read.namestr << " " << indexIt->pos << std::endl;
		curIt = indexIt;
		++indexIt;
		if (curIt == positions.end() || !(*curIt == *indexIt)) {
			// just print this one read.
			//			read.PrintSeq(std::cout);
			ParseName(read.namestr, pos, strand);
			namestrm.str("");
			namestrm << read.namestr << " pos=" << pos << " strand=" << strand;
			read.namestr = namestrm.str();
			read.PrintlnSeq(std::cout);
			continue;
		}
		ssize_t copy = 1;
		if (*indexIt == *curIt) {
			ParseName(read.namestr, pos, strand);
			if (curIt->pos >= genome.length) {
				pos = curIt->pos - genome.length - read.length;
				strand = 1;
				strandStr = "REV";
			}
			else {
				pos = curIt->pos;
				strand = 0;
				strandStr = "FOR";
			}
			namestrm.str("");
			namestrm << read.namestr << " pos=" << pos << " strand=" << strand << " copy=0";
			
			read.namestr = namestrm.str();
			read.PrintlnSeq(std::cout);

			while (indexIt != positions.end() and *indexIt == *curIt) {

				if (indexIt->pos >= genome.length) {
					pos = indexIt->pos - genome.length - read.length;
					strand = 1;
					strandStr = "REV";
				}
				else {
					pos = indexIt->pos;
					strand = 0;
					strandStr = "FOR";
				}
				namestrm.str("");
				
				namestrm << readIndex << "_" << strandStr
								 << "_" << pos << " pos=" << pos << " strand=" << strand << " copy=" << copy;
				read.namestr = namestrm.str();
				read.PrintlnSeq(std::cout);
				
				indexIt++;
				++copy;
			}

		}
	}
		
	return 0;
}

ssize_t LookupRead(std::vector<SeqPos>&genome, DNASequence &read) {
	SeqPos readSeq;
	readSeq.str = (char*) &read.seq[0];
	readSeq.pos = 0;
	ssize_t beg = 0;
	ssize_t end = genome.size();
	ssize_t mid = (end + beg) / 2;
	
	while (mid > beg && mid <= end) {
		if (readSeq == genome[mid]){ 
			return mid;
		}
		if (readSeq < genome[mid]) {
			end = mid;
		}
		if (readSeq > genome[mid]) {
			beg = mid + 1;
		}
		mid = (end + beg)/2;
	}
	if (readSeq == genome[mid]) {
		return mid;
	}
	else {
		return -1;
	}
}

