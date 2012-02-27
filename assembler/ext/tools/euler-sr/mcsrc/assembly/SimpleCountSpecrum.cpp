/***************************************************************************
 * Title:          SimpleCountSpecrum.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>

char rc(char c) {
	if (c == 'a' or c == 'A') 
		return 'T';
	if (c == 't' or c == 'T') 
		return 'A';
	if (c == 'c' or c == 'C') 
		return 'G';
	if (c == 'g' or c == 'G') 
		return 'C';
}
	
	
typedef char Tuple[17] ;
								
int main(int argc, char* argv[]) {
	std::string inFile = argv[1];
	
	std::ifstream in;

	in.open(inFile.c_str());
	ssize_t numReads = 0;
	std::string line;
	ssize_t readLength = -1;
	while(in) {
		std::getline(in, line);
		if (line[0] == '>')
			numReads++;
		else {
			if (readLength == -1) {
				readLength = line.size();
			}
		}
	}
	std::vector<std::string> tuples;
	
	tuples.resize(numReads*2*(readLength - 16 + 1));

	ssize_t i;
	ssize_t cur = 0;
	in.close(); in.clear();
	in.open(inFile.c_str());
	ssize_t p;
	ssize_t t;
	ssize_t bad;
	ssize_t numPos;
	while(in) {
		if (!std::getline(in, line))
			break;
		if (line[0] != '>') {
			for (p = 0; p < readLength -16 + 1; p++) {
				bad = 0;
				tuples[cur] = "                ";
				tuples[cur+1] = "                ";
				for (t = 0; t < 16; t++) {
					if (line[p+t] == '.') bad = 1;
					tuples[cur][t] = line[p+t];
					tuples[cur+1][16 - t - 1] = rc(line[p+t]);
				}
				if (!bad) {
					cur+=2;
				}
			}
		}
	}
		std::cout << "got " << cur << " tuples " << " of " << numReads*2*(readLength-16+1) 
						<< " sorting . " << std::endl;
		exit(0);
	tuples.resize(cur);
	std::sort(tuples.begin(), tuples.end());
	// now output the multiplicities of tuples.
	i = 1;
	ssize_t curMult = 1;
	while (i < cur) {
		if (tuples[i] == tuples[i-1])
			curMult++;
		else {
			std::cout << tuples[i-1] << " " << curMult << std::endl;
			curMult = 1;
		}
		i++;
	}				
	return 0;
}
	
