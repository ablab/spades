/*
 * assembler.cpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include <cassert>
#include <iostream>
#include <list>
#include <cstdio>
#include <vector>
#include <ctime>
#include <array>
#include <string>
#include "ireadstream.hpp"

using namespace std;

#define MPSIZE 100
#define K 11

int main(int argc, char *argv[]) {
	std::cerr << "Hello, I am assembler!" << std::endl;
	time_t now = time(NULL);

	//vector<MatePair<MPSIZE> > *mate_pairs = FASTQParser<MPSIZE>::readAll("./data/short/s_6_1.fastq.gz", "./data/short/s_6_2.fastq.gz");

	cerr << "Time: " << (time(NULL) - now) << " sec." << endl;
	return 0;
}
