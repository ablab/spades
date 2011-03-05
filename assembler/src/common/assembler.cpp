/*
 * assembler.cpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include "ireadstream.hpp"
#include <cassert>
#include <iostream>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
using namespace std;

#define MPSIZE 100
#define K 11
#define filename1 "./data/MG1655-K12_emul1.fasta.gz"
#define filename2 "./data/MG1655-K12_emul2.fasta.gz"
//#define filename1 "./data/s_6_1.fastq.gz"
//#define filename2 "./data/s_6_2.fastq.gz"

int main(int argc, char *argv[]) {
	cerr << "Hello, I am assembler!" << endl;
	time_t now = time(NULL);

	// read all 'read's


	ireadstream<100,2,int> irs(filename1, filename2);
	vector<mate_read<100,int>::type> *v = irs.readAll();
	irs.close();
	cerr << "Total reads (mate, without Ns): " << v->size() << endl;
	cerr << "Current time: " << (time(NULL) - now) << " sec." << endl;

	// construct graph



	// simplify graph


	// output graph

	cerr << "Total Time: " << (time(NULL) - now) << " sec." << endl;
	return 0;
}
