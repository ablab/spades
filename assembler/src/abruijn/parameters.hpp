#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <sstream>
#include "ireadstream.hpp"

#define MPSIZE 100
#define K 25
#define HASH_SEED 1845724623

#define CUT 1000000

using namespace std;

static string itoa(int n) {
	std::stringstream ss;
	ss << n;
	return ss.str();
}

//static pair<string,string> inputFileNames = make_pair("./data/MG1655-K12_emul1.fasta.gz", "./data/MG1655-K12_emul2.fasta.gz");
static pair<string,string> inputFileNames = make_pair("./data/s_6_1.fastq.gz", "./data/s_6_2.fastq.gz");
static string outputFileName = "./output/a" + itoa(CUT);

static ireadstream<MPSIZE, 2> irs(inputFileNames.first.c_str(), inputFileNames.second.c_str());

#endif /* PARAMETERS_H_ */
