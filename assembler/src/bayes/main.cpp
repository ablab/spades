#include <iostream>
#include <numeric>
#include "bayes_quality.hpp"
#include "ireadstream.hpp"
#include "ifaststream.hpp"
#include "read.hpp"
#include "quality.hpp"

LOGGER("b");

#define SKIP_READS 6
#define PROCESS_READS 25

#define READ_FILENAME "/home/student/nikolenko/python/bayesQuality/s_6_1.fastq.gz"

using namespace bayes_quality;

void processQualityReads(const char *filename, BayesQualityGenome & bqg) {
		
	vector<Read> *vec = ireadstream::readAll(filename, PROCESS_READS);
	vector<Read> tmpvec; for (size_t i=0; i < PROCESS_READS; ++i) tmpvec.push_back(vec->at(i));
	bqg.ProcessReads(tmpvec);
	delete vec;
	
	return;
}

int main(int argc, char* argv[]) {
	INFO("Hello, Bayes!");
	
	string readfilename = "";
	if (argc > 1) readfilename = argv[1];
	else readfilename = READ_FILENAME;

	//ireadstream ifs("/home/student/nikolenko/python/bayesQuality/biggenome.fasta.gz");
	ifaststream ifs("./data/bayes/biggenome.fasta");
	string name, genome;
	ifs >> name >> genome;
	INFO("!" << name);
	BayesQualityGenome bqg(genome.data());

	processQualityReads(readfilename.data(), bqg);
	
	INFO("total good positions: " << bqg.TotalGoodPositions() << "/" << bqg.TotalPositions());

	return 0;
}
