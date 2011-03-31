#include <iostream>
#include <numeric>
#include "bayes_quality.hpp"
//#include "quality_read_stream.hpp"
#include "ireadstream.hpp"
#include "ifaststream.hpp"
#include "read.hpp"
#include "quality.hpp"

LOGGER("b");

#define SKIP_READS 6
#define PROCESS_READS 25

using namespace bayes_quality;

void processQualityReads(const char *filename, BayesQualityGenome & bqg) {
		
	vector<Read> *vec = ireadstream::readAll(filename, PROCESS_READS);
	vector<Read> tmpvec; for (size_t i=0; i < PROCESS_READS; ++i) tmpvec.push_back(vec->at(i));
	bqg.ProcessReads(tmpvec);
	delete vec;
	
	return;
}

int main() {
	INFO("Hello, Bayes!");
	

	//ireadstream ifs("/home/student/nikolenko/python/bayesQuality/biggenome.fasta.gz");
	ifaststream ifs("./data/bayes/biggenome.fasta");
	string name, genome;
	ifs >> name >> genome;
	INFO("!" << name);
	BayesQualityGenome bqg(genome.data());

	processQualityReads("/home/student/nikolenko/python/bayesQuality/s_6_1.fastq.gz", bqg);
	
	INFO("total good positions: " << bqg.TotalGoodPositions() << "/" << bqg.TotalPositions());

	return 0;
}
