#include <omp.h>
#include <iostream>
#include <numeric>
#include "bayes_quality.hpp"
#include "ireadstream.hpp"
#include "read.hpp"
#include "quality.hpp"

#include "common/logging.hpp"
DECL_PROJECT_LOGGER("b")

#define PROCESS_READS 25
#define READ_FILENAME "/home/student/nikolenko/python/bayesQuality/s_6_1.fastq.gz"
#define GENOME_FILENAME "data/bayes/biggenome.fasta"
#define BOWTIE_COMMAND "./src/libs/bowtie-0.12.7/bowtie -a -v 2 data/bayes/biggenome"

using namespace bayes_quality;

int main(int argc, char* argv[]) {

	#pragma omp parallel
	{
	INFO("Hello from thread " << omp_get_thread_num() << " out of " << omp_get_num_threads());
	INFO("Hello, Bayes!");
	}
	
	string readfilename = "";
	string genomefilename = "";
	string bowtieindex = "";
	string bowtiecmd = "";
	size_t toskip = 0;
	for (int i=1; i < argc; ++i) {
		if (!strcmp(argv[i], "--bowtie") || !strcmp(argv[i], "-b")) {
			bowtiecmd = argv[i+1];
		}
		if (!strcmp(argv[i], "--bindex") || !strcmp(argv[i], "-bi")) {
			bowtieindex = argv[i+1];
		}
		if (!strcmp(argv[i], "--reads") || !strcmp(argv[i], "-r")) {
			readfilename = argv[i+1];
		}
		if (!strcmp(argv[i], "--genome") || !strcmp(argv[i], "-g")) {
			genomefilename = argv[i+1];
		}
		if (!strcmp(argv[i], "--skip") || !strcmp(argv[i], "-s")) {
			toskip = atoi(argv[i+1]);
		}
	}
	if (readfilename == "") readfilename = READ_FILENAME;
	if (genomefilename == "") genomefilename = GENOME_FILENAME;
	if (bowtiecmd == "") bowtiecmd = BOWTIE_COMMAND;
	if (toskip == 0) toskip = SKIP_READS;


	ireadstream ifs(genomefilename.data());
	Read r1;
	ifs >> r1;
	INFO("!" << r1.getName());
	BayesQualityGenome bqg(r1.getSequenceString().data());
	
	#ifdef USE_BOWTIE
		bqg.setBowtie(bowtiecmd, bowtieindex);
	#endif

	ireadstream irs(readfilename.data());
	Read r;
	irs >> r;
	r.trimNs();
	INFO(r.getSequence());
	irs.close();

	bqg.ProcessReads(readfilename.data(), toskip);
	
	INFO("total good positions: " << bqg.TotalGoodPositions() << "/" << bqg.TotalPositions());

	return 0;
}


/*
void processQualityReads(const char *filename, BayesQualityGenome & bqg) {
	ireadstream ifs(filename);
	if (PROCESS_READS > 0) {
		INFO("processing " << PROCESS_READS << " reads");
		vector<Read> tmpvec;
		Read r;
		for (size_t i=0; i < PROCESS_READS; ++i) {
			ifs >> r; r.trimNs();
			if (r.size() < MIN_READ_SIZE) {
				--i; continue;
			}
			INFO(r.getSequenceString());
			tmpvec.push_back(r);
		}
		bqg.ProcessReads(tmpvec);
	} else {
		vector<Read> *vec = ireadstream::readAll(filename, PROCESS_READS);
		delete vec;
	}
	
	return;
}
*/
