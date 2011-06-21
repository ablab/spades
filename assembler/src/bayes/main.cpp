#include <omp.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>
#include "bayes_quality.hpp"
#include "ireadstream.hpp"
#include "read.hpp"
#include "quality.hpp"
#include "strobe_read.hpp"

#include "common/logging.hpp"
DECL_PROJECT_LOGGER("b")

#define PROCESS_READS 25
#define READ_FILENAME "/home/student/nikolenko/python/bayesQuality/s_6_1.fastq.gz"
#define GENOME_FILENAME "data/bayes/biggenome.fasta"
#define BOWTIE_COMMAND "./src/libs/bowtie-0.12.7/bowtie -a -v 2 tmpgenome.fasta"
#define TMP_GENOME_FILE "tmpgenome.fasta"

using namespace bayes_quality;

/**
  * splits a string by delimiter tracking char numbers
  * omits zero size 
  */
vector<pair<string, pair<size_t, size_t> > > splitString(string s, char delim) {
	stringstream temp (stringstream::in | stringstream::out);
	vector<pair<string, pair<size_t, size_t> > > elems(0);
	if(s.size() == 0 || delim == 0) return elems;
	size_t index = 0;
	for (size_t i=0; i<s.length(); ++i) {
		if(s[i] == delim) {
			if (temp.str().length() > 0) {
				elems.push_back(make_pair(temp.str(), make_pair(index, i)));
			}
			temp.str("");
			index = i+1;
		} else temp << s[i];
	}
	if(temp.str().size() > 0) {
		elems.push_back(make_pair(temp.str(), make_pair(index, s.size())));
	}
	return elems;
}



/**
  * splits a sequence (read) into subsequences by Ns
  */
vector<Read> splitReadByN(const Read r) {
	vector<pair<string, pair<size_t, size_t> > > v = splitString(r.getSequenceString(), 'N');
	vector<Read> res;
	for (size_t i=0; i<v.size(); ++i) {
		res.push_back(r.getName() + '.' + (char)((int)('a') + i), v[i].first);						//quality
	}
	return res;
}


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


	cout << genomefilename.data() << "\n";

	ireadstream ifs(genomefilename.data());
	Read r1;
//	string genome;
//	while (!ifs.eof()) {
		ifs >> r1;
//		vector<Read> v = splitReadByN(r1);
//		for (size_t j = 0; j < v.size(); ++j) {
//			genome += v[j].getSequenceString();
//		}
//	}
	INFO(r1.getSequenceString().length());
	BayesQualityGenome bqg(r1.getSequenceString().data());
	
	#ifdef USE_BOWTIE
/*		ofstream os;
		os.open(TMP_GENOME_FILE);
		os << genome.data();
		os.close();*/
		bqg.setBowtie(bowtiecmd, bowtieindex);
	#endif

	cout << readfilename.data() << "\n";
	const string strings[1] = { readfilename };
	typedef SingleReader<Read, ireadstream>::type SingleStream;
	SingleStream irs( strings );
	vector<Read> r(1);
	irs >> r;
	cout << r[0].getSequenceString() << "\n";
	r[0].trimNs();
	INFO(r[0].getSequenceString());
	cout << r[0].getSequenceString() << "\n";
	irs.close();

	bqg.ProcessReads(readfilename.data(), toskip);
	
	INFO("total good positions: " << bqg.TotalGoodPositions() << "/" << bqg.TotalPositions());

	string command = "rm -rf ";
	command += TMP_GENOME_FILE;
	system(command.data());

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
