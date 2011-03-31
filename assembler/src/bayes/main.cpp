#include <iostream>
#include <numeric>
#include "bayes_quality.hpp"
//#include "quality_read_stream.hpp"
#include "ireadstream.hpp"
#include "read.hpp"
#include "quality.hpp"

LOGGER("b");

#define SKIP_READS 6
#define PROCESS_READS 5000

using namespace bayes_quality;

void processQualityReads(const char *filename, BayesQualityGenome & bqg) {
	
	vector<Read> *vec = ireadstream::readAll(filename, PROCESS_READS);
	vector<Read> tmpvec; for (size_t i=0; i < PROCESS_READS; ++i) tmpvec.push_back(vec->at(i));
	bqg.ProcessReads(tmpvec);
	delete vec;
	
	return;
	
	/*while (!qrs.eof()) {
		QRead qr(qrs.Next());
		double res = bqg.ReadBQ(qr);
		INFO(bqg.LastMatchReadString());
		INFO(bqg.LastMatchPrettyString());
		INFO(bqg.LastMatchString());
		ostringstream m; for (size_t i=0; i< bqg.LastMatch().size(); ++i) m << bqg.LastMatch()[i]; INFO(m.str());
		INFO(res << ", best: " << bqg.LastMatchQ() << "/" << bqg.LastTotalQ() << " at " << bqg.LastMatchIndex() << " with " << bqg.LastMatchInserts() << " inserts and " << bqg.LastMatchDeletes() << " deletes");
	}*/
}

int main() {
	INFO("Hello, Bayes!");
	

	//ireadstream ifs("/home/student/nikolenko/python/bayesQuality/biggenome.fasta.gz");
	ireadstream ifs("/home/student/nikolenko/python/bayesQuality/genome.fasta");
	Read r;
	ifs >> r;
	INFO("hello");
	BayesQualityGenome bqg(r.getSequence().data());

	processQualityReads("/home/student/nikolenko/python/bayesQuality/s_6_1.fastq.gz", bqg);
	
	INFO("total good positions: " << bqg.TotalGoodPositions() << "/" << bqg.TotalPositions());

	return 0;
}
