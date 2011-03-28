#include <iostream>
#include <numeric>
#include "bayes_quality.hpp"
#include "quality_read_stream.hpp"
#include "ifaststream.hpp"

LOGGER("b");

#define SKIP_READS 6
#define PROCESS_READS 5

using namespace bayes_quality;

void processQualityReads(const char *filename, BayesQualityGenome & bqg) {
	QualityReadStream qrs(filename);
	
	vector<QRead *> vec;
	for (size_t i=0; i<PROCESS_READS; ++i) {
		QRead *qr = new QRead(qrs.Next());
		vec.push_back( qr );
	}
	
	bqg.ProcessReads(vec);
	
	for (size_t i=0; i<PROCESS_READS; ++i) {
		delete vec[i];
	}
	vec.clear();
	
	return;
	
	while (!qrs.eof()) {
		QRead qr(qrs.Next());
		double res = bqg.ReadBQ(qr);
		INFO(bqg.LastMatchReadString());
		INFO(bqg.LastMatchPrettyString());
		INFO(bqg.LastMatchString());
		ostringstream m; for (size_t i=0; i< bqg.LastMatch().size(); ++i) m << bqg.LastMatch()[i]; INFO(m.str());
		INFO(res << ", best: " << bqg.LastMatchQ() << "/" << bqg.LastTotalQ() << " at " << bqg.LastMatchIndex() << " with " << bqg.LastMatchInserts() << " inserts and " << bqg.LastMatchDeletes() << " deletes");
	}
}

int main() {
	INFO("Hello, Bayes!");
	

	ifaststream ifs("/home/student/nikolenko/python/bayesQuality/biggenome.fasta");
	// ifaststream ifs("/home/student/nikolenko/python/bayesQuality/genome.fasta");
	string name, seq;
	ifs >> name >> seq;
	INFO(name.data());
	BayesQualityGenome bqg(seq.data());

	processQualityReads("/home/student/nikolenko/python/bayesQuality/s_6_1.fastq.gz", bqg);

	// processQualityReads("/home/student/nikolenko/python/bayesQuality/smalltest.fastq", bqg);
	
	return 0;
}
