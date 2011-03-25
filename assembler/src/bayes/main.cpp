#include <iostream>
#include "bayes_quality.hpp"
#include "quality_read_stream.hpp"
#include "ifaststream.hpp"

LOGGER("b");

using namespace bayes_quality;

int main() {
	INFO("Hello, Bayes!");
	
	QualityReadStream qrs("/home/student/nikolenko/python/bayesQuality/test.fastq");
	std::vector<QRead *> qrv;
	while (!qrs.eof()) {
		QRead *qr = new QRead(qrs.Next());
		qrv.push_back(qr);
	}
	
	ifaststream ifs("/home/student/nikolenko/python/bayesQuality/biggenome.fasta");
	string name, seq;
	ifs >> name >> seq;
	cout << name.data() << "\n";

	BayesQualityGenome bqg(seq.data());
	for (size_t i=0; i<qrv.size(); ++i) {
		INFO(qrv[i]->first.str());
		double result = bqg.ReadBQInt(*qrv[i]);
		INFO(result);
	}
	
	for (size_t i=0; i<qrv.size(); ++i) delete qrv[i];
	
	return 0;
}
