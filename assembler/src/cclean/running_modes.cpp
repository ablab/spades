#include "running_modes.hpp"
#include "QcException.hpp"
#include "aho_corasick.hpp"
#include "output.hpp"
#include "io/read_processor.hpp"
#include "job_wrappers.hpp"
#include "logger/log_writers.hpp"

void exactMatch(std::ostream& output, std::ostream& bed, ireadstream * input, const Database * data) {
	INFO("Create Aho-Corasick pattern automata ... ");

	std::map<std::string *, std::string *>::const_iterator it = data->get_data_iterator();
	AhoCorasick ahoCorasick;
	for (int i = 0; i < data->get_sequences_amount(); ++i) {
		ahoCorasick.addString(it->second);
		it++;
	}
	ahoCorasick.init();

	INFO("Done");

	ExactMatchJobWrapper filler(data, output, bed, ahoCorasick);
	hammer::ReadProcessor rp(omp_get_max_threads());
	rp.Run(*input, filler);
	VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

	ahoCorasick.cleanup();
}

void exactAndAlign(std::ostream& output, std::ostream& bed, ireadstream * input, const Database * data) {
	INFO("Create Aho-Corasick automata for kmers... ");

	std::map<std::string *, std::vector<std::string *> >::const_iterator it = data->get_kmer_iterator();
	AhoCorasick ahoCorasick;
	for (int i = 0; i < data->get_kmers_amount(); ++i) {
		ahoCorasick.addString(it->first);
		it++;
	}
	ahoCorasick.init();

	INFO("Done");

	ExactAndAlignJobWrapper filler(data, output, bed, ahoCorasick);
	hammer::ReadProcessor rp(omp_get_max_threads());
	rp.Run(*input, filler);
	VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

	ahoCorasick.cleanup();
}

void alignment(std::ostream& output, std::ostream& bed, ireadstream * input, const Database * data) {
	AlignmentJobWrapper filler(data, output, bed);
	hammer::ReadProcessor rp(omp_get_max_threads());
	rp.Run(*input, filler);
	VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
}
