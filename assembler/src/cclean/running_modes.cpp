#include "running_modes.hpp"
#include "QcException.hpp"
#include "aho_corasick.hpp"
#include "output.hpp"
#include "io/read_processor.hpp"
#include "job_wrappers.hpp"
#include "logger/log_writers.hpp"

void getDbAhoCorasick(const Database * data, AhoCorasick& ahoCorasick) {
	INFO("Create Aho-Corasick automata ... ");
	auto it = data->get_data_iterator();
	const int amount = data->get_sequences_amount();

	for (int i = 0; i < amount; ++i) {
		ahoCorasick.addString(it->second);
		it++;
	}
	ahoCorasick.init();

	INFO("Done");
}

void getKmersAhoCorasick(const Database * data, AhoCorasick& ahoCorasick) {
	INFO("Create Aho-Corasick automata for kmers... ");

	auto it = data->get_kmer_iterator();
	const int amount = data->get_kmers_amount();

	for (int i = 0; i < amount; ++i) {
		ahoCorasick.addString(it->first);
		it++;
	}
	ahoCorasick.init();

	INFO("Done");
}

void exactMatch(std::ostream& output, std::ostream& bed, ireadstream * input, const Database * data) {
	AhoCorasick ahoCorasick;
	getDbAhoCorasick(data, ahoCorasick);

	ExactMatchJobWrapper filler(data, output, bed, ahoCorasick);
	hammer::ReadProcessor rp(cfg::get().nthreads);
	rp.Run(*input, filler);
	VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

	ahoCorasick.cleanup();
}

void exactAndAlign(std::ostream& output, std::ostream& bed, ireadstream * input, const Database * data) {
	AhoCorasick dbAhoCorasick, kmersAhoCorasick;
	getDbAhoCorasick(data, dbAhoCorasick);
	getKmersAhoCorasick(data, kmersAhoCorasick);

	ExactAndAlignJobWrapper filler(data, output, bed, dbAhoCorasick, kmersAhoCorasick);
	hammer::ReadProcessor rp(cfg::get().nthreads);
	rp.Run(*input, filler);
	VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

	dbAhoCorasick.cleanup();
	kmersAhoCorasick.cleanup();
}

void alignment(std::ostream& output, std::ostream& bed, ireadstream * input, const Database * data) {
	AlignmentJobWrapper filler(data, output, bed);
	hammer::ReadProcessor rp(cfg::get().nthreads);
	rp.Run(*input, filler);
	VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
}
