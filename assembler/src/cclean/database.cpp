#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <exception>

#include "database.hpp"
#include "QcException.hpp"
#include "utils.hpp"
#include "io/read_processor.hpp"
#include "config_struct_cclean.hpp"
#include "valid_kmer_generator.hpp"

void DatabaseFiller::insert2db(const cclean::KMer seq, std::string * sequence) {
	std::string * kmer = new std::string(seq.str());

#pragma omp critical
	{
		std::string str = seq.str();
		if (kmer2listOfSeq->end() == kmer2listOfSeq->find(&str)) {
			std::vector<std::string * > source;
			source.push_back(sequence);
			kmer2listOfSeq->insert(std::make_pair(kmer, source));
		} else {
			(*kmer2listOfSeq)[kmer].push_back(sequence);
		}
	}
}

bool DatabaseFiller::operator()(const Read &r) {
	try {
		std::string * name = new std::string(r.getName());
		std::string * sequence = new std::string(r.getSequenceString());
		//complement strings:
		std::string * name_c = new std::string(std::string(r.getName()) + " (complementary)");
		std::string * sequence_c = new std::string(reverseComplement(*sequence));

#pragma omp critical
		{
			name2seq->insert(std::make_pair(name, sequence));
			seq2name->insert(std::make_pair(sequence, name));
			name2seq->insert(std::make_pair(name_c, sequence_c));
			seq2name->insert(std::make_pair(sequence_c, name_c));
		}
		ValidKMerGenerator<cclean::K> gen(r);
		while (gen.HasMore()) {
			cclean::KMer seq = gen.kmer();
			insert2db(seq, sequence);

			seq = !seq;
			insert2db(seq, sequence_c);

			gen.Next();
		}

	} catch (std::exception& e) {
	}
	return false;
}

Database::Database(const std::string& filename) {
	ireadstream * input = new ireadstream(filename);
	DatabaseFiller filler;
	hammer::ReadProcessor rp(std::min(omp_get_max_threads(), cclean_cfg::get().nthreads));
	rp.Run(*input, filler);
	delete input;

	name2seq = filler.getName2seq();
	seq2name = filler.getSeq2name();
	kmer2listOfSeq = filler.getKmer2listOfSeq();
}

Database::~Database() {
	for (std::map<std::string *, std::string *>::const_iterator it = name2seq->begin(); it != name2seq->end(); ++it) {
		delete it->first;
		delete it->second;
	}

	for (std::map<std::string *, std::vector<std::string *>, Compare>::const_iterator it = kmer2listOfSeq->begin(); it != kmer2listOfSeq->end(); ++it) {
		delete it->first;
	}

	delete name2seq;
	delete seq2name;
	delete kmer2listOfSeq;
}

void Database::get_sequence_by_name(const std::string& name, std::string& out_seq) const {
	std::map<std::string *, std::string *>::const_iterator it = name2seq->find(const_cast<std::string *>(&name));
	if (name2seq->end() == it) {
		throw QcException("Element not found: " + name);
	}
	out_seq.assign(*(it->second));
}

void Database::get_name_by_sequence(const std::string& seq, std::string& out_name) const {
	std::map<std::string *, std::string *>::const_iterator it = seq2name->find(const_cast<std::string *>(&seq));
	if (name2seq->end() == it) {
		throw QcException("Element not found: " + seq);
	}
	out_name.assign(*(it->second));
}

void Database::get_sequences_for_kmer(const std::string& kmer, std::vector<std::string *>& out_seq) const {
	std::map<std::string *, std::vector<std::string *>, Compare>::const_iterator it = kmer2listOfSeq->find(const_cast<std::string *>(&kmer));
	if (kmer2listOfSeq->end() != it) {
		out_seq.assign(it->second.begin(), it->second.end());
	}
}

int Database::get_sequences_amount() const {
	return name2seq->size();
}

int Database::get_kmers_amount() const {
	return kmer2listOfSeq->size();
}

std::map<std::string *, std::string *>::const_iterator Database::get_data_iterator() const {
	return name2seq->begin();
}

std::map<std::string *, std::vector<std::string *>, Compare>::const_iterator Database::get_kmer_iterator() const {
	return kmer2listOfSeq->begin();
}
