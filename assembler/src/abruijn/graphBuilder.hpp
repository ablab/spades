#ifndef GRAPHBUILDER_H_
#define GRAPHBUILDER_H_

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <set>
#include "hash.hpp"
#include "parameters.hpp"
#include "logging.hpp"
#include "abruijngraph.hpp"
#include "ireadstream.hpp"

namespace abruijn {

using namespace std;
using namespace __gnu_cxx;
using hashing::hash_t;

class GraphBuilder
{
	//typedef hash_map< Sequence, int, HashSym<Sequence>, EqSym<Sequence> > SeqCount;
	//SeqCount seqCount;
	/**
	 * Stores bitmask:
	 * 1 = the lexicographically lesser k-mer has been seen as non-rightmost k-mer in read
	 * 2 = the greater k-mer ...
	 */
	hashing::HashSym<Sequence> hashSym;
	typedef vector<hash_t> hash_vector;
	hash_vector ha;

public:
	void findMinimizers(Sequence s);
	void findLocalMinimizers(Sequence s, size_t window_size);
	void findSecondMinimizer(Sequence s);
	void revealTips(Sequence s);
	void findTipExtensions(Sequence s);
	void lookRight(Sequence s);
	void addToGraph(Sequence s);

	map<hash_t, char> has_right;
	map<hash_t, char> tips;
	set<hash_t> earmarked_hashes;
	map<hash_t, set<hash_t> > tip_extensions;
	hash_vector hbest;
	abruijn::Graph graph_;
	size_t htake_;
};

template <typename Reader>
class GraphBuilderMaster {
	Reader reader_;
	int mode_;
	GraphBuilder gb_;
public:
	GraphBuilderMaster(Reader reader, size_t htake, int mode) :
		reader_(reader), mode_(mode), gb_() {
		gb_.hbest.reserve(htake);
		gb_.htake_ = htake;
	}

	void build() {
		Read r;

		INFO("===== Finding " << gb_.htake_ << " minimizers in each read... =====");
		reader_.reset();
		for (size_t i = 0; !reader_.eof(); ++i) {
			reader_ >> r;
			if (mode_ & 1) {
				gb_.findMinimizers(r.getSequence());
			} else {
				gb_.findLocalMinimizers(r.getSequence(), 51);
			}
			VERBOSE(i, " single reads");
		}
		INFO("Done: " << gb_.earmarked_hashes.size() << " earmarked hashes");

		if ((mode_ & 2) && (gb_.htake_ == 1)) {
			INFO("===== Finding second minimizers... =====");
			reader_.reset();
			for (size_t i = 0; !reader_.eof(); ++i) {
				reader_ >> r;
				gb_.findSecondMinimizer(r.getSequence());
				VERBOSE(i, " single reads");
			}
			INFO("Done: " << gb_.earmarked_hashes.size() << " earmarked hashes");
		}

		for(;;) {
			size_t eh = gb_.earmarked_hashes.size();
			gb_.has_right.clear();
			gb_.tips.clear();
			gb_.tip_extensions.clear();

			INFO("===== Revealing tips... =====");
			reader_.reset();
			for (size_t i = 0; !reader_.eof(); ++i) {
				reader_ >> r;
				gb_.revealTips(r.getSequence());
				VERBOSE(i, " single reads");
			}
			for (map<hash_t, char>::iterator it = gb_.has_right.begin(); it != gb_.has_right.end(); ++it) {
				if (it->second != 3) {
					TRACE(it->first << " " << (int) it->second);
					gb_.tips.insert(*it);
				}
			}
			gb_.has_right.clear();
			INFO("Done: " << gb_.tips.size() << " tips.");

			INFO("===== Finding tip extensions... =====");
			reader_.reset();
			for (size_t i = 0; !reader_.eof(); ++i) {
				reader_ >> r;
				gb_.findTipExtensions(r.getSequence());
				VERBOSE(i, " single reads");
			}
			INFO("Done: " << gb_.has_right.size() << " possible tip extensions");

			INFO("===== Looking to the right... =====");
			reader_.reset();
			for (size_t i = 0; !reader_.eof(); ++i) {
				reader_ >> r;
				gb_.lookRight(r.getSequence());
				VERBOSE(i, " single reads");
			}
			for (map<hash_t, set<hash_t> >::iterator it = gb_.tip_extensions.begin(); it != gb_.tip_extensions.end(); ++it) {
				bool ok = false;
				hash_t found = hashing::kMax;
				for (set<hash_t>::iterator ext = it->second.begin(); ext != it->second.end(); ext++) {
					if (gb_.earmarked_hashes.count(*ext)) {
						ok = true;
						break;
					}
					if (gb_.has_right[*ext] == 3 && found == hashing::kMax) {
						found = *ext;
					}
				}
				if (ok) {
					continue;
				}
				gb_.earmarked_hashes.insert(found);
			}
			INFO("Done: " << eh << " -> " << gb_.earmarked_hashes.size() << " earmarked hashes");
			if (eh == gb_.earmarked_hashes.size()) {
				break;
			}
		}

		INFO("===== Adding reads to graph as paths... =====");
		reader_.reset();
		for (size_t i = 0; !reader_.eof(); ++i) {
			reader_ >> r;
			gb_.addToGraph(r.getSequence());
			VERBOSE(i, " single reads");
		}
		INFO("Done: " << gb_.graph_.vertices.size() << " vertices");

		INFO("===== Condensing-A graph... =====");
		gb_.graph_.Condense();
		INFO("Done: " << gb_.graph_.vertices.size() << " vertices");

		return;
	}

	abruijn::Graph* graph() {
		return &gb_.graph_;
	}
};

}

#endif /* GRAPHBUILDER_H_ */
