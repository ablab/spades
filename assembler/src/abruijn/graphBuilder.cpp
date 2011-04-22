#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <set>
#include "abruijngraph.hpp"
#include "hash.hpp"
#include "graphBuilder.hpp"
#include "parameters.hpp"
#include <ext/hash_map>
#include "strobe_reader.hpp"
#include "logging.hpp"

using namespace std;
using namespace __gnu_cxx;

namespace abruijn {

//typedef hash_map< Sequence, int, HashSym<Sequence>, EqSym<Sequence> > SeqCount;
//SeqCount seqCount;
/**
 * Stores bitmask:
 * 1 = the lexicographically lesser k-mer has been seen as non-rightmost k-mer in read
 * 2 = the greater k-mer ...
 */
map<hash_t, char> has_right;
map<hash_t, char> tips;
set<hash_t> earmarked_hashes;
map<hash_t, set<hash_t> > tip_extensions;

abruijn::Graph graph;

HashSym<Sequence> hashSym;
vector<hash_t> ha;
hash_t hbest[HTAKE];

bool isTrusted(hash_t hash) {
	return true;
}

/**
 * Calculates HTAKE (distinct) minimal hash-values of all k-mers in this read,
 * and puts them into earmarked_hashes
 */
void findMinimizers(Sequence s) {
	ha.reserve(s.size());
	hashSym.kmers(s, ha);
	for (size_t i = 0; i < HTAKE; i++) {
		hbest[i] = maxHash;
	}
	for (size_t i = 0; i + K <= s.size(); i++) {
		hash_t hi = ha[i];
		if (!isTrusted(hi)) {
			continue;
		}
		int j = HTAKE;
		while (j > 0 && hi < hbest[j - 1]) {
			j--;
		}
		if (j == HTAKE || hi == hbest[j]) {
			continue;
		}
		for (int k = HTAKE - 1; k > j; k--) {
			hbest[k] = hbest[k - 1];
		}
		hbest[j] = hi;
	}
	for (size_t i = 0; i < HTAKE && hbest[i] < maxHash; i++) {
		earmarked_hashes.insert(hbest[i]);
	}
}

/**
 * If only one k-mer from s is earmarked,
 * this method earmarks 1 more k-mer (with second minimal hash-value)
 */
void findSecondMinimizer(Sequence s) {
	hashSym.kmers(s, ha);
	hash_t he = maxHash;
	hash_t hb = maxHash;
	for (size_t i = 0; i + K <= s.size(); i++) {
		hash_t hi = ha[i];
		if (earmarked_hashes.count(hi)) {
			if (he == maxHash) {
				he = hi;
			} else {
				return;
			}
		} else if (hi < hb) {
			if (!isTrusted(hi)) {
				continue;
			}
			hb = hi;
		}
	}
	assert(hb < maxHash);
	earmarked_hashes.insert(hb);
}

bool Lesser(const Sequence& s) {
	return s < !s;
}

bool LesserK(const Sequence& s, size_t i) {
	return Lesser(s.Subseq(i, i + K));
}

void revealTips(Sequence s) {
	hashSym.kmers(s, ha);
	vector<int> index;
	for (size_t i = 0; i + K <= s.size(); i++) {
		hash_t hi = ha[i];
		if (earmarked_hashes.count(hi)) {
			index.push_back(i);
		}
	}
	for (size_t i = 0; i < index.size(); ++i) {
		hash_t hi = ha[index[i]];
		bool lesser = LesserK(s, index[i]);
		if (i < index.size() - 1) {
			has_right[hi] |= lesser ? 1 : 2;
		}
		if (i > 0) {
			has_right[hi] |= lesser ? 2 : 1;
		}
	}
}

void findTipExtensions(Sequence s) {
	hashSym.kmers(s, ha);
	vector<int> index;
	for (size_t i = 0; i + K <= s.size(); i++) {
		hash_t hi = ha[i];
		if (tips.count(hi)) {
			index.push_back(i);
		}
	}
	for (size_t i = 0; i < index.size(); i++) {
		bool lesser = LesserK(s, index[i]);
		int l = lesser ? 1 : 2;
		hash_t hi = ha[index[i]];
		TRACE(index[i] << " " << l << " " << (int) tips[hi]);
		size_t low, high;
		if (tips[hi] == l) {
			low = 0;
			high = index[i];
		} else {
			low = index[i] + 1;
			high = s.size() + 1 - K;
		}
		int x = 0;
		for (size_t j = low; j < high; j++) {
			hash_t hj = ha[j];
			if (!isTrusted(hj)) {
				continue;
			}
			assert(!earmarked_hashes.count(hj));
			has_right[hj] = 0;
			tip_extensions[hi].insert(hj);
			x++;
		}
		TRACE(x << " possible continuations");
	}
}

void lookRight(Sequence s) {
	hashSym.kmers(s, ha);
	vector<size_t> index;
	for (size_t i = 0; i + K <= s.size(); i++) {
		hash_t hi = ha[i];
		if (has_right.count(hi)) {
			index.push_back(i);
		}
	}
	if (index.size() == 0) {
		return;
	}
	size_t low = 0;
	size_t high = s.size() - K;
	while (!earmarked_hashes.count(ha[low])) low++;
	while (!earmarked_hashes.count(ha[high])) high--;
	for (size_t i = 0; i < index.size(); ++i) {
		hash_t hi = ha[index[i]];
		bool lesser = LesserK(s, i);
		if (index[i] < high) {
			has_right[hi] |= lesser ? 1 : 2;
		}
		if (index[i] > low) {
			has_right[hi] |= lesser ? 2 : 1;
		}
	}
}

//void extendTip(Sequence s) {
//
//}

void addToGraph(Sequence s) {
	vector<abruijn::Vertex*> vs;
	vector<size_t> index;
	hashSym.kmers(s, ha);
	for (size_t i = 0; i + K <= s.size(); i++) {
		if (earmarked_hashes.find(ha[i]) != earmarked_hashes.end()) {
			vs.push_back(graph.getVertex(s.Subseq(i, i + K)));
			index.push_back(i);
		}
	}
	for (size_t i = 0; i + 1 < vs.size(); ++i) {
		graph.addEdge(vs[i], vs[i + 1], s.Subseq(index[i], index[i + 1] + K));
	}
}

void GraphBuilder::build() {
	std::string file_names[2] = {INPUT_FILES};
//	vector<Read> *v1 = ireadstream::readAll(file_names[0], CUT);

	StrobeReader<2, Read, ireadstream> sr(file_names);
	SimpleReaderWrapper<StrobeReader<2, Read, ireadstream> > srw(sr);
	vector<Read> v;
	Read r;
#ifdef CUT
	size_t cut = CUT;
	size_t cut2 = cut2;
#endif
#ifndef CUT
	size_t cut = -1;
	size_t cut2 = -1;
#endif

	INFO("===== Finding " << HTAKE << " minimizers in each read... =====");
	srw.reset();
	for (size_t i = 0; !srw.eof() && i < cut2; ++i) {
		srw >> r;
		findMinimizers(r.getSequence());
		VERBOSE(i, " single reads");
	}
	INFO("Done: " << earmarked_hashes.size() << " earmarked hashes");

//	if (HTAKE == 1) {
//		INFO("===== Finding second minimizers... =====");
//		srw.reset();
//		for (size_t i = 0; !srw.eof() && i < cut2; ++i) {
//			srw >> r;
//			findSecondMinimizer(r.getSequence());
//			VERBOSE(i, " single reads");
//		}
//		INFO("Done: " << earmarked_hashes.size() << " earmarked hashes");
//	}

for(;;) {
	size_t eh = earmarked_hashes.size();
	has_right.clear();
	tips.clear();
	tip_extensions.clear();

	INFO("===== Revealing tips... =====");
	srw.reset();
	for (size_t i = 0; !srw.eof() && i < cut2; ++i) {
		srw >> r;
		revealTips(r.getSequence());
		VERBOSE(i, " single reads");
	}
	for (map<hash_t, char>::iterator it = has_right.begin(); it != has_right.end(); ++it) {
		if (it->second != 3) {
			TRACE(it->first << " " << (int) it->second);
			tips.insert(*it);
		}
	}
	has_right.clear();
	INFO("Done: " << tips.size() << " tips.");

	INFO("===== Finding tip extensions... =====");
	srw.reset();
	for (size_t i = 0; !srw.eof() && i < cut2; ++i) {
		srw >> r;
		findTipExtensions(r.getSequence());
		VERBOSE(i, " single reads");
	}
	INFO("Done: " << has_right.size() << " possible tip extensions");

	INFO("===== Looking to the right... =====");
	srw.reset();
	for (size_t i = 0; !srw.eof() && i < cut2; ++i) {
		srw >> r;
		lookRight(r.getSequence());
		VERBOSE(i, " single reads");
	}
	for (map<hash_t, set<hash_t> >::iterator it = tip_extensions.begin(); it != tip_extensions.end(); ++it) {
		bool ok = false;
		hash_t found = maxHash;
		for (set<hash_t>::iterator ext = it->second.begin(); ext != it->second.end(); ext++) {
			if (earmarked_hashes.count(*ext)) {
				ok = true;
				break;
			}
			if (has_right[*ext] == 3 && found == maxHash) {
				found = *ext;
			}
		}
		if (ok) {
			continue;
		}
		earmarked_hashes.insert(found);
	}
	INFO("Done: " << eh << " -> " << earmarked_hashes.size() << " earmarked hashes");
	if (eh == earmarked_hashes.size()) {
		break;
	}
}

	INFO("Adding reads to graph as paths...");
	srw.reset();
	for (size_t i = 0; !srw.eof() && i < cut2; ++i) {
		srw >> r;
		addToGraph(r.getSequence());
		VERBOSE(i, " single reads");
	}
	INFO("Done: " << graph.vertices.size() << " vertices");

//	INFO("Condensing-A graph...");
//	graph.condenseA();

	INFO("Getting statistics...");
	graph.stats();

	INFO("Outputting graph to " << OUTPUT_FILE);
	ofstream outputStream((OUTPUT_FILES + ".dot").c_str(), ios::out);
	graph.output(outputStream);
	outputStream.close();
}

}
