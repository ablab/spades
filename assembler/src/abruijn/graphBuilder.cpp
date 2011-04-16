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

LOGGER("a.graphBuilder");

using namespace std;
using namespace __gnu_cxx;
using namespace abruijn;

//typedef hash_map< Sequence, int, HashSym<Sequence>, EqSym<Sequence> > SeqCount;
//SeqCount seqCount;
set<hash_t> earmarkedHashes;

abruijn::Graph graph;

LoggerPtr logger(Logger::getLogger("a.graphBuilder"));
HashSym<Sequence> hashSym;
hash_t ha[MPSIZE - K + 1];
hash_t hbest[HTAKE];

bool isTrusted(hash_t hash) {
	return true;
}

/**
 * Calculates two (distinct) minimal hash-values of all k-mers in this read,
 * and puts them into goodHashes
 */
void processReadA(Sequence s) {
	hashSym.kmers(s, ha);
	for (int i = 0; i < HTAKE; i++) {
		hbest[i] = maxHash;
	}
	for (int i = 0; i + K <= MPSIZE; i++) {
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
	for (int i = 0; i < HTAKE && hbest[i] < maxHash; i++) {
		earmarkedHashes.insert(hbest[i]);
	}
}

void processReadAA(Sequence s) {
	hashSym.kmers(s, ha);
	hash_t he = maxHash;
	hash_t hb = maxHash;
	for (int i = 0; i + K <= MPSIZE; i++) {
		hash_t hi = ha[i];
		if (!isTrusted(hi)) {
			continue;
		}
		if (earmarkedHashes.find(hi) != earmarkedHashes.end()) {
			if (he == maxHash) {
				he = hi;
			} else {
				return;
			}
		} else if (hi < hb) {
			hb = hi;
		}
	}
	earmarkedHashes.insert(hb);
}

void processReadE(Sequence s) {
	vector<abruijn::Vertex*> vs;
	vector<int> index;
	hashSym.kmers(s, ha);
	for (int i = 0; i + K <= MPSIZE; i++) {
		if (earmarkedHashes.find(ha[i]) != earmarkedHashes.end()) {
			vs.push_back(graph.getVertex(s.Subseq(i, i + K)));
			index.push_back(i);
		}
	}
	for (size_t i = 0; i + 1 < vs.size(); i++) {
		graph.addEdge(vs[i], vs[i + 1], s.Subseq(index[i], index[i + 1] + K));
	}
}

void GraphBuilder::build() {
	std::string file_names[2] = {INPUT_FILES};
//	vector<Read> *v1 = ireadstream::readAll(file_names[0], CUT);

	StrobeReader<2, Read, ireadstream> sr(file_names);
	INFO("Processing-A...");
	vector<Read> v;
	for (size_t i = 0; !sr.eof() && i < CUT; i++ ) {
		sr >> v;
		processReadA(v[0].getSequence());
		processReadA(v[1].getSequence());
		VERBOSE(i, " reads");
	}
	INFO("Done: " << earmarkedHashes.size() << " earmarked hashes");

	if (HTAKE == 1) {
		sr.reset();
		INFO("Processing-AA...");
		for (size_t i = 0; !sr.eof() && i < CUT; i++ ) {
			sr >> v;
			processReadAA(v[0].getSequence());
			processReadAA(v[1].getSequence());
			VERBOSE(i, " reads");
		}
		INFO("Done: " << earmarkedHashes.size() << " earmarked hashes");
	}

	sr.reset();
	INFO("Processing-E...");
	for (size_t i = 0; !sr.eof() && i < CUT; i++ ) {
		sr >> v;
		processReadE(v[0].getSequence());
		processReadE(v[1].getSequence());
		VERBOSE(i, " reads");
	}
	INFO("Done: " << graph.vertices.size() << " vertices");

//	INFO("Processing-E...");
//	for (size_t i = 0; i < v1->size(); i++) {
//		processReadE(Seq<MPSIZE>((*v1)[i]));
//		processReadE(Seq<MPSIZE>((*v2)[i]));
//		VERBOSE(i, " reads processed");
//	}
//	INFO(graph.vertices.size() << " vertices");

//	INFO("Condensing-A graph...");
//	graph.condenseA();

	INFO("Getting statistics...");
	graph.stats();

	INFO("Outputting graph to " << OUTPUT_FILE);
	ofstream outputStream((OUTPUT_FILES + ".dot").c_str(), ios::out);
	graph.output(outputStream);
	outputStream.close();
}
