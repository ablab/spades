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
#include "logging.hpp"

LOGGER("a.graphBuilder");

using namespace std;
using namespace __gnu_cxx;
using namespace abruijn;

//typedef hash_map< Sequence, int, HashSym<Sequence>, EqSym<Sequence> > SeqCount;
//SeqCount seqCount;
set<hash_t> goodHashes;

abruijn::Graph graph;

LoggerPtr logger(Logger::getLogger("a.graphBuilder"));
HashSym<Sequence> hashSym;
hash_t ha[MPSIZE - K + 1];
hash_t hbest[HTAKE];

/**
 * Calculates two (distinct) minimal hash-values of all k-mers in this read,
 * and puts them into goodHashes
 */
void processReadA(Seq<MPSIZE> r) {
	hashSym.kmers(Sequence(r), ha);
	for (int i = 0; i < HTAKE; i++) {
		hbest[i] = maxHash;
	}
	for (int i = 0; i + K <= MPSIZE; i++) {
		int j = HTAKE;
		while (j > 0 && ha[i] < hbest[j - 1]) {
			j--;
		}
		if (j == HTAKE || ha[i] == hbest[j]) {
			continue;
		}
		for (int k = HTAKE - 1; k > j; k--) {
			hbest[k] = hbest[k - 1];
		}
		hbest[j] = ha[i];
	}
	for (int i = 0; i < HTAKE && hbest[i] < maxHash; i++) {
		goodHashes.insert(hbest[i]);
	}
}

void selectGood() {
}

void processReadB(Seq<MPSIZE> r) {
	vector<abruijn::Vertex*> vs;
	vector<int> index;
	Sequence s(r);
	hashSym.kmers(s, ha);
	for (int i = 0; i + K <= MPSIZE; i++) {
		if (goodHashes.find(ha[i]) != goodHashes.end()) {
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
	vector<Read> *v1 = ireadstream::readAll(file_names[0], CUT);
	vector<Read> *v2 = ireadstream::readAll(file_names[1], CUT);
//	ireadstream<MPSIZE, 2> inputStream(file_names);

	INFO("Processing-A...");
	//mate_read<MPSIZE>::type mp;
	for (size_t i = 0; i < v1->size(); i++) {
//		TODO:
//		DEBUG(((*v1)[i]).isValid());
		//inputStream >> mp;
		processReadA(Seq<MPSIZE>((*v1)[i]));
		processReadA(Seq<MPSIZE>((*v2)[i]));
		VERBOSE(i, " reads read");
	}
	//inputStream.reset();
	INFO("processReadA done: " << goodHashes.size() << " earmarked hashes");

	INFO("Selecting good kmers...");
	selectGood();

	INFO("Processing-B...");
	for (size_t i = 0; i < v1->size(); i++) {
		//inputStream >> mp;
		processReadB(Seq<MPSIZE>((*v1)[i]));
		processReadB(Seq<MPSIZE>((*v2)[i]));
		VERBOSE(i, " reads processed");
	}
	//inputStream.close();
	delete v1;
	delete v2;
	INFO(graph.vertices.size() << " vertices");

//	INFO("Condensing-A graph...");
//	graph.condenseA();

	INFO("Getting statistics...");
	graph.stats();

	INFO("Outputting graph to " << OUTPUT_FILE);
	ofstream outputStream((OUTPUT_FILES + ".dot").c_str(), ios::out);
	graph.output(outputStream);
	outputStream.close();
}
