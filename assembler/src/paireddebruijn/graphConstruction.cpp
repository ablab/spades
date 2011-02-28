#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"

using namespace paired_assembler;
void constructGraph() {
	readsToPairs(parsed_reads, parsed_k_l_mers);
	pairsToSequences(parsed_k_l_mers, parsed_k_sequence);
	edgesMap edges = sequencesToMap(parsed_k_sequence);
	Graph *g = new Graph();
	createVertices(g, edges);
}

edgesMap sequencesToMap(string parsed_k_sequence) {
	FILE *inFile = fopen(parsed_k_sequence.c_str(), "r");
	vector<VertexPrototype *> prototypes;
	edgesMap res;
	prototypes.reserve(maxSeqLength);
	while (1) {
		char s[maxSeqLength];
		int size, scanf_res;
		ll kmer;
		if ((scanf_res = fscanf(inFile, "%lld %d", &kmer, &size)) != 2) {
			if (scanf_res == 0)
				cerr << "sequencesToMap finished reading";
			else {
				cerr << "sequencesToMap error in reading headers";
			}
		}
		prototypes.clear();
		forn(i, size) {
			scanf_res = fscanf(inFile, "%s", s);
			if (!scanf_res)
				cerr << "sequencesToMap error in reading sequences";
			Sequence *seq = new Sequence(s);
			VertexPrototype *v = new VertexPrototype();
			v->lower = seq;
			v->start = 0;
			v->finish = 0;
			prototypes.pb(v);
		}
		res.insert(mp(kmer, prototypes));
	}
	return res;
}

void createVertices(Graph *g, edgesMap &edges){
	for(edgesMap::iterator iter = edges.begin(); iter != edges.end(); ++iter){
		int size = iter->second.size();
		forn(i, size) {
			if (((iter->se)[i]->start == 0) && ((iter->se)[i]->finish == 0)){
				g->addVertex(new Sequence(decompress(k>>2, k - 1)), (iter->se)[i]->lower);//(lower->shiftRight)
			}
		}
	}
}

