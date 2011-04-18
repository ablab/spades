#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"
#include "vardConstruction.hpp"
#include "constdConstruction.hpp"

LOGGER("p.graphConstruction");


using namespace paired_assembler;


void constructGraph(PairedGraph &graph) {
	INFO ("Read edges");
	edgesMap edges = sequencesToMap(parsed_k_sequence);
	appendLmers(parsed_l_mers, edges);
	graph.VertexCount = 0;
	if (distance_type == "const")
		constd::createVertices(edges, graph);
	else {
		vard::createVertices(edges, graph);
		vard::clearUseOfEdgePrototypes(edges);
		vard::createEdges(edges, graph, true);

	}
		//assert(0);
	INFO ("End create vertices");
}


edgesMap sequencesToMap(string parsed_k_sequence) {
	FILE *inFile = fopen(parsed_k_sequence.c_str(), "r");

	edgesMap res;
	vector<EdgePrototype *> prototypes;
	prototypes.reserve(maxSeqLength);
	int count = 0;
	while (1) {
		int size, scanf_res;
		char s[maxSeqLength];
		int coverage;

		ll kmer;
		count++;
		if (!(count & ((1 << 16) - 1))) {
			cerr << count << "k-seq readed" << endl;
		}
		scanf_res = fscanf(inFile, "%lld %d", &kmer, &size);
		//		cerr<<scanf_res;
		if ((scanf_res) != 2) {

			if (scanf_res == -1) {
				cerr << "sequencesToMap finished reading";
				break;
			} else {
				cerr << "sequencesToMap error in reading headers";
				continue;
			}
		}
		prototypes.clear();
		if (size == 0) {
			Sequence *seq;
			seq = new Sequence("");
		}
		else
		forn(i, size) {
			scanf_res = fscanf(inFile, "%s %d", s, &coverage);
			if (!scanf_res) {
				cerr << "sequencesToMap error in reading sequences";
			}
			//			cerr <<s;
			Sequence *seq;
			seq = new Sequence(s);
			EdgePrototype *v = new EdgePrototype(seq, 0);
			v->coverage = coverage;
			if (!i)
				prototypes.pb(v);
		}

		if (size>0) res.insert(mp(kmer, prototypes));
	}
	return res;
}

void appendLmers(string parsed_l_mers, edgesMap &edges) {
	FILE *inFile = fopen(parsed_l_mers.c_str(), "r");
	vector<EdgePrototype *> prototypes;
	prototypes.resize(1);
	prototypes[0] = NULL;
	int count = 0;

	while (1) {
		char s[maxSeqLength];
		int coverage;
		int size, scanf_res;
		ll kmer;
		count++;
		if (!(count & ((1 << 16) - 1))) {
			cerr << count << "lmer - readed" << endl;
		}
		scanf_res = fscanf(inFile, "%lld %d", &kmer, &coverage);
//		cerr << kmer << " readed "<< endl;
		//		cerr<<scanf_res;
		if ((scanf_res) != 2) {

			if (scanf_res <= 0) {
				cerr << "appendLmers finished reading";
				break;
			} else {
				cerr << "appendLmers error in reading";
				continue;
			}
		}
		if (edges.find(kmer) == edges.end()) {
			Sequence *seq = new Sequence("A");
			EdgePrototype *v = new EdgePrototype(seq, 0);
			v->coverage = 1 + coverage * range_variating * 1.5;
			prototypes[0] = v;
			cerr <<" kmer inserting: " << kmer << endl;
			edges.insert(mp(kmer, prototypes));
		}
	}
}
