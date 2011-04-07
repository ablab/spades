#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"
#include "vardConstruction.hpp"
#include "constdConstruction.hpp"

LOGGER("p.graphConstruction");


using namespace paired_assembler;


void constructGraph(PairedGraph &graph) {
	INFO ("Read edges");
	edgesMap edges = sequencesToMap(parsed_k_sequence, true);
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


edgesMap sequencesToMap(string parsed_k_sequence, bool usePaired) {
	FILE *inFile = fopen(parsed_k_sequence.c_str(), "r");

	vector<EdgePrototype *> prototypes;
	edgesMap res;
	prototypes.reserve(maxSeqLength);
	int count = 0;
	while (1) {
		char s[maxSeqLength];
		int coverage;
		int size, scanf_res;
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
		forn(i, size) {
			scanf_res = fscanf(inFile, "%s %d", s, &coverage);
			if (!scanf_res) {
				cerr << "sequencesToMap error in reading sequences";
			}
			//			cerr <<s;
			Sequence *seq;
			if (usePaired) {
				seq = new Sequence(s);
			} else
				seq = new Sequence(auxilary_lmer);
			EdgePrototype *v = new EdgePrototype(seq, 0);
			v->coverage = coverage;
			if (usePaired || !i)
				prototypes.pb(v);
		}
		res.insert(mp(kmer, prototypes));
	}
	return res;
}
