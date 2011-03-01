#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"

#define similar(a,b,c) false

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
			v->used = 0;
			prototypes.pb(v);
		}
		res.insert(mp(kmer, prototypes));
	}
	return res;
}

void createVertices(Graph *g, edgesMap &edges) {
	vertecesMap verts;
	for (edgesMap::iterator iter = edges.begin(); iter != edges.end();) {
		int size = iter->second.size();
		forn(i, size) {
			ll kmer = iter->fi;
			if ((!(iter->se)[i]->used)) {
				ll finishKmer = kmer & (~((ll) 3 << (2 * (k - 1))));
				Sequence *finishSeq = new Sequence((iter->se)[i]->lower->Str());
				ll startKmer = kmer >> 2;
				expandDown(edges, verts, finishKmer, finishSeq);
			}
		}
		edges.erase(iter++);
	}
}

void expandDown(edgesMap &edges, vertecesMap &verts, ll finishKmer,
		Sequence *finishSeq) {
	while (1) {
		vertecesMap::iterator iter = verts.find(finishKmer);
		if (iter != verts.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (similar(finishSeq, (iter->se)[i]->lower, k))
					return;
			}
		}
		if (!CheckUnuqueWayUp(edges, finishKmer, finishSeq))
			return;
		if (!GoUnuqueWayDown(edges, finishKmer, finishSeq))
			return;
	}
}

int CheckUnuqueWayUp(edgesMap &edges, ll finishKmer, Sequence *finishSeq) {
	int count = 0;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl << (2 * (k - 1)) | finishKmer;
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (similar(finishSeq, (iter->se)[i]->lower, k)) {
					count++;
					if (count > 1)
						return 0;
				}
			}
		}

	}
	return count;
}

int GoUnuqueWayUp(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq) {
	int count = 0;
	Sequence *PossibleSequence;
	ll PossibleKmer;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl << (2 * (k - 1)) | finishKmer;
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (similar(finishSeq, (iter->se)[i]->lower, k)) {
					count++;
					if (count > 1)
						return 0;
					PossibleKmer = tmpKmer;
					PossibleSequence = (iter->se)[i]->lower;

				}
			}
		}
	}
	if (count == 1) {
		finishKmer = PossibleKmer >> 2;
		finishSeq = PossibleSequence;
		return 1;
	}
	return 0;
}

int GoUnuqueWayDown(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq) {
	int count = 0;
	Sequence *PossibleSequence;
	ll PossibleKmer;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl | (finishKmer << (2));
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (similar(finishSeq, (iter->se)[i]->lower, k)) {
					count++;
					if (count > 1)
						return 0;
					PossibleKmer = tmpKmer;
					PossibleSequence = (iter->se)[i]->lower;

				}
			}
		}
	}
	if (count == 1) {
		finishKmer = (PossibleKmer) & (~(((ll) 3) << (2 * (k - 1))));
		finishSeq = PossibleSequence;
		return 1;
	}
	else return 0;
}


int CheckUnuqueWayDown(edgesMap &edges, ll finishKmer, Sequence* finishSeq) {
	int count = 0;
	Sequence *PossibleSequence;
	ll PossibleKmer;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl | (finishKmer << (2));
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (similar(finishSeq, (iter->se)[i]->lower, k)) {
					count++;
					if (count > 1)
						return 0;
					PossibleKmer = tmpKmer;
					PossibleSequence = (iter->se)[i]->lower;

				}
			}
		}
	}
	if (count == 1) {
		return 1;
	}
	else return 0;
}

