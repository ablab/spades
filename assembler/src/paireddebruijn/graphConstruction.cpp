#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"
#include "../seq.hpp"
#include "../graphVisualizer.hpp"

int VertexCount;

using namespace paired_assembler;

void constructGraph() {
	//	readsToPairs(parsed_reads, parsed_k_l_mers);
	//	pairsToSequences(parsed_k_l_mers, parsed_k_sequence);
	cerr << "Read edges" << endl;
	edgesMap edges = sequencesToMap(parsed_k_sequence, false);
	cerr << "go to graph" << endl;
	gvis::GraphScheme<int> g("Paired");
	cerr << "Start vertices" << endl;
	VertexCount = 0;
	createVertices(g, edges);
}

edgesMap sequencesToMap(string parsed_k_sequence, bool usePaired) {
	FILE *inFile = fopen(parsed_k_sequence.c_str(), "r");

	vector<VertexPrototype *> prototypes;
	edgesMap res;
	prototypes.reserve(maxSeqLength);
	while (1) {
		char s[maxSeqLength];
		int size, scanf_res;
		ll kmer;
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
			scanf_res = fscanf(inFile, "%s", s);
			if (!scanf_res) {
				cerr << "sequencesToMap error in reading sequences";
			}
			//			cerr <<s;
			Sequence *seq;
			if (usePaired) {
				seq = new Sequence(s);
			} else
				seq = new Sequence(auxilary_lmer);
			//			cerr <<"seq = "<<seq->size(); //Magic
			//			cerr <<endl; //Magic
			VertexPrototype *v = new VertexPrototype();
			v->lower = seq;
			v->start = 0;
			v->finish = 0;
			v->used = 0;
			if (usePaired || !i)
				prototypes.pb(v);
		}
		res.insert(mp(kmer, prototypes));
	}
	return res;
}

void createVertices(gvis::GraphScheme<int> &g, edgesMap &edges) {
	int inD[MAX_VERT_NUMBER], outD[MAX_VERT_NUMBER];
	char Buffer[20];
	vertecesMap verts;
	cerr << "Start createVertices " << edges.size() << endl;
	forn(i,MAX_VERT_NUMBER) {
		inD[i] = 0;
		outD[i] = 0;
	}

	for (edgesMap::iterator iter = edges.begin(); iter != edges.end();) {
		int size = iter->second.size();
		ll kmer = iter->fi;
		//cerr << "kmer " << kmer << " " << decompress(kmer, k) << " Pairs "		<< size << endl;
		forn(i, size) {
			if ((!(iter->se)[i]->used)) {
				int length = 1;
				(iter->se)[i]->used = 1;
				cerr << "seq " << (iter->se)[i]->lower->Str() << endl;
				//	cerr<<(iter->se)[i]->lower->Str()<<" "<<(iter->se)[i]->lower->size();
				ll finishKmer = kmer & (~((ll) 3 << (2 * (k - 1))));
				Sequence *finishSeq = new Sequence(
						(iter->se)[i]->lower->Subseq(1,
								(iter->se)[i]->lower->size()));
				;
				ll startKmer = kmer >> 2;
				Sequence *startSeq = new Sequence(
						(iter->se)[i]->lower->Subseq(0,
								(iter->se)[i]->lower->size() - 1));
				cerr << "expandRight " << finishKmer << " " << finishSeq->Str()
						<< endl;
				length += expandRight(edges, verts, finishKmer, finishSeq);
				int toVert = storeVertex(g, verts, finishKmer, finishSeq);

				//				if (toVert == VertexCount - 1)
	//			g.addVertex(toVert, decompress(finishKmer, k - 1));

				cerr << "expandLeft " << startKmer << " " << startSeq->Str()
						<< endl;
				length += expandLeft(edges, verts, startKmer, startSeq);
				int fromVert = storeVertex(g, verts, startKmer, startSeq);
//				if (fromVert == VertexCount - 1)
	//				g.addVertex(fromVert, decompress(startKmer, k - 1));
				cerr << "from " << fromVert << " to " << toVert << " length "
						<< length << endl;
				if (fromVert - toVert != 1)
					cerr << "GOOD" << endl;
				sprintf(Buffer, "%d", length);
				g.addEdge(fromVert, toVert, Buffer);
				outD[fromVert]++;
				inD[toVert]++;

			}
		}
		//		edges.erase(iter++);
		iter++;
	}
	g.output();
	forn(i, VertexCount) {
		cerr << i << " +" << inD[i] << " -" << outD[i] << endl;
	}
}

int expandRight(edgesMap &edges, vertecesMap &verts, ll &finishKmer,
		Sequence* &finishSeq) {
	int length = 0;
	vertecesMap::iterator iter;
	while (1) {
		cerr << "expandRight: process " << decompress(finishKmer, k - 1) << " "
				<< finishSeq->Str() << endl;
		iter = verts.find(finishKmer);
		if (iter != verts.end()) {

			int size = iter->second.size();
			//			cerr<<"expandRight Such kMer exist in vertices"<<size<< " find "<< finishSeq->Str()<< endl;
			forn(i, size) {
				//				cerr<<"posible "<<((iter->se)[i]->lower)->Str()<<endl;
				if (finishSeq->similar(*((iter->se)[i]->lower), l - 1, 0)) {
					cerr << "expandRight: vertex presented "
							<< (iter->se)[i]->start << endl;
					return length;
				}
			}
		}

		if (!checkUniqueWayLeft(edges, finishKmer, finishSeq)) {
			cerr << "expandRight: way Left not unique" << endl;
			return length;
		}

		if (!goUniqueWayRight(edges, finishKmer, finishSeq)) {
			cerr << "expandRight: way Right not unique" << endl;
			return length;
		} else
			length++;
	}
}

int expandLeft(edgesMap &edges, vertecesMap &verts, ll &startKmer,
		Sequence* &startSeq) {
	int length = 0;

	while (1) {
		cerr << "expandLeft: process " << decompress(startKmer, k - 1) << " "
				<< startSeq->Str() << endl;
		vertecesMap::iterator iter = verts.find(startKmer);
		if (iter != verts.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (startSeq->similar(*((iter->se)[i]->lower), l - 1, 0)) {
					cerr << "expandLeft: vertex presented "
							<< (iter->se)[i]->start << endl;

					return length;
				}
			}
		}
		if (!checkUniqueWayRight(edges, startKmer, startSeq)) {
			cerr << "expandLeft: way Right not unique" << endl;
			return length;
		}
		if (!goUniqueWayLeft(edges, startKmer, startSeq)) {
			cerr << "expandLeft: way Left not unique" << endl;
			return length;
		} else
			length++;
	}
}

int checkUniqueWayLeft(edgesMap &edges, ll finishKmer, Sequence *finishSeq) {
	int count = 0;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl << (2 * (k - 1)) | finishKmer;
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (finishSeq->similar(*(iter->se)[i]->lower, l - 1, -1)) {
					count++;
					//					cerr<<"posible way Left "<<count<<" "<<decompress(tmpKmer,k) <<" "<<((iter->se)[i]->lower)->Str()<< " similar for "<<finishSeq->Str()<<" with dir -1 and param "<<l-1<< endl;
					if (count > 1) {
						cerr << "count > 1" << endl;
						return 0;
					}
				}
			}
		}
	}
	if (count == 0)
		cerr << "checkUniqueWayLeft count = 0" << endl;
	return count;
}

int goUniqueWayLeft(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq) {
	int count = 0;
	Sequence *PossibleSequence;
	ll PossibleKmer;
	int seqIndex;
	edgesMap::iterator PossibleIter;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl << (2 * (k - 1)) | finishKmer;
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (finishSeq->similar(*((iter->se)[i]->lower), l - 1, -1)) {
					count++;
					if (count > 1)
						return 0;
					PossibleKmer = tmpKmer;
					PossibleSequence = (iter->se)[i]->lower;
					seqIndex = i;
					PossibleIter = iter;

				}
			}
		}
	}
	if (count == 1) {
		finishKmer = PossibleKmer >> 2;
		finishSeq = new Sequence(
				(PossibleIter->se)[seqIndex]->lower->Subseq(0,
						(PossibleIter->se)[seqIndex]->lower->size() - 1));//PossibleSequence;
		(PossibleIter->se)[seqIndex]->used = 1;
		return 1;
	}
	cerr << "goUniqueWayLeft: no way Left exist" << endl;
	return 0;
}

int goUniqueWayRight(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq) {
	int count = 0;
	Sequence *PossibleSequence;
	ll PossibleKmer;
	int seqIndex;
	edgesMap::iterator PossibleIter;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl | (finishKmer << (2));
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (finishSeq->similar(*((iter->se)[i]->lower), l - 1, 1)) {
					if ((iter->se)[i]->used)
						return 0;
					count++;
					if (count > 1)
						return 0;
					PossibleKmer = tmpKmer;
					PossibleSequence = (iter->se)[i]->lower;
					seqIndex = i;
					PossibleIter = iter;

				}
			}
		}
	}
	if (count == 1) {
		finishKmer = (PossibleKmer) & (~(((ll) 3) << (2 * (k - 1))));
		finishSeq = new Sequence(
				(PossibleIter->se)[seqIndex]->lower->Subseq(1,
						(PossibleIter->se)[seqIndex]->lower->size()));//PossibleSequence;
		//		finishSeq = PossibleSequence;
		(PossibleIter->se)[seqIndex]->used = 1;
		return 1;
	} else {
		cerr << "goUniqueWayRight: no way right exist" << endl;
		return 0;
	}
}

int checkUniqueWayRight(edgesMap &edges, ll finishKmer, Sequence* finishSeq) {
	int count = 0;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl | (finishKmer << (2));
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (finishSeq->similar(*((iter->se)[i]->lower), l - 1, 1)) {
					count++;
					if (count > 1) {
						cerr << "count > 1" << endl;
						return 0;
					}
				}
			}
		}
	}
	if (count == 1) {
		return 1;
	} else {
		cerr << "checkUniqueWayLeft count = 0" << endl;
		return 0;
	}
}
int storeVertex(gvis::GraphScheme<int> &g, vertecesMap &verts, ll newKmer,
		Sequence* newSeq) {
	vertecesMap::iterator iter = verts.find(newKmer);
	if (iter != verts.end()) {
		int size = iter->second.size();
		forn(i, size) {
			if (newSeq->similar(*((iter->se)[i]->lower), l - 1, 0))
				return (iter->se)[i]->start;
		}
		VertexPrototype *v = new VertexPrototype();
		v->lower = newSeq;
		v->start = VertexCount;
		v->finish = 0;
		v->used = 0;
		VertexCount++;

		(iter->se).pb(v);
		g.addVertex(VertexCount - 1, decompress(newKmer, k - 1));
		return VertexCount - 1;
	} else {
		vector<VertexPrototype *> prototypes;
		VertexPrototype *v = new VertexPrototype();
		v->lower = newSeq;
		v->start = VertexCount;
		v->finish = 0;
		v->used = 0;
		VertexCount++;
		prototypes.pb(v);
		verts.insert(mp(newKmer, prototypes));
		g.addVertex(VertexCount - 1, decompress(newKmer, k - 1));
		return VertexCount - 1;
	}

}

