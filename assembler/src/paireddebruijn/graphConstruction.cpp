#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"
#include "seq.hpp"
#include "graphVisualizer.hpp"

int VertexCount;
char EdgeStr[1000000];
char EdgeStrLo[1000000];
int inD[MAX_VERT_NUMBER], outD[MAX_VERT_NUMBER];
int outputEdges[MAX_VERT_NUMBER][MAX_DEGREE];
int inputEdges[MAX_VERT_NUMBER][MAX_DEGREE];

const int minIntersect = l - 1;
int EdgeId;

using namespace paired_assembler;

void constructGraph() {
	//	readsToPairs(parsed_reads, parsed_k_l_mers);
	//	pairsToSequences(parsed_k_l_mers, parsed_k_sequence);
	cerr << "Read edges" << endl;
	edgesMap edges = sequencesToMap(parsed_k_sequence, true);
	cerr << "go to graph" << endl;
	gvis::GraphScheme<int> g("Paired");
	cerr << "Start vertices" << endl;
	VertexCount = 0;
	verticesMap verts;
	longEdgesMap longEdges;
	createVertices(g, edges, verts, longEdges);
	expandDefinite(verts, longEdges);
	outputLongEdges(longEdges);
}

edgesMap sequencesToMap(string parsed_k_sequence, bool usePaired) {
	FILE *inFile = fopen(parsed_k_sequence.c_str(), "r");

	vector<VertexPrototype *> prototypes;
	edgesMap res;
	prototypes.reserve(maxSeqLength);
	int count = 0;
	while (1) {
		char s[maxSeqLength];
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
			VertexPrototype *v = new VertexPrototype(seq, 0);
			if (usePaired || !i)
				prototypes.pb(v);
		}
		res.insert(mp(kmer, prototypes));
	}
	return res;
}

void createVertices(gvis::GraphScheme<int> &g, edgesMap &edges, verticesMap &verts, longEdgesMap &longEdges) {
	char Buffer[2000];
	EdgeId = 0;
	cerr << "Start createVertices " << edges.size() << endl;
	forn(i,MAX_VERT_NUMBER) {
		inD[i] = 0;
		outD[i] = 0;
	}

	for (edgesMap::iterator iter = edges.begin(); iter != edges.end();) {
		int size = iter->second.size();
		ll kmer = iter->fi;
		forn(i, size) {
			if ((!(iter->se)[i]->used)) {
				int length = 1;
				sprintf(EdgeStr + 500000, "%s", decompress(kmer, k).c_str());
				sprintf(EdgeStrLo + 500000, "%s", ((iter->se)[i])->lower->str().c_str());
				(iter->se)[i]->used = 1;
				ll finishKmer = kmer & (~((ll) 3 << (2 * (k - 1))));
				Sequence *finishSeq = new Sequence(
						(iter->se)[i]->lower->Subseq(1,
								(iter->se)[i]->lower->size()));
				ll startKmer = kmer >> 2;
				Sequence *startSeq = new Sequence(
						(iter->se)[i]->lower->Subseq(0,
								(iter->se)[i]->lower->size() - 1));
				length += expandRight(edges, verts, finishKmer, finishSeq);
				int toVert = storeVertex(g, verts, finishKmer, finishSeq);

				int toleft = expandLeft(edges, verts, startKmer, startSeq);
				length += toleft;
				int fromVert = storeVertex(g, verts, startKmer, startSeq);

				cerr << EdgeId << ": (" << length << ") " << ((char*) (EdgeStr
						+ 500000 - toleft)) << endl;
				cerr << ((char*) (EdgeStrLo
						+ 500000 - toleft)) << endl;
				Sequence* UpperSeq = new Sequence(((char*) (EdgeStr
						+ 500000 - toleft)));
				Sequence* LowerSeq = new Sequence(((char*) (EdgeStrLo
						+ 500000 - toleft)));
				Edge* newEdge = new Edge(UpperSeq, LowerSeq, fromVert, toVert, length);
				longEdges.insert(make_pair(EdgeId, newEdge));
				if ((length < 300) && (length > k - 1)) {
					EdgeStr[500000 - toleft + length] = 0;
					sprintf(Buffer, "\"%d: (%d) %s\"", EdgeId, length,
							EdgeStr + 500000 - toleft + k - 1);
				} else {
					sprintf(Buffer, "\"%d: (%d)\"", EdgeId, length);
				}

				g.addEdge(fromVert, toVert, Buffer);
				outputEdges[fromVert][outD[fromVert]] = EdgeId;
				inputEdges[toVert][inD[toVert]] = EdgeId;
				outD[fromVert]++;
				inD[toVert]++;

				EdgeId++;

			}
		}
		(iter->second).clear();
		edges.erase(iter++);
	}
	g.output();
	forn(i, VertexCount) {
		cerr << i << " +" << inD[i] << " -" << outD[i] << endl;
	}
}

int expandRight(edgesMap &edges, verticesMap &verts, ll &finishKmer,
		Sequence* &finishSeq) {
	int length = 0;
	verticesMap::iterator iter;
	while (1) {
		iter = verts.find(finishKmer);
		if (iter != verts.end()) {

			int size = iter->second.size();
			forn(i, size) {
				if (finishSeq->similar(*((iter->se)[i]->lower), minIntersect, 0)) {
					return length;
				}
			}
		}

		if (!checkUniqueWayLeft(edges, finishKmer, finishSeq)) {
			return length;
		}

		if (!goUniqueWayRight(edges, finishKmer, finishSeq)) {
			return length;
		} else {
			EdgeStr[500000 + k + length] = nucl(finishKmer & 3);
			EdgeStrLo[500000 + k + length] = nucl((*finishSeq)[finishSeq->size()-1]);
			length++;
			EdgeStr[500000 + k + length] = 0;
			EdgeStrLo[500000 + k + length] = 0;
		}
	}
}

int expandLeft(edgesMap &edges, verticesMap &verts, ll &startKmer,
		Sequence* &startSeq) {
	int length = 0;

	while (1) {
		verticesMap::iterator iter = verts.find(startKmer);
		if (iter != verts.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (startSeq->similar(*((iter->se)[i]->lower), minIntersect, 0)) {

					return length;
				}
			}
		}
		if (!checkUniqueWayRight(edges, startKmer, startSeq)) {
			return length;
		}
		if (!goUniqueWayLeft(edges, startKmer, startSeq)) {
			return length;
		} else {
			length++;
			EdgeStr[500000 - length] = nucl((startKmer >> ((k - 2) * 2)) & 3);
			EdgeStrLo[500000 - length] = nucl((*startSeq)[0]);
		}
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
				if (finishSeq->similar(*(iter->se)[i]->lower, minIntersect, -1)) {
					count++;
					if (count > 1) {
						return 0;
					}
				}
			}
		}
	}
	return count;
}

int goUniqueWayLeft(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq) {
	int count = 0;
	Sequence *PossibleSequence;
	ll PossibleKmer = 0;
	int seqIndex = 0;
	edgesMap::iterator PossibleIter;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl << (2 * (k - 1)) | finishKmer;
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (finishSeq->similar(*((iter->se)[i]->lower), minIntersect, -1)) {
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
	return 0;
}

int goUniqueWayRight(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq) {
	int count = 0;
	Sequence *PossibleSequence;
	ll PossibleKmer = 0;
	int seqIndex = 0;
	edgesMap::iterator PossibleIter;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl | (finishKmer << (2));
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				if (finishSeq->similar(*((iter->se)[i]->lower), minIntersect, 1)) {
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
		(PossibleIter->se)[seqIndex]->used = 1;
		return 1;
	} else {
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
				if (finishSeq->similar(*((iter->se)[i]->lower), minIntersect, 1)) {
					count++;
					if (count > 1) {
						return 0;
					}
				}
			}
		}
	}
	if (count == 1) {
		return 1;
	} else {
		return 0;
	}
}

verticesMap::iterator addKmerToMap(verticesMap &verts, ll kmer) {
	verticesMap::iterator position = verts.find(kmer);
	if (position == verts.end()) {
		vector<VertexPrototype *> prototypes;
		return verts.insert(make_pair(kmer, prototypes)).first;
	} else {
		return position;
	}
}

/*First argument of result is id of the vertex. Second argument is true if new entry was created and false otherwise
 *
 */
pair<int, bool> addVertexToMap(verticesMap &verts, ll newKmer, Sequence* newSeq) {
	verticesMap::iterator position = addKmerToMap(verts, newKmer);
	vector<VertexPrototype *> *sequences = &position->second;
	for (vector<VertexPrototype *>::iterator it = sequences->begin(); it
			!= sequences->end(); ++it) {
		if (newSeq->similar(*((*it)->lower), minIntersect, 0)) {
			return make_pair((*it)->VertexId, false);
		}
	}
	sequences->push_back(new VertexPrototype(newSeq, VertexCount));
	VertexCount++;
	return make_pair(VertexCount - 1, true);
}

int storeVertex(gvis::GraphScheme<int> &g, verticesMap &verts, ll newKmer,
		Sequence* newSeq) {
	pair<int, bool> addResult = addVertexToMap(verts, newKmer, newSeq);
	if (addResult.second) {
		g.addVertex(addResult.first, decompress(newKmer, k - 1));
	}
	return addResult.first;
}

void resetVertexCount() {
	VertexCount = 0;
}

void expandDefinite(verticesMap &verts, longEdgesMap &longEdges){
	longEdgesMap::iterator it;
	int expandEdgeIndex;
	forn(i,VertexCount){
		if (outD[i]==1){
			expandEdgeIndex=outputEdges[i][0];
			int DestVertex = longEdges[expandEdgeIndex]->ToVertex;
			int a=0;
			while ((inputEdges[DestVertex][a]!=expandEdgeIndex))a++;
			assert(a<inD[DestVertex]);
			while (a<inD[DestVertex]-1){
				inputEdges[DestVertex][a]=inputEdges[DestVertex][a+1];
			}
			inD[DestVertex]--;
			forn(j,inD[i]){
				longEdges[inputEdges[i][j]]->ExpandRight(*(longEdges[expandEdgeIndex]));
				inputEdges[DestVertex][inD[DestVertex]] = inputEdges[i][j];
				inD[DestVertex]++;
			}
			it=longEdges.find(expandEdgeIndex);
			longEdges.erase(it);
			outD[i]=0; inD[i]=0;
		}
	}

}


void outputLongEdges(longEdgesMap &longEdges){
	char Buffer[100];
	gvis::GraphScheme<int> g("Paired_ext");
	freopen("data/graph2.dot", "w",stdout);
	for (longEdgesMap::iterator it=longEdges.begin(); it!=longEdges.end();++it){
		sprintf(Buffer,"%i (%i)",it->first, it->second->length);
		g.addEdge(it->second->FromVertex, it->second->ToVertex, Buffer);
	}
	g.output();
}

/*
void traceReads(verticesMap &verts, longEdgesMap &longEdges){

	map<int, vector<int>> InEdges;
	map<int, vector<int>> OutEdges;

	char UpperRead[1000];
	char LowerRead[1000];
	FILE * inFile = fopen(parsed_reads.c_str(),"r");
	while (fscanf(inFile, "%s %s",UpperRead, LowerRead)) {
//		ProcessRead();
//		ProcessVerticess();
	}
}
*/

