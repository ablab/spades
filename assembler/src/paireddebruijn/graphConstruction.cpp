#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"
#include "seq.hpp"
#include "graphVisualizer.hpp"
#include "pairedGraph.hpp"
#include "graphio.hpp"
#include "readTracing.hpp"
#include "graphSimplification.hpp"

LOGGER("p.graphConstruction");
int VertexCount;
char EdgeStr[1000000];
char EdgeStrLo[1000000];
//int inD[MAX_VERT_NUMBER], outD[MAX_VERT_NUMBER];
//int outputEdges[MAX_VERT_NUMBER][MAX_DEGREE];
//int inputEdges[MAX_VERT_NUMBER][MAX_DEGREE];
//int fakeOutputVertices[MAX_VERT_NUMBER][MAX_DEGREE];
//int fakeInputVertices[MAX_VERT_NUMBER][MAX_DEGREE];

int minIntersect ;
int EdgeId;
using namespace paired_assembler;
PairedGraph graph;
longEdgesMap longEdges;

void constructGraph() {

	minIntersect = l - 1;
	INFO ("Read edges");
	edgesMap edges = sequencesToMap(parsed_k_sequence, true);
	INFO ("go to graph");
	gvis::GraphPrinter<int> g("Paired");
	INFO ("Start vertices");
	VertexCount = 0;
	verticesMap verts;
	createVertices(g, edges, verts, longEdges, graph);
//	vertexDist(longEdges, graph, 172);
	freopen("data/beforeExpand.dot", "w",stdout);
	outputLongEdges(longEdges, graph);

	expandDefinite(longEdges , graph, VertexCount, true);

	freopen("data/afterExpand.dot", "w",stdout);
	outputLongEdges(longEdges, graph);
	freopen("data/afterExpand_g.dot", "w",stdout);
	outputLongEdgesThroughGenome(longEdges,graph,VertexCount);
//	freopen(graph.c_str(), "w",stdout);
	cerr << endl << "End vertices" <<endl;
//	return;

//	freopen(graph2.c_str(), "w",stdout);
//	outputLongEdges(longEdges);
//	cerr << "TraceReads" << endl;

	traceReads(verts, longEdges, graph, VertexCount, EdgeId);
	freopen("data/ReadsTraced.dot", "w", stdout);
	outputLongEdges(longEdges);
	freopen("data/ReadsTraced_g.dot", "w", stdout);
	outputLongEdgesThroughGenome(longEdges,graph,VertexCount);
	graph.recreateVerticesInfo(VertexCount, longEdges);

	while (processLowerSequence(longEdges, graph, VertexCount))
	{
		graph.recreateVerticesInfo(VertexCount, longEdges);
		expandDefinite(longEdges , graph, VertexCount);
	}
	freopen("data/afterLowers.dot", "w",stdout);
	outputLongEdges(longEdges);
	freopen("data/afterLowers_g.dot", "w",stdout);
	outputLongEdgesThroughGenome(longEdges,graph,VertexCount);

	graph.recreateVerticesInfo(VertexCount, longEdges);
	freopen("data/afterLowers_info.dot", "w",stdout);
	outputLongEdges(longEdges, graph);

	//freopen("data/LowerProcessed.dot", "w", stdout);
///	outputLongEdges(longEdges);

	freopen("data/afterExtraxtDefinite1.dot", "w",stdout);
	extractDefinite(longEdges , graph, VertexCount,1);

	outputLongEdges(longEdges, graph);
	freopen("data/afterExtraxtDefinite2.dot", "w",stdout);
	extractDefinite(longEdges , graph, VertexCount,0);

	outputLongEdges(longEdges, graph);

	cerr << "\n Finished";


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

void createVertices(gvis::GraphPrinter<int> &g, edgesMap &edges,
		verticesMap &verts, longEdgesMap &longEdges, PairedGraph &graph) {
	char Buffer[2000];
	EdgeId = 0;
	cerr << "Start createVertices " << edges.size() << endl;
	forn(i, MAX_VERT_NUMBER) {
		graph.degrees[i][0] = 0;
		graph.degrees[i][1] = 0;

	}
	int count = 0;
	for (edgesMap::iterator iter = edges.begin(); iter != edges.end();) {
		int size = iter->second.size();
		ll kmer = iter->fi;
		forn(i, size) {
			if ((!(iter->se)[i]->used)) {
				int length = 1;
				count++;
//				cerr << count << endl;
				if (((iter->se)[i])->lower->size() < l) {
					cerr<<"Bad edge: "<<((iter->se)[i])->lower->size()<<" "<<((iter->se)[i])->lower->str()<<endl;
				}
				assert (((iter->se)[i])->lower->size() >= l);
				sprintf(EdgeStr + 500000, "%s", decompress(kmer, k).c_str());
				sprintf(EdgeStrLo + 500000, "%s",
						((iter->se)[i])->lower->str().c_str());
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
				//cerr<<endl << "TO LEFT" << toleft << endl;
				if (verts.size() > 20000) {
					cerr << endl <<" Too much vertices";
					assert(0);
				}
				length += toleft;
				int fromVert = storeVertex(g, verts, startKmer, startSeq);
				//				if (! (count && ((1<<14) -1 ))) {
//				cerr << EdgeId << ": (" << length << ") " << ((char*) (EdgeStr
//						+ 500000 - toleft)) << endl;
//				}
				Sequence* UpperSeq = new Sequence(((char*) (EdgeStr
						+ 500000 - toleft)));
				Sequence* LowerSeq = new Sequence(((char*) (EdgeStrLo
						+ 500000 - toleft)));
				if (LowerSeq->size() < l) {
					cerr << endl << LowerSeq->str() << endl << UpperSeq->str()
							<< endl;
					assert(0);
				}
				Edge* newEdge = new Edge(UpperSeq, LowerSeq, fromVert, toVert,
						length, EdgeId);
				longEdges.insert(make_pair(EdgeId, newEdge));
				if ((length < 300) && (length > k - 1)) {
					EdgeStr[500000 - toleft + length] = 0;
					sprintf(Buffer, "%d: (%d) %s", EdgeId, length,
							EdgeStr + 500000 - toleft + k - 1);
				} else {
					sprintf(Buffer, "%d: (%d)", EdgeId, length);
				}

				g.addEdge(fromVert, toVert, Buffer);
				graph.edgeIds[fromVert][graph.degrees[fromVert][1]][OUT_EDGE] = EdgeId;
				graph.edgeIds[toVert][graph.degrees[toVert][0]][IN_EDGE] = EdgeId;
				graph.degrees[fromVert][1]++;
				graph.degrees[toVert][0]++;

				EdgeId++;

			}
		}
		(iter->second).clear();
		edges.erase(iter++);
	}
	g.output();
	forn(i, VertexCount) {
		cerr << i << " +" << graph.degrees[i][0] << " -" << graph.degrees[i][1] << endl;
	}
}

int expandRight(edgesMap &edges, verticesMap &verts, ll &finishKmer,
		Sequence* &finishSeq) {
	int length = 0;
	verticesMap::iterator iter;
//	cerr << "expand_right"<<endl;
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
//		cerr << "before checkUniqueWayLeft" << "<<" << finishKmer << " " << finishSeq->str();
		if (!checkUniqueWayLeft(edges, finishKmer, finishSeq)) {
			return length;
		}
//		cerr << "after checkUniqueWayLeft";
		int go_res = 0;
		if ((go_res = goUniqueWayRight(edges, finishKmer, finishSeq)) == 0) {
			return length;
		} else {
//			cerr << "expanding right" << endl;
			if (go_res == 1) {
				EdgeStr[500000 + k + length] = nucl(finishKmer & 3);
				EdgeStrLo[500000 + l + length] = nucl((*finishSeq)[finishSeq->size()-1]);
				length++;
				EdgeStr[500000 + k + length] = 0;
				EdgeStrLo[500000 + l + length] = 0;
			} else {
//				cerr << endl << "expanded down_right" <<endl;
			}
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
		int go_res;
//		cerr << endl<<"trying to go left"<<endl;
		if ( 0 == (go_res = goUniqueWayLeft(edges, startKmer, startSeq))) {
			return length;
		} else {
			if (go_res != 2) {
				length++;
				EdgeStr[500000 - length] = nucl((startKmer >> ((k - 2) * 2)) & 3);
				EdgeStrLo[500000 - length] = nucl((*startSeq)[0]);
			}
			else {
				cerr << endl<<startKmer <<" expanded down_left"<<endl;
			}
		}
	}
}

int goUniqueWayLeft(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq) {
	int count = 0;
	Sequence *PossibleSequence;
	ll PossibleKmer = 0;
	int seqIndex = 0;
	edgesMap::iterator PossibleIter;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (((ll) Nucl) << (2 * (k - 1))) | finishKmer;
		edgesMap::iterator iter = edges.find(tmpKmer);
		//		cerr << endl << tmpKmer;
		//		assert(0);
		if (iter != edges.end()) {

			//			cerr<< endl <<"found left kmer" <<endl;
			int size = iter->second.size();
			forn(i, size) {
				if (finishSeq->similar(*((iter->se)[i]->lower), minIntersect,
						-1)) {
					count++;
					if (count > 1)
						return 0;

					//					cerr<< endl <<"found something" <<endl;
					PossibleKmer = tmpKmer;
					PossibleSequence = (iter->se)[i]->lower;
					seqIndex = i;
					PossibleIter = iter;
				}
			}
		}
	}
	bool sameK = false;
	int tcount = 0;
/*
	if (count <= 1) {
		edgesMap::iterator iter = edges.find(finishKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				Sequence *Ps = (iter->se)[i]->lower;
				if (finishSeq->similar(*Ps, minIntersect, -1)) {
					tcount++;
					if (tcount > 1)
						break;
					PossibleKmer = finishKmer;
					PossibleSequence = (iter->se)[i]->lower;
					seqIndex = i;
					PossibleIter = iter;
					sameK = true;
				}
			}
		}
		sameK = true;
		if (count == 1 && sameK) {
			//			assert("1 == 0");
			//			cerr << endl << "something" << endl;
		}
	}
*/
	if (count == 1 || tcount == 1) {
		finishKmer = PossibleKmer >> 2;
		finishSeq = new Sequence(
				(PossibleIter->se)[seqIndex]->lower->Subseq(0,
						(PossibleIter->se)[seqIndex]->lower->size() - 1));//PossibleSequence;
		(PossibleIter->se)[seqIndex]->used = 1;
		if (sameK)
			return 2;
		else
			return 1;
	}
	cerr << " suspend";
	return 0;
}

int goUniqueWayRight(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq) {
	int count = 0;
	Sequence *PossibleSequence;
	ll PossibleKmer = 0;
	int seqIndex = 0;

//	cerr << " goUniqueWayRight" << endl;
	edgesMap::iterator PossibleIter;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = (ll) Nucl | (finishKmer << (2));
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				Sequence *Ps = (iter->se)[i]->lower;
				if (finishSeq->similar(*Ps, minIntersect, 1)) {
					if ((iter->se)[i]->used)
						break;
					count++;
					if (count > 1)
						return 0;
					PossibleKmer = tmpKmer;
					PossibleSequence = Ps;
					seqIndex = i;
					PossibleIter = iter;

				}
			}
		}
	}
	int tcount = 0;

	bool sameK = false;

	/*
	 	if (count <= 1) {
		edgesMap::iterator iter = edges.find(finishKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				Sequence *Ps = (iter->se)[i]->lower;
				if (finishSeq->similar(*Ps, minIntersect, 1))
				{
					if ((iter->se)[i]->used)
						break;
					tcount++;
					if (tcount > 1)
						break;
					PossibleKmer = finishKmer;
					PossibleSequence = (iter->se)[i]->lower;
					seqIndex = i;
					PossibleIter = iter;
					sameK = true;
				}
			}
		}
		//		assert(1 == 0);
		sameK = true;
	}
*/
	if (count == 1 || tcount == 1) {
		int tcount = 0;
		finishKmer = (PossibleKmer) & (~(((ll) 3) << (2 * (k - 1))));
		finishSeq = new Sequence(
				(PossibleIter->se)[seqIndex]->lower->Subseq(1,
						(PossibleIter->se)[seqIndex]->lower->size()));//PossibleSequence;
		(PossibleIter->se)[seqIndex]->used = 1;
		if (sameK)
			return 2;
		else
			return 1;
	} else {
		return 0;
	}
}

int countWays(vector<VertexPrototype *> &v, Sequence *finishSeq, int direction) {
	int count = 0;
//	cerr <<" countWays started"<< endl;
	for (vector<VertexPrototype *>::iterator it = v.begin(); it != v.end(); ++it) {
		if (finishSeq->similar(*((*it)->lower), minIntersect, direction)) {
			count++;
			if (count > 1) {
				return count;
			}
		}
	}
	return count;
}

/*
 * Method adds nucleotide to the side of kMer defined by direction
 */ll pushNucleotide(ll kMer, int length, int direction, int nucl) {
	if (direction == RIGHT) {
		return (ll) nucl | (kMer << (2));
	} else {
		return (ll) nucl << (2 * length) | kMer;
	}
}

int checkUniqueWay(edgesMap &edges, ll finishKmer, Sequence *finishSeq,
		int direction) {
	int count = 0;
//	cerr << "checkUniqueWay" << endl;
	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = pushNucleotide(finishKmer, k - 1, direction, Nucl);
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			count += countWays(iter->second, finishSeq, direction);
			if (count > 1)
				return 0;
		}
	}
	return count == 1;
}

int checkUniqueWayLeft(edgesMap &edges, ll finishKmer, Sequence *finishSeq) {
	return checkUniqueWay(edges, finishKmer, finishSeq, LEFT);
}

int checkUniqueWayRight(edgesMap &edges, ll finishKmer, Sequence* finishSeq) {
	return checkUniqueWay(edges, finishKmer, finishSeq, RIGHT);
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

/*First argument of result is id of the vertex. Second argument is true if new entry
 * was created and false otherwise
 *
 */
pair<int, bool> addVertexToMap(verticesMap &verts, ll newKmer, Sequence* newSeq) {
	verticesMap::iterator position = addKmerToMap(verts, newKmer);
	vector<VertexPrototype *> *sequences = &position->second;
	for (vector<VertexPrototype *>::iterator it = sequences->begin(); it
			!= sequences->end(); ++it) {
		Sequence *otherSequence = (*it)->lower;
		if (newSeq->similar(*otherSequence, minIntersect, 0)) {
			return make_pair((*it)->VertexId, false);
		}
	}
	sequences->push_back(new VertexPrototype(newSeq, VertexCount));
	VertexCount++;
	return make_pair(VertexCount - 1, true);
}

int storeVertex(gvis::GraphPrinter<int> &g, verticesMap &verts, ll newKmer,
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


