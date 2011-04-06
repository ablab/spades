
#include "constdConstruction.hpp"

LOGGER("p.constdConstruction");


namespace constd {

char EdgeStr[1000000];
char EdgeStrLo[1000000];
void createVertices(edgesMap &edges, PairedGraph &graph) {
	char Buffer[2000];
	graph.EdgeId = 0;
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
				int EdgeCoverage;
				count++;
//				cerr << count << endl;
				if (((iter->se)[i])->lower->size() < l) {
					DEBUG("Bad edge: "<<((iter->se)[i])->lower->size()<<" "<<((iter->se)[i])->lower->str());
				}
				assert (((iter->se)[i])->lower->size() >= l);
				EdgeCoverage = (iter->se)[i]->coverage;
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
				length += expandRight(edges, graph.verts, finishKmer, finishSeq, EdgeCoverage);
				int toVert = storeVertex(graph, finishKmer, finishSeq);
				int toleft = expandLeft(edges, graph.verts, startKmer, startSeq, EdgeCoverage);
				//cerr<<endl << "TO LEFT" << toleft << endl;
				if (graph.verts.size() > 20000) {
					cerr << endl <<" Too much vertices";
					assert(0);
				}
				length += toleft;
				int fromVert = storeVertex(graph, startKmer, startSeq);
				//				if (! (count && ((1<<14) -1 ))) {
//				cerr << graph.EdgeId << ": (" << length << ") " << ((char*) (EdgeStr
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
						length, graph.EdgeId, EdgeCoverage);
				graph.longEdges.insert(make_pair(graph.EdgeId, newEdge));
				if ((length < 300) && (length > k - 1)) {
					EdgeStr[500000 - toleft + length] = 0;
					sprintf(Buffer, "%d: (%d) %s", graph.EdgeId, length,
							EdgeStr + 500000 - toleft + k - 1);
				} else {
					sprintf(Buffer, "%d: (%d)", graph.EdgeId, length);
				}

				graph.edgeIds[fromVert][graph.degrees[fromVert][1]][OUT_EDGE] = graph.EdgeId;
				graph.edgeIds[toVert][graph.degrees[toVert][0]][IN_EDGE] = graph.EdgeId;
				graph.degrees[fromVert][1]++;
				graph.degrees[toVert][0]++;

				graph.EdgeId++;

			}
		}
		(iter->second).clear();
		edges.erase(iter++);
	}
	forn(i, graph.VertexCount) {
		cerr << i << " +" << graph.degrees[i][0] << " -" << graph.degrees[i][1] << endl;
	}
}



int expandRight(edgesMap &edges, verticesMap &verts, ll &finishKmer,
		Sequence* &finishSeq, int &EdgeCoverage) {
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
		if ((go_res = goUniqueWayRight(edges, finishKmer, finishSeq, EdgeCoverage)) == 0) {
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
		Sequence* &startSeq, int &EdgeCoverage) {
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
		if ( 0 == (go_res = goUniqueWayLeft(edges, startKmer, startSeq, EdgeCoverage))) {
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

int goUniqueWayLeft(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage) {
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

	if (count == 1 || tcount == 1) {

		if (EdgeCoverage < (PossibleIter->se)[seqIndex]->coverage)
			EdgeCoverage = (PossibleIter->se)[seqIndex]->coverage;
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

int goUniqueWayRight(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage) {
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

	bool sameK = false;

	if (count == 1 ) {
		if (EdgeCoverage < (PossibleIter->se)[seqIndex]->coverage)
			EdgeCoverage = (PossibleIter->se)[seqIndex]->coverage;
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

int countWays(vector<EdgePrototype *> &v, Sequence *finishSeq, int direction) {
	int count = 0;
//	cerr <<" countWays started"<< endl;
	for (vector<EdgePrototype *>::iterator it = v.begin(); it != v.end(); ++it) {
		if (finishSeq->similar(*((*it)->lower), minIntersect, direction)) {
			count++;
			if (count > 1) {
				return count;
			}
		}
	}
	return count;
}


int checkUniqueWay(edgesMap &edges, ll finishKmer, Sequence *finishSeq,
		int direction) {
	int count = 0;
//	cerr << "findUniqueWay" << endl;
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

}
