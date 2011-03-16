#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"
#include "seq.hpp"
#include "graphVisualizer.hpp"
#include "pairedGraph.hpp"
#define RIGHT 1
#define LEFT -1

int VertexCount;
char EdgeStr[1000000];
char EdgeStrLo[1000000];
//int inD[MAX_VERT_NUMBER], outD[MAX_VERT_NUMBER];
//int outputEdges[MAX_VERT_NUMBER][MAX_DEGREE];
//int inputEdges[MAX_VERT_NUMBER][MAX_DEGREE];
//int fakeOutputVertices[MAX_VERT_NUMBER][MAX_DEGREE];
//int fakeInputVertices[MAX_VERT_NUMBER][MAX_DEGREE];

const int minIntersect = l - 1;
int EdgeId;
using namespace paired_assembler;
PairedGraph graph;

void constructGraph() {

	//		readsToPairs(parsed_reads, parsed_k_l_mers);
	//		pairsToSequences(parsed_k_l_mers, parsed_k_sequence);
	cerr << "Read edges" << endl;
	edgesMap edges = sequencesToMap(parsed_k_sequence, true);
	cerr << "go to graph" << endl;
	gvis::GraphPrinter<int> g("Paired");
	cerr << "Start vertices" << endl;
	VertexCount = 0;
	verticesMap verts;
	longEdgesMap longEdges;
	createVertices(g, edges, verts, longEdges, graph);
	expandDefinite(longEdges , graph);
//	freopen(graph.c_str(), "w",stdout);
	cerr << endl << "End vertices" <<endl;
//	return;

	freopen(graph2.c_str(), "w",stdout);
	outputLongEdges(longEdges);
	cerr << "TraceReads" << endl;

	traceReads(verts, longEdges, graph);
	freopen(threaded_graph.c_str(), "w", stdout);
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

void createVertices(gvis::GraphPrinter<int> &g, edgesMap &edges,
		verticesMap &verts, longEdgesMap &longEdges, PairedGraph &graph) {
	char Buffer[2000];
	EdgeId = 0;
	cerr << "Start createVertices " << edges.size() << endl;
	forn(i,MAX_VERT_NUMBER) {
		graph.inD[i] = 0;
		graph.outD[i] = 0;
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
				cerr << EdgeId << ": (" << length << ") " << ((char*) (EdgeStr
						+ 500000 - toleft)) << endl;
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
				graph.outputEdges[fromVert][graph.outD[fromVert]] = EdgeId;
				graph.inputEdges[toVert][graph.inD[toVert]] = EdgeId;
				graph.outD[fromVert]++;
				graph.inD[toVert]++;

				EdgeId++;

			}
		}
		(iter->second).clear();
		edges.erase(iter++);
	}
	g.output();
	forn(i, VertexCount) {
		cerr << i << " +" << graph.inD[i] << " -" << graph.outD[i] << endl;
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
				cerr << endl << "expanded down_right" <<endl;
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

	if (count <= 1) {
		edgesMap::iterator iter = edges.find(finishKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				Sequence *Ps = (iter->se)[i]->lower;
				if (finishSeq->similar(*Ps, minIntersect, -1)) {
					/*if ((iter->se)[i]->used)
						return 0;
					if (*Ps == *finishSeq) {
						cerr << endl << "sameSeq";
						continue;
					}*/
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

	int tcount = 0;
	if (count <= 1) {
		edgesMap::iterator iter = edges.find(finishKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				Sequence *Ps = (iter->se)[i]->lower;
				if (finishSeq->similar(*Ps, minIntersect, 1)) {
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
	for (vector<VertexPrototype *>::iterator it = v.begin(); it != v.end(); ++it)
		if (finishSeq->similar(*((*it)->lower), minIntersect, direction)) {
			count++;
			if (count > 1) {
				return count;
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

void expandDefinite(longEdgesMap &longEdges, PairedGraph &graph) {
	longEdgesMap::iterator it;
	int expandEdgeIndex;
	cerr << "expandDefiniteStart" << endl;
	forn(i,VertexCount) {
		if ((graph.outD[i] == 1) && (graph.inD[i] > 0)) {
			cerr << i << endl;
			expandEdgeIndex = edgeRealId(graph.outputEdges[i][0], longEdges);
			int DestVertex = longEdges[expandEdgeIndex]->ToVertex;
			int a = 0;
			while ((edgeRealId(graph.inputEdges[DestVertex][a], longEdges)
					!= expandEdgeIndex))
				a++;
			assert(a<graph.inD[DestVertex]);
			while (a < graph.inD[DestVertex] - 1) {
				graph.inputEdges[DestVertex][a] = graph.inputEdges[DestVertex][a + 1];
				a++;
			}
			graph.inD[DestVertex]--;
			forn(j,graph.inD[i]) {
				longEdges[graph.inputEdges[i][j]]->ExpandRight(
						*(longEdges[expandEdgeIndex]));
				graph.inputEdges[DestVertex][graph.inD[DestVertex]] = graph.inputEdges[i][j];
				graph.inD[DestVertex]++;
			}
			it = longEdges.find(expandEdgeIndex);
			longEdges.erase(it);
			graph.outD[i] = 0;
			graph.inD[i] = 0;
		}
	}

	forn(i,VertexCount) {
		if ((graph.inD[i] == 1) && (graph.outD[i] > 0)) {
			cerr << i << endl;
			expandEdgeIndex = edgeRealId(graph.inputEdges[i][0], longEdges);
			int SourceVertex = longEdges[expandEdgeIndex]->FromVertex;
			int a = 0;
			while (edgeRealId(graph.outputEdges[SourceVertex][a], longEdges)
					!= expandEdgeIndex)
				a++;
			assert(a<graph.outD[SourceVertex]);
			while (a < graph.outD[SourceVertex] - 1) {
				graph.outputEdges[SourceVertex][a] = graph.outputEdges[SourceVertex][a + 1];
				a++;
			}
			graph.outD[SourceVertex]--;
			forn(j,graph.outD[i]) {
				longEdges[graph.outputEdges[i][j]]->ExpandLeft(
						*(longEdges[expandEdgeIndex]));
				graph.outputEdges[SourceVertex][graph.outD[SourceVertex]]
						= graph.outputEdges[i][j];
				graph.outD[SourceVertex]++;
			}
			it = longEdges.find(expandEdgeIndex);
			longEdges.erase(it);
			graph.outD[i] = 0;
			graph.inD[i] = 0;
		}
	}
}

void outputLongEdges(longEdgesMap &longEdges) {
	char Buffer[100];
	gvis::GraphPrinter<int> g("Paired_ext");
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
			sprintf(Buffer, "%i (%i)", it->first, it->second->length);
			//		else sprintf(Buffer,"%i (%i) FAKE now it is %d",it->first, it->second->length,it->second->EdgeId);

			g.addEdge(it->second->FromVertex, it->second->ToVertex, Buffer);
			cerr << it->first << " (" << it->second->length << "):" << endl;
			cerr << it->second->upper->str() << endl;
			cerr << it->second->lower->str() << endl;
		}
	}
	g.output();
}

void traceReads(verticesMap &verts, longEdgesMap &longEdges, PairedGraph &graph) {

	map<int, vector<pair<int, int>>> EdgePairs;
	freopen(parsed_reads.c_str(), "r", stdin);
	char *upperNuclRead = new char[readLength + 2];
	char *lowerNuclRead = new char[readLength + 2];
	char *upperRead = new char[readLength + 2];
	char *lowerRead = new char[readLength + 2];
	ll count=0;
	ll upperMask = (((ll) 1) << (2 * (k - 1))) - 1;
	ll lowerMask = (((ll) 1) << (2 * (l - 1))) - 1;
	//	FILE* fout = fopen("data/filtered_reads","w");
	while (nextReadPair(upperNuclRead, lowerNuclRead)) {
		if (!(count & (1024*128 - 1)))
		cerr<<"read number "<<count<<" processed"<<endl;
		count++;
		codeRead(upperNuclRead, upperRead);
		codeRead(lowerNuclRead, lowerRead);
		Sequence* upRead = new Sequence(upperNuclRead);
		Sequence* loRead = new Sequence(lowerNuclRead);
		int shift = (l - k) / 2;
		ll upper = extractMer(upperRead, shift+1, k-1);
		for (int j = 0; j + l < readLength; j++) {
			verticesMap::iterator vertIter = verts.find(upper);
			if (vertIter!=verts.end()) {
				//	cerr<<"kmer found for j="<<j<<endl;
				for (vector<VertexPrototype *>::iterator it = vertIter->second.begin(); it
						!= vertIter->second.end(); ++it) {
					if ((*it)->lower->similar(loRead->Subseq(1+j, l+j),l-1)) {
						//				cerr<<"vertex found for lower "<<(*it)->lower->str()<<endl;
						//						fprintf(fout,"%s %s\n",upperNuclRead,lowerNuclRead);
						int VertId = (*it)->VertexId;
						if ((graph.inD[VertId]!=0)&&(graph.outD[VertId]!=0)) {
							int tmpIn = -1;
							forn(i,graph.inD[VertId]) {
								int CurEdge=graph.inputEdges[VertId][i];
								int subsizeup = longEdges[CurEdge]->upper->size();
								int subsizelo = longEdges[CurEdge]->lower->size();
								if (subsizeup>j+k+shift) subsizeup = j+k+shift;
								if (subsizelo>j+l) subsizelo = j+l;
								//cerr<<"CheckUp "<<longEdges[CurEdge]->upper->str()<<" VS "<<upRead->Subseq(0,j+k+shift).str()<<" with "<<subsizeup<<endl;
								//cerr<<"CheckLo "<<longEdges[CurEdge]->lower->str()<<" VS "<<loRead->Subseq(0,j+l).str()<<" with "<<subsizelo<<endl;
								if ((longEdges[CurEdge]->upper->similar(upRead->Subseq(0,j+k+shift),subsizeup,RIGHT))&&
										(longEdges[CurEdge]->lower->similar(loRead->Subseq(0,j+l),subsizelo,RIGHT))) {
									if (tmpIn ==-1) {
										tmpIn = CurEdge;
									}
									else {
										tmpIn = -1;
										break;
									}
								}
							}
							if (tmpIn == -1) break;
							int tmpOut = -1;

							forn(i,graph.outD[VertId]) {
								int CurEdge=graph.outputEdges[VertId][i];
								int subsizeup = longEdges[CurEdge]->upper->size();
								int subsizelo = longEdges[CurEdge]->lower->size();
								if (subsizeup>readLength-j-1-shift) subsizeup = readLength-j-1-shift;
								if (subsizelo>readLength-j-l) subsizelo = readLength-j-l;
								if ((longEdges[CurEdge]->upper->similar(upRead->Subseq(j+shift+1,readLength),subsizeup,LEFT))&&
										(longEdges[CurEdge]->lower->similar(loRead->Subseq(j+1,readLength),subsizelo,LEFT))) {
									if (tmpOut ==-1) {
										tmpOut = CurEdge;
									}
									else {
										tmpOut = -1;
										break;
									}
								}
							}
							if (tmpOut != -1) {
								map<int, vector<pair<int,int>>>::iterator epIter = EdgePairs.find(VertId);
								if (epIter == EdgePairs.end()) {
									vector<pair<int,int>> pairVect;
									pairVect.pb(make_pair(tmpIn,tmpOut));
									EdgePairs.insert(mp(VertId,pairVect));
									cerr<<"vert "<<VertId<<" "<<tmpIn<<" "<<tmpOut<<endl;
								}
								else
								{
									if (find(epIter->second.begin(),epIter->second.end(),make_pair(tmpIn,tmpOut)) == epIter->second.end()) {
										epIter->second.pb(make_pair(tmpIn,tmpOut));
									}
								}
							}
						}
						break;
					}
				}
			}

			upper <<= 2;
			upper += upperRead[j + l - shift];
			upper &= upperMask;

		}
		delete upRead;
		delete loRead;
	}
	//	fclose(fout);
	forn(curVertId,VertexCount) {
		if ((graph.inD[curVertId]!=0)&&(graph.outD[curVertId]!=0)) {
			cerr<<"Vertex "<<curVertId<<" connect edges:"<< endl;
			forn(i,(EdgePairs[curVertId]).size()) {
				cerr<<(EdgePairs[curVertId])[i].first<<" ("<<longEdges[(EdgePairs[curVertId])[i].first]->FromVertex<<", "<<longEdges[(EdgePairs[curVertId])[i].first]->ToVertex<<")   ";
				cerr<<(EdgePairs[curVertId])[i].second<<" ("<<longEdges[(EdgePairs[curVertId])[i].second]->FromVertex<<", "<<longEdges[(EdgePairs[curVertId])[i].second]->ToVertex<<")"<<endl;
			}
		}
	}

	//resolving simple case InD=OutD and all edge is unique
	forn(curVertId,VertexCount) {
		if ((graph.inD[curVertId]!=0)&&(graph.outD[curVertId]!=0))
		if ((graph.inD[curVertId]==EdgePairs[curVertId].size())&&(graph.outD[curVertId]==EdgePairs[curVertId].size()))
		{
			cerr<<"Process vertex "<<curVertId<<" IN "<<graph.inD[curVertId]<<" OUT "<<graph.outD[curVertId]<<" unique ways "<<EdgePairs[curVertId].size()<<endl;
			bool allPairUnique = true;
			forn(i,(EdgePairs[curVertId]).size()) {
				int CurIn = (EdgePairs[curVertId])[i].first;
				int CurOut = (EdgePairs[curVertId])[i].second;
				forn(j,(EdgePairs[curVertId]).size()) {
					if ((EdgePairs[curVertId])[j].first==CurIn)
					if ((EdgePairs[curVertId])[j].second!=CurOut) allPairUnique = false;
					if ((EdgePairs[curVertId])[j].second==CurOut)
					if ((EdgePairs[curVertId])[j].first!=CurIn) allPairUnique = false;
				}
			}
			forn(i,(EdgePairs[curVertId]).size()) {
				if (allPairUnique) {
					int CurIn = edgeRealId((EdgePairs[curVertId])[i].first, longEdges);
					int CurOut = edgeRealId((EdgePairs[curVertId])[i].second, longEdges);
					cerr<<"New inclusion "<<CurIn<<"("<<longEdges[CurIn]->FromVertex<<","<<longEdges[CurIn]->ToVertex<<") <- "<<CurOut<<"("<<longEdges[CurOut]->FromVertex<<","<<longEdges[CurOut]->ToVertex<<")"<<endl;
					longEdges[CurIn]->ExpandRight(*longEdges[CurOut]);
					longEdges[CurOut] = longEdges[CurIn];
					cerr<<"New edges "<<CurIn<<"("<<longEdges[CurIn]->FromVertex<<","<<longEdges[CurIn]->ToVertex<<") <- "<<CurOut<<"("<<longEdges[CurOut]->FromVertex<<","<<longEdges[CurOut]->ToVertex<<")"<<endl;
					graph.inD[curVertId]--;
					graph.outD[curVertId]--;
				}
			}
		}
	}

	freopen("data/graph_after_obvious.dot", "w", stdout);
	outputLongEdges(longEdges);

	//resolve multi case;
	int FakeVertexCount = VertexCount;
	int FakeVertexStart = VertexCount;
	map<int, int> FakeVertexToReal;
	forn(curVertId,VertexCount) {
		if ((graph.inD[curVertId]!=0)&&(graph.outD[curVertId]!=0)) {
			if ((graph.inD[curVertId]<=EdgePairs[curVertId].size())&&(graph.outD[curVertId]<=EdgePairs[curVertId].size())) {
				bool allIns = false;
				bool allOuts = false;
				forn(i,graph.inD[curVertId]) {
					allIns = false;
					forn (j,EdgePairs[curVertId].size()) {
						if ((EdgePairs[curVertId])[j].first ==graph.inputEdges[curVertId][i]) {
							allIns = true;
							break;
						}
					}
					if (!allIns) break;
				}
				forn(i,graph.outD[curVertId]) {
					allOuts = false;
					forn (j,EdgePairs[curVertId].size()) {
						if ((EdgePairs[curVertId])[j].second ==graph.outputEdges[curVertId][i]) {
							allOuts = true;
							break;
						}
					}
					if (!allOuts) break;
				}
				if (allIns&&allOuts) {
					int tmpCurOut = edgeRealId(graph.outputEdges[curVertId][0],longEdges);
					string tmpLoSeq = longEdges[tmpCurOut]->upper->Subseq(0,k-1).str();
					string tmpUpSeq = longEdges[tmpCurOut]->upper->Subseq(0,l-1).str();

					//create FakeVerices and Fake edges;
					//create fake Vertices for in edges;
					int tmpFictStartIn = FakeVertexCount;
					forn(i,graph.inD[curVertId]) {
						int CurIn = edgeRealId(graph.inputEdges[curVertId][i],longEdges);
						graph.inputEdges[curVertId][i] = CurIn;
						FakeVertexToReal.insert(make_pair(FakeVertexCount,curVertId));
						longEdges[CurIn]->ToVertex = FakeVertexCount;
						graph.inD[FakeVertexCount]=1;
						graph.outD[FakeVertexCount]=0;
						graph.inputEdges[FakeVertexCount][0]=CurIn;
						FakeVertexCount++;
					}
					//create fake Vertices for out edges;
					int tmpFictStartOut = FakeVertexCount;
					forn(i,graph.outD[curVertId]) {
						int CurOut = edgeRealId(graph.outputEdges[curVertId][i],longEdges);
						graph.outputEdges[curVertId][i] = CurOut;
						FakeVertexToReal.insert(make_pair(FakeVertexCount,curVertId));
						longEdges[CurOut]->FromVertex = FakeVertexCount;
						graph.inD[FakeVertexCount]=0;
						graph.outD[FakeVertexCount]=1;
						graph.outputEdges[FakeVertexCount][0]=CurOut;
						FakeVertexCount++;
					}
					//create fake edges

					forn (tmpEdgePair,EdgePairs[curVertId].size()) {
						int tmpFrom = 0;
						int tmpTo = 0;
						while (edgeRealId(graph.inputEdges[curVertId][tmpFrom], longEdges)!=edgeRealId((EdgePairs[curVertId])[tmpEdgePair].first,longEdges)) {
							tmpFrom++;
							assert(tmpFrom<graph.inD[curVertId]);
						}
						while (edgeRealId(graph.outputEdges[curVertId][tmpTo],longEdges)!=edgeRealId((EdgePairs[curVertId])[tmpEdgePair].second, longEdges)) {
							tmpTo++;
							assert(tmpTo<graph.outD[curVertId]);
						}
						Sequence *UpSeq = new Sequence(tmpLoSeq);
						Sequence *LoSeq = new Sequence(tmpUpSeq);

						Edge *tmpEdge = new Edge(UpSeq, LoSeq, tmpFictStartIn + tmpFrom, tmpFictStartOut + tmpTo,0, EdgeId);
						longEdges.insert(make_pair(EdgeId,tmpEdge));
						graph.outputEdges[tmpFictStartIn + tmpFrom][graph.outD[tmpFictStartIn + tmpFrom]] = EdgeId;
						graph.outD[tmpFictStartIn + tmpFrom]++;
						graph.inputEdges[tmpFictStartOut + tmpTo][graph.inD[tmpFictStartOut + tmpTo]] = EdgeId;
						graph.inD[tmpFictStartOut + tmpTo]++;
						EdgeId++;
					}

					graph.inD[curVertId]=0;
					graph.outD[curVertId]=0;
				}

			}
		}
	}

	VertexCount = FakeVertexCount;
	expandDefinite(longEdges, graph);
	for(longEdgesMap::iterator it= longEdges.begin(); it !=longEdges.end(); ++it) {
		if (it->second->FromVertex>=FakeVertexStart) it->second->FromVertex = FakeVertexToReal[it->second->FromVertex];
		if (it->second->ToVertex>=FakeVertexStart) it->second->ToVertex = FakeVertexToReal[it->second->ToVertex];
	}

	freopen("data/graph_after_fake.dot", "w", stdout);
	outputLongEdges(longEdges);

}

