#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"
#include "seq.hpp"
#include "graphVisualizer.hpp"
#define RIGHT 1
#define LEFT -1

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
	createVertices(g, edges, verts, longEdges);
	cerr << endl << "End vertices" <<endl;
//	return;
	expandDefinite(verts, longEdges);
	freopen(graph.c_str(), "w",stdout);
	outputLongEdges(longEdges);
	cerr << "TraceReads" << endl;

	traceReads(verts, longEdges);
	freopen(threaded_graph.c_str(), "w",stdout);
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

void createVertices(gvis::GraphPrinter<int> &g, edgesMap &edges, verticesMap &verts, longEdgesMap &longEdges) {
	char Buffer[2000];
	EdgeId = 0;
	cerr << "Start createVertices " << edges.size() << endl;
	forn(i,MAX_VERT_NUMBER) {
		inD[i] = 0;
		outD[i] = 0;
	}
	int count = 0;
	for (edgesMap::iterator iter = edges.begin(); iter != edges.end();) {
		int size = iter->second.size();
		ll kmer = iter->fi;
		forn(i, size) {
			if ((!(iter->se)[i]->used)) {
				int length = 1;
				count ++;
//				cerr << count << endl;
				assert (((iter->se)[i])->lower->size() >= l);
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
				//cerr<<endl << "TO LEFT" << toleft << endl;
				length += toleft;
				int fromVert = storeVertex(g, verts, startKmer, startSeq);
//				if (! (count && ((1<<14) -1 ))) {
					cerr << EdgeId << ": (" << length << ") " << ((char*) (EdgeStr
						+ 500000 - toleft)) << endl;
					cerr << ((char*) (EdgeStrLo
						+ 500000 - toleft)) << endl;
	//			}
				Sequence* UpperSeq = new Sequence(((char*) (EdgeStr
						+ 500000 - toleft)));
				Sequence* LowerSeq = new Sequence(((char*) (EdgeStrLo
						+ 500000 - toleft)));
				if (LowerSeq->size() < l) {
					cerr<<endl<< LowerSeq->str()<<endl << UpperSeq->str()<<endl;
					assert(0);
				}
				Edge* newEdge = new Edge(UpperSeq, LowerSeq, fromVert, toVert, length, EdgeId);
				longEdges.insert(make_pair(EdgeId, newEdge));
				if ((length < 300) && (length > k - 1)) {
					EdgeStr[500000 - toleft + length] = 0;
					sprintf(Buffer, "%d: (%d) %s", EdgeId, length,
							EdgeStr + 500000 - toleft + k - 1);
				} else {
					sprintf(Buffer, "%d: (%d)", EdgeId, length);
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

		if (goUniqueWayRight(edges, finishKmer, finishSeq) != 1) {
			return length;
		} else {
//			cerr << "expanding right" << endl;
			EdgeStr[500000 + k + length] = nucl(finishKmer & 3);
			EdgeStrLo[500000 + l + length] = nucl((*finishSeq)[finishSeq->size()-1]);
			length++;
			EdgeStr[500000 + k + length] = 0;
			EdgeStrLo[500000 + l + length] = 0;
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
		if (1 != (go_res = goUniqueWayLeft(edges, startKmer, startSeq))) {
			return length;
		} else {
//			if (go_res != 2) {
				length++;
				EdgeStr[500000 - length] = nucl((startKmer >> ((k - 2) * 2)) & 3);
				EdgeStrLo[500000 - length] = nucl((*startSeq)[0]);
//			}
//			else {
//				cerr << endl<<"kmer_expanding left"<<endl;
//			}
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
				if (finishSeq->similar(*((iter->se)[i]->lower), minIntersect, -1)) {
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
	if (count == -1) {
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
					count++;
					if (count > 1)
						return 0;
					PossibleKmer = finishKmer;
					PossibleSequence = (iter->se)[i]->lower;
					seqIndex = i;
					PossibleIter = iter;

				}
			}
		}
//		assert(1 == 0);
		sameK = true;
		if (count == 1 && sameK) {
//			assert("1 == 0");
//			cerr << endl << "something" << endl;
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
				Sequence *Ps = (iter->se)[i]->lower;
				if (finishSeq->similar(*Ps, minIntersect, 1)) {
					if ((iter->se)[i]->used)
						return 0;
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
	if (count == -1) {
		edgesMap::iterator iter = edges.find(finishKmer);
		if (iter != edges.end()) {
			int size = iter->second.size();
			forn(i, size) {
				Sequence *Ps = (iter->se)[i]->lower;
				if (finishSeq->similar(*Ps, minIntersect, 1)) {
					if ((iter->se)[i]->used)
						return 0;
					count++;
					if (count > 1)
						return 0;
					PossibleKmer = finishKmer;
					PossibleSequence = (iter->se)[i]->lower;
					seqIndex = i;
					PossibleIter = iter;

				}
			}
		}
//		assert(1 == 0);
		sameK = true;
	}
	if (count == 1) {
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
 */
ll pushNucleotide(ll kMer, int length, int direction, int nucl) {
	if(direction == RIGHT) {
		return (ll) nucl | (kMer << (2));
	} else {
		return (ll) nucl << (2 * length) | kMer;
	}
}

int checkUniqueWay(edgesMap &edges, ll finishKmer, Sequence *finishSeq, int direction) {
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

void expandDefinite(verticesMap &verts, longEdgesMap &longEdges){
	longEdgesMap::iterator it;
	int expandEdgeIndex;
	cerr << "expandDefiniteStart"<<endl;
	forn(i,VertexCount){
		if ((outD[i]==1)&&(inD[i]>0)){
			cerr << i << endl;
			expandEdgeIndex=outputEdges[i][0];
			int DestVertex = longEdges[expandEdgeIndex]->ToVertex;
			int a=0;
			while ((inputEdges[DestVertex][a]!=expandEdgeIndex))a++;
			assert(a<inD[DestVertex]);
			while (a<inD[DestVertex]-1){
				inputEdges[DestVertex][a]=inputEdges[DestVertex][a+1];
				a++;
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


	forn(i,VertexCount){
		if ((inD[i]==1)&&(outD[i]>0)){
			cerr << i << endl;
			expandEdgeIndex=inputEdges[i][0];
			int SourceVertex = longEdges[expandEdgeIndex]->FromVertex;
			int a=0;
			while ((outputEdges[SourceVertex][a]!=expandEdgeIndex))a++;
			assert(a<outD[SourceVertex]);
			while (a<outD[SourceVertex]-1){
				outputEdges[SourceVertex][a]=outputEdges[SourceVertex][a+1];
				a++;
			}
			outD[SourceVertex]--;
			forn(j,outD[i]){
				longEdges[outputEdges[i][j]]->ExpandLeft(*(longEdges[expandEdgeIndex]));
				outputEdges[SourceVertex][outD[SourceVertex]] = outputEdges[i][j];
				outD[SourceVertex]++;
			}
			it=longEdges.find(expandEdgeIndex);
			longEdges.erase(it);
			outD[i]=0; inD[i]=0;
		}
	}
}


void outputLongEdges(longEdgesMap &longEdges){
	char Buffer[100];
	gvis::GraphPrinter<int> g("Paired_ext");
	for (longEdgesMap::iterator it=longEdges.begin(); it!=longEdges.end();++it){
		if (it->second->EdgeId == it->first)
			{
			sprintf(Buffer,"%i (%i)",it->first, it->second->length);
	//		else sprintf(Buffer,"%i (%i) FAKE now it is %d",it->first, it->second->length,it->second->EdgeId);

			g.addEdge(it->second->FromVertex, it->second->ToVertex, Buffer);
			cerr<<it->first<<" ("<<it->second->length<<"):"<<endl;
			cerr<<it->second->upper->str()<<endl;
			cerr<<it->second->lower->str()<<endl;
		}
	}
	g.output();
}

void traceReads(verticesMap &verts, longEdgesMap &longEdges){

	map<int, vector<pair<int,int>>> EdgePairs;
	freopen(parsed_reads.c_str(), "r", stdin);
	char *upperNuclRead = new char[readLength + 2];
	char *lowerNuclRead = new char[readLength + 2];
	char *upperRead = new char[readLength + 2];
	char *lowerRead = new char[readLength + 2];
	ll count=0;
	ll upperMask = (((ll) 1) << (2 * (k - 1))) - 1;
	ll lowerMask = (((ll) 1) << (2 * (l - 1))) - 1;
	FILE* fout = fopen("data/filtered_reads","w");
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
					if ((*it)->lower->similar(loRead->Subseq(1+j, l+j),l-1)){
		//				cerr<<"vertex found for lower "<<(*it)->lower->str()<<endl;
						fprintf(fout,"%s %s\n",upperNuclRead,lowerNuclRead);
						int VertId = (*it)->VertexId;
						if ((inD[VertId]!=0)&&(outD[VertId]!=0)) {
							int tmpIn = -1;
							forn(i,inD[VertId]){
								int CurEdge=inputEdges[VertId][i];
								int subsizeup = longEdges[CurEdge]->upper->size();
								int subsizelo = longEdges[CurEdge]->lower->size();
								if (subsizeup>j+k+shift) subsizeup = j+k+shift;
								if (subsizelo>j+l) subsizelo = j+l;
								//cerr<<"CheckUp "<<longEdges[CurEdge]->upper->str()<<" VS "<<upRead->Subseq(0,j+k+shift).str()<<" with "<<subsizeup<<endl;
								//cerr<<"CheckLo "<<longEdges[CurEdge]->lower->str()<<" VS "<<loRead->Subseq(0,j+l).str()<<" with "<<subsizelo<<endl;
								if ((longEdges[CurEdge]->upper->similar(upRead->Subseq(0,j+k+shift),subsizeup,RIGHT))&&
									(longEdges[CurEdge]->lower->similar(loRead->Subseq(0,j+l),subsizelo,RIGHT))){
									if (tmpIn ==-1){
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

							forn(i,outD[VertId]){
								int CurEdge=outputEdges[VertId][i];
								int subsizeup = longEdges[CurEdge]->upper->size();
								int subsizelo = longEdges[CurEdge]->lower->size();
								if (subsizeup>readLength-j-1-shift) subsizeup = readLength-j-1-shift;
								if (subsizelo>readLength-j-l) subsizelo = readLength-j-l;
								if ((longEdges[CurEdge]->upper->similar(upRead->Subseq(j+shift+1,readLength),subsizeup,LEFT))&&
									(longEdges[CurEdge]->lower->similar(loRead->Subseq(j+1,readLength),subsizelo,LEFT))){
									if (tmpOut ==-1){
										tmpOut = CurEdge;
									}
									else {
										tmpOut = -1;
										break;
									}
								}
							}
							if (tmpOut != -1){
								map<int, vector<pair<int,int>>>::iterator epIter = EdgePairs.find(VertId);
								if (epIter == EdgePairs.end()){
									vector<pair<int,int>> pairVect;
									pairVect.pb(make_pair(tmpIn,tmpOut));
									EdgePairs.insert(mp(VertId,pairVect));
									cerr<<"vert "<<VertId<<" "<<tmpIn<<" "<<tmpOut<<endl;
								}
								else
								{
									if (find(epIter->second.begin(),epIter->second.end(),make_pair(tmpIn,tmpOut)) == epIter->second.end()){
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
	fclose(fout);
	forn(curVertId,VertexCount){
		if ((inD[curVertId]!=0)&&(outD[curVertId]!=0)) {
			cerr<<"Vertex "<<curVertId<<" connect edges:"<< endl;
			forn(i,(EdgePairs[curVertId]).size()){
				cerr<<(EdgePairs[curVertId])[i].first<<" ("<<longEdges[(EdgePairs[curVertId])[i].first]->FromVertex<<", "<<longEdges[(EdgePairs[curVertId])[i].first]->ToVertex<<")   ";
				cerr<<(EdgePairs[curVertId])[i].second<<" ("<<longEdges[(EdgePairs[curVertId])[i].second]->FromVertex<<", "<<longEdges[(EdgePairs[curVertId])[i].second]->ToVertex<<")"<<endl;
			}
		}
	}


	forn(curVertId,VertexCount){
		if ((inD[curVertId]!=0)&&(outD[curVertId]!=0))
		if ((inD[curVertId]<=EdgePairs[curVertId].size())&&(outD[curVertId]<=EdgePairs[curVertId].size()))
		{
			cerr<<"Process vertex "<<curVertId<<" IN "<<inD[curVertId]<<" OUT "<<outD[curVertId]<<" unique ways "<<EdgePairs[curVertId].size()<<endl;
			forn(i,(EdgePairs[curVertId]).size()){
				int CurIn = (EdgePairs[curVertId])[i].first;
				int CurOut = (EdgePairs[curVertId])[i].second;
				bool singlePair = true;
				forn(j,(EdgePairs[curVertId]).size()){
					if ((EdgePairs[curVertId])[j].first==CurIn)
						if ((EdgePairs[curVertId])[j].second!=CurOut) singlePair = false;
					if ((EdgePairs[curVertId])[j].second==CurOut)
						if ((EdgePairs[curVertId])[j].first!=CurIn) singlePair = false;
				}
				if (singlePair){
					cerr<<"Search in edge "<<CurIn<<" position"<<endl;
					while (longEdges[CurIn]->EdgeId !=CurIn) {
						cerr<<"Edge "<<CurIn<<" included into "<<longEdges[CurIn]->EdgeId<<endl;
						CurIn = longEdges[CurIn]->EdgeId;
					}
					cerr<<"Search out edge "<<CurOut<<" position"<<endl;
					while (longEdges[CurOut]->EdgeId !=CurOut){
						cerr<<"Edge "<<CurOut<<" included into "<<longEdges[CurOut]->EdgeId<<endl;
						CurOut = longEdges[CurOut]->EdgeId;
					}
					cerr<<"New inclusion "<<CurIn<<"("<<longEdges[CurIn]->FromVertex<<","<<longEdges[CurIn]->ToVertex<<") <- "<<CurOut<<"("<<longEdges[CurOut]->FromVertex<<","<<longEdges[CurOut]->ToVertex<<")"<<endl;

					longEdges[CurIn]->ExpandRight(*longEdges[CurOut]);
					longEdges[CurOut] = longEdges[CurIn];
					cerr<<"New edges "<<CurIn<<"("<<longEdges[CurIn]->FromVertex<<","<<longEdges[CurIn]->ToVertex<<") <- "<<CurOut<<"("<<longEdges[CurOut]->FromVertex<<","<<longEdges[CurOut]->ToVertex<<")"<<endl;
					inD[curVertId]--;
					outD[curVertId]--;
				}
			}
		}
	}
}






