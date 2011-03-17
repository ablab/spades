#include "common.hpp"
#include "pairedGraph.hpp"
#include "graphio.hpp"

using namespace paired_assembler;

void expandDefinite(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount) {
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
				graph.inputEdges[DestVertex][a]
						= graph.inputEdges[DestVertex][a + 1];
				a++;
			}
			graph.inD[DestVertex]--;
			forn(j,graph.inD[i]) {
				longEdges[graph.inputEdges[i][j]]->ExpandRight(
						*(longEdges[expandEdgeIndex]));
				graph.inputEdges[DestVertex][graph.inD[DestVertex]]
						= graph.inputEdges[i][j];
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
				graph.outputEdges[SourceVertex][a]
						= graph.outputEdges[SourceVertex][a + 1];
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

void traceReads(verticesMap &verts, longEdgesMap &longEdges,
		PairedGraph &graph, int &VertexCount, int &EdgeId) {

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
	//included in next step :)
	/*	forn(curVertId,VertexCount) {
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
	 */
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
	expandDefinite(longEdges, graph, VertexCount);
	for(longEdgesMap::iterator it= longEdges.begin(); it !=longEdges.end(); ++it) {
		if (it->second->FromVertex>=FakeVertexStart) it->second->FromVertex = FakeVertexToReal[it->second->FromVertex];
		if (it->second->ToVertex>=FakeVertexStart) it->second->ToVertex = FakeVertexToReal[it->second->ToVertex];
	}

	freopen("data/graph_after_fake.dot", "w", stdout);
	outputLongEdges(longEdges);

}

bool isPath(Edge &e1, Edge &e2) {
	if (e1.ToVertex != e2.FromVertex)
		return false;
	int sts = readLength + insertLength;
	if (e1.upper->size() + e2.upper->size() >= sts)
		return true;
	int left1 = max(0, e1.length - sts);
	int right1 = min(e1.upper->size(), e1.upper->size() - sts + e1.length);
	int left2 = left1 + sts - e1.length;
	int right2 = right1 + sts - e1.length;
	return e1.lower->Subseq(left1, right1) == e2.upper->Subseq(left2, right2);
}
