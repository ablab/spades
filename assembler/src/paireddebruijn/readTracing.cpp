#include "common.hpp"
#include "pairedGraph.hpp"
#include "graphSimplification.hpp"
#include "graphio.hpp"

LOGGER ("p.readTracing");

using namespace paired_assembler;

int edgeIdToLocalId(int dir, PairedGraph &graph, int edgeId) {
	int vertId;
	edgeId = edgeRealId( edgeId , graph.longEdges);
	if (dir == 0)
		vertId = graph.longEdges[edgeId]->ToVertex;
	else
		vertId = graph.longEdges[edgeId]->FromVertex;
//	TRACE("edge connects "<< graph.longEdges[edgeId]->ToVertex << " " << graph.longEdges[edgeId]->FromVertex);
	forn(i, MAX_DEGREE) {
//		TRACE("in edgeIdTo.." << graph.edgeIds[vertId][i][dir] << "  " << edgeId);
		if (graph.edgeIds[vertId][i][dir] == edgeId)
			return i;
	}
	return -1;
}
void traceReads(verticesMap &verts, longEdgesMap &longEdges,
		PairedGraph &graph, int &VertexCount, int &EdgeId) {

	map<int, vector<pair<int, int>>> EdgePairs;
	int table[MAX_DEGREE][MAX_DEGREE];
	double fake_cov[MAX_DEGREE][MAX_DEGREE];
	double in_cov[MAX_DEGREE];
	double out_cov[MAX_DEGREE];
	double eps_out_cov[MAX_DEGREE];
	double eps = 10000;
	INFO("traceReads started");
	INFO(parsed_reads);
	FILE * inFile = fopen(parsed_reads.c_str(), "r");
	char *upperNuclRead = new char[readLength + 2];
	char *lowerNuclRead = new char[readLength + 2];
	char *upperRead = new char[readLength + 2];
	char *lowerRead = new char[readLength + 2];
	ll count=0;
	ll upperMask = (((ll) 1) << (2 * (k - 1))) - 1;
	ll lowerMask = (((ll) 1) << (2 * (l - 1))) - 1;
	//	FILE* fout = fopen("data/filtered_reads","w");
	while (nextReadPair(inFile, upperNuclRead, lowerNuclRead)) {
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
					cerr<<"kmer found for j="<<j<<endl;
				for (vector<VertexPrototype *>::iterator it =vertIter->second.begin(); it
						!= vertIter->second.end(); ++it) {

					if ((*it)->lower->similar(loRead->Subseq(1+j, l+j),l-1)) {
										cerr<<"vertex found for lower "<<(*it)->lower->str()<<endl;
//												fprintf(fout,"%s %s\n",upperNuclRead,lowerNuclRead);
						int VertId = (*it)->VertexId;
						pair<int, int> vDist = vertexDist(longEdges, graph,VertId);

						if ((graph.degrees[VertId][0]!=0)&&(graph.degrees[VertId][1]!=0)) {
							int tmpIn = -1;
							forn(i,graph.degrees[VertId][0]) {
								int CurEdge=graph.edgeIds[VertId][i][IN_EDGE];
								int subsizeup = longEdges[CurEdge]->upper->size();
								int subsizelo = longEdges[CurEdge]->lower->size();
								if (subsizeup>j+k+shift) subsizeup = j+k+shift;
								if (subsizelo>j+l) subsizelo = j+l;
								int startup;
								int startlo;
								//								if (vDist.first+vDist.second < readLength+k){
								//									cerr<<"CheckUp "<<longEdges[CurEdge]->upper->Subseq(longEdges[CurEdge]->upper->size()-subsizeup,longEdges[CurEdge]->upper->size()).str()<<" VS "<<upRead->Subseq(j+k+shift-subsizeup,j+k+shift).str()<<" with "<<subsizeup<<endl;
								//									cerr<<"CheckLo "<<longEdges[CurEdge]->lower->Subseq(longEdges[CurEdge]->lower->size()-subsizelo,longEdges[CurEdge]->lower->size()).str()<<" VS "<<loRead->Subseq(j+l-subsizelo,j+l).str()<<" with "<<subsizelo<<endl;
								//								}
								if ((longEdges[CurEdge]->upper->Subseq(longEdges[CurEdge]->upper->size()-subsizeup,longEdges[CurEdge]->upper->size()).similar(upRead->Subseq(j+k+shift-subsizeup,j+k+shift),subsizeup,RIGHT))&&
										(longEdges[CurEdge]->lower->Subseq(longEdges[CurEdge]->lower->size()-subsizelo,longEdges[CurEdge]->lower->size()).similar(loRead->Subseq(j+l-subsizelo,j+l),subsizelo,RIGHT))) {
									//									cerr<<"in EQU ";
									if (tmpIn ==-1) {
										//										cerr<<"first time"<<endl;
										tmpIn = CurEdge;
									}
									else {
										//cerr<<"second time"<<endl;
										tmpIn = -1;
										break;
									}
								}
							}
							if (tmpIn == -1) break;
							int tmpOut = -1;

							forn(i,graph.degrees[VertId][1]) {
								int CurEdge=graph.edgeIds[VertId][i][OUT_EDGE];
								int subsizeup = longEdges[CurEdge]->upper->size();
								int subsizelo = longEdges[CurEdge]->lower->size();
								if (subsizeup>readLength-j-1-shift) subsizeup = readLength-j-1-shift;
								if (subsizelo>readLength-j-1) subsizelo = readLength-j-1;
								//								if (vDist.first+vDist.second < readLength+k){
								//									cerr<<"CheckUp "<<longEdges[CurEdge]->upper->Subseq(0,subsizeup).str()<<" VS "<<upRead->Subseq(j+shift+1,j+shift+1+subsizeup).str()<<" with "<<subsizeup<<endl;
								//									cerr<<"CheckLo "<<longEdges[CurEdge]->lower->Subseq(0,subsizelo).str()<<" VS "<<loRead->Subseq(j+1,j+1+subsizelo).str()<<" with "<<subsizelo<<endl;
								//								}
								if ((longEdges[CurEdge]->upper->Subseq(0,subsizeup).similar(upRead->Subseq(j+shift+1,j+shift+1+subsizeup),subsizeup,LEFT))&&
										(longEdges[CurEdge]->lower->Subseq(0,subsizelo).similar(loRead->Subseq(j+1,j+1+subsizelo),subsizelo,LEFT))) {
									//									cerr<<"out EQU ";
									if (tmpOut ==-1) {
										//										cerr<<"first time"<<endl;
										tmpOut = CurEdge;
									}
									else {
										//										cerr<<"second time"<<endl;
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
	INFO ("reads processed, started tracing");
	TRACE ("Fake Tracing");
	forn(curVertId,VertexCount) {
		if ((graph.degrees[curVertId][0]!=0)&&(graph.degrees[curVertId][1]!=0)) {
			cerr<<"Vertex "<<curVertId<<" connect edges:"<< endl;
			forn(i,(EdgePairs[curVertId]).size()) {
				cerr<<(EdgePairs[curVertId])[i].first<<" ("<<longEdges[(EdgePairs[curVertId])[i].first]->FromVertex<<", "<<longEdges[(EdgePairs[curVertId])[i].first]->ToVertex<<")   ";
				cerr<<(EdgePairs[curVertId])[i].second<<" ("<<longEdges[(EdgePairs[curVertId])[i].second]->FromVertex<<", "<<longEdges[(EdgePairs[curVertId])[i].second]->ToVertex<<")"<<endl;
			}
		}
	}
	/*
	 //resolving simple case InD=OutD and all edge is unique
	 //included in next step :)

	 forn(curVertId,VertexCount) {
	 if ((graph.inD[curVertId]!=0)&&(graph.degrees[curVertId][1!=0))
	 if ((graph.inD[curVertId]==EdgePairs[curVertId].size())&&(graph.degrees[curVertId][1]==EdgePairs[curVertId].size()))
	 {
	 cerr<<"Process vertex "<<curVertId<<" IN "<<graph.inD[curVertId]<<" OUT "<<graph.degrees[curVertId][1<<" unique ways "<<EdgePairs[curVertId].size()<<endl;
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
	 graph.degrees[curVertId][1]--;
	 }
	 }
	 }
	 }

	 outputLongEdges(longEdges);
	 */
	//resolve multi case;
	int FakeVertexCount = VertexCount;
	int FakeVertexStart = VertexCount;
	map<int, int> FakeVertexToReal;
	forn(curVertId,VertexCount) {
		if ((graph.degrees[curVertId][0]!=0)&&(graph.degrees[curVertId][1]!=0)) {
			pair<int, int> vDist = vertexDist(longEdges,graph,curVertId);
			//			if (vDist.first+vDist.second-k+1<=readLength)
			{
				cerr<<"vertex "<<curVertId<<" dist <"<<vDist.first<<", "<<vDist.second<<">"<<endl;
				if ((graph.degrees[curVertId][0]<=EdgePairs[curVertId].size())&&(graph.degrees[curVertId][1]<=EdgePairs[curVertId].size())) {
					bool allIns = false;
					bool allOuts = false;
					forn(i,graph.degrees[curVertId][0]) {
						allIns = false;
						forn (j,EdgePairs[curVertId].size()) {
							if ((EdgePairs[curVertId])[j].first ==graph.edgeIds[curVertId][i][IN_EDGE]) {
								allIns = true;
								break;
							}
						}
						if (!allIns) break;
					}
					forn(i,graph.degrees[curVertId][1]) {
						allOuts = false;
						forn (j,EdgePairs[curVertId].size()) {
							if ((EdgePairs[curVertId])[j].second ==graph.edgeIds[curVertId][i][OUT_EDGE]) {
								allOuts = true;
								break;
							}
						}
						if (!allOuts) break;
					}
					if (allIns&&allOuts) {
						TRACE ("Faking vert " << curVertId);
						int tmpCurOut = edgeRealId(graph.edgeIds[curVertId][0][OUT_EDGE],longEdges);
						string tmpUpSeq = longEdges[tmpCurOut]->upper->Subseq(0,k-1).str();
						string tmpLoSeq = longEdges[tmpCurOut]->lower->Subseq(0,l-1).str();

						//create FakeVerices and Fake edges;

						memset(table, 0, sizeof(table));
						memset(fake_cov, 0, sizeof(fake_cov));
						TRACE ("there are "<< EdgePairs[curVertId].size() << " pairs");
						forn (tmp, EdgePairs[curVertId].size()){
	//						TRACE((EdgePairs[curVertId])[tmp].first);
	//						TRACE(edgeRealId((EdgePairs[curVertId])[tmp].first, graph.longEdges));
							TRACE(edgeIdToLocalId(0, graph, (EdgePairs[curVertId])[tmp].first) <<" "<< edgeIdToLocalId(1, graph, (EdgePairs[curVertId])[tmp].second)<< " "<<EdgePairs[curVertId][tmp].first << " " << EdgePairs[curVertId][tmp].second);
							table[edgeIdToLocalId(0, graph, (EdgePairs[curVertId])[tmp].first)][edgeIdToLocalId(1, graph, (EdgePairs[curVertId])[tmp].second)] = 1;
							//fake_deg[tmp.second] ++;
						}
						//create fake Vertices for in edges;
						int tmpFictStartIn = FakeVertexCount;
						memset(in_cov, 0, sizeof(in_cov));
						forn(i,graph.degrees[curVertId][0]) {
							int CurIn = edgeRealId(graph.edgeIds[curVertId][i][IN_EDGE],longEdges);
							graph.edgeIds[curVertId][i][IN_EDGE] = CurIn;
							FakeVertexToReal.insert(make_pair(FakeVertexCount,curVertId));
							cerr<<"Edge "<<CurIn<<" ("<<longEdges[CurIn]->FromVertex<<", "<<longEdges[CurIn]->ToVertex<<") maped to ";
							longEdges[CurIn]->ToVertex = FakeVertexCount;
							graph.degrees[FakeVertexCount][0]=1;
							graph.degrees[FakeVertexCount][1]=0;
							in_cov[i] = graph.longEdges[CurIn]->coverage;
							graph.edgeIds[FakeVertexCount][0][IN_EDGE]=CurIn;
							cerr<<"edge "<<CurIn<<" ("<<longEdges[CurIn]->FromVertex<<", "<<longEdges[CurIn]->ToVertex<<")"<<endl;
							FakeVertexCount++;
						}



						//create fake Vertices for out edges;
						int tmpFictStartOut = FakeVertexCount;
						memset(out_cov, 0, sizeof(out_cov));
						forn(i,graph.degrees[curVertId][1]) {
							int CurOut = edgeRealId(graph.edgeIds[curVertId][i][OUT_EDGE],longEdges);
							graph.edgeIds[curVertId][i][OUT_EDGE] = CurOut;
							cerr<<"Edge "<<CurOut<<" ("<<longEdges[CurOut]->FromVertex<<", "<<longEdges[CurOut]->ToVertex<<") maped to ";
							FakeVertexToReal.insert(make_pair(FakeVertexCount,curVertId));
							longEdges[CurOut]->FromVertex = FakeVertexCount;
							graph.degrees[FakeVertexCount][0]=0;
							graph.degrees[FakeVertexCount][1]=1;
							graph.edgeIds[FakeVertexCount][0][OUT_EDGE]=CurOut;

							out_cov[i] = graph.longEdges[CurOut]->coverage;
							cerr<<"edge "<<CurOut<<" ("<<longEdges[CurOut]->FromVertex<<", "<<longEdges[CurOut]->ToVertex<<")"<<endl;
							FakeVertexCount++;
						}


						//compute fake coverage

						TRACE("in_cov:")
						forn(i, MAX_DEGREE){
							TRACE(in_cov[i]);
						}
						TRACE("out_cov:")
						forn(i, MAX_DEGREE){
							TRACE(out_cov[i]);
						}

						TRACE("FAKE_EDGES");
						forn(i, MAX_DEGREE)
							forn(j, MAX_DEGREE)
								if (table[i][j] > eps)
									TRACE(i << " "<< j);
						while (1) {
							double sum = 0;
							forn (i, MAX_DEGREE) {
								sum += out_cov[i];
							}
							if (sum < eps) break;
							TRACE("SUM: "<<sum);
							forn (i, MAX_DEGREE) {
								eps_out_cov[i] = out_cov[i] * (eps / sum) ;
								out_cov[i] -= eps_out_cov[i];
							}
							forn(j, MAX_DEGREE) {
								double tdeg = 0;
								forn(i, MAX_DEGREE)
									if (table[i][j])
										tdeg += in_cov[i];
								if (tdeg >eps/10)
								forn(i, MAX_DEGREE) {
									if (table[i][j]) {
//										assert(tdeg>eps/10);
										double add = (in_cov[i]/tdeg) * eps_out_cov[j];
										if (add > in_cov[i]) add = in_cov[i];
										fake_cov[i][j] +=add;
										in_cov[i] -= add;
									}
								}
							}
						}
						TRACE("FAKE_COV");
						forn(i, MAX_DEGREE)
							forn(j, MAX_DEGREE)
								if (fake_cov[i][j] > eps/10)
									TRACE(i << " "<< j << " " << fake_cov[i][j]);
						//create fake edges
						TRACE("FAKE_COV traced");
						forn (tmpEdgePair,EdgePairs[curVertId].size()) {
							int tmpFrom = 0;
							int tmpTo = 0;
							while (edgeRealId(graph.edgeIds[curVertId][tmpFrom][IN_EDGE], longEdges)!=edgeRealId((EdgePairs[curVertId])[tmpEdgePair].first,longEdges)) {
								tmpFrom++;
								assert(tmpFrom<graph.degrees[curVertId][0]);
							}
							while (edgeRealId(graph.edgeIds[curVertId][tmpTo][OUT_EDGE],longEdges)!=edgeRealId((EdgePairs[curVertId])[tmpEdgePair].second, longEdges)) {
								tmpTo++;
								assert(tmpTo<graph.degrees[curVertId][1]);
							}
							Sequence *UpSeq = new Sequence(tmpUpSeq);
							Sequence *LoSeq = new Sequence(tmpLoSeq);

//							Edge *tmpEdge = new Edge(UpSeq, LoSeq, tmpFictStartIn + tmpFrom, tmpFictStartOut + tmpTo,0, EdgeId, fake_cov[edgeIdToLocalId(0, graph, EdgePairs[curVertId][tmpEdgePair].first)][edgeIdToLocalId(1, graph, EdgePairs[curVertId][tmpEdgePair].second)]);
							Edge *tmpEdge = new Edge(UpSeq, LoSeq, tmpFictStartIn + tmpFrom, tmpFictStartOut + tmpTo,0, EdgeId, fake_cov[tmpFrom][tmpTo]);
							longEdges.insert(make_pair(EdgeId,tmpEdge));
							//							cerr<<"Virtual edge "<<EdgeId<<
							cerr<<"Virtual edge "<<EdgeId<<" ("<<longEdges[EdgeId]->FromVertex<<", "<<longEdges[EdgeId]->ToVertex<<") ["<<tmpFrom<<", "<<tmpTo<<"] cov "<< fake_cov[tmpFrom][tmpTo]<<endl;
							graph.edgeIds[tmpFictStartIn + tmpFrom][graph.degrees[tmpFictStartIn + tmpFrom][1]][OUT_EDGE] = EdgeId;
							graph.degrees[tmpFictStartIn + tmpFrom][1]++;
							graph.edgeIds[tmpFictStartOut + tmpTo][graph.degrees[tmpFictStartOut + tmpTo][0]][OUT_EDGE] = EdgeId;
							graph.degrees[tmpFictStartOut + tmpTo][0]++;
							EdgeId++;

						}

						graph.degrees[curVertId][0] = 0;
						graph.degrees[curVertId][1] = 0;
					}
				}
			}
		}
	}

	VertexCount = FakeVertexCount;
	graph.recreateVerticesInfo(VertexCount, longEdges);
	expandDefinite(longEdges, graph, VertexCount);
//	for(longEdgesMap::iterator it= longEdges.begin(); it !=longEdges.end(); ++it) {
//		if (it->second->FromVertex>=FakeVertexStart) it->second->FromVertex = FakeVertexToReal[it->second->FromVertex];
//		if (it->second->ToVertex>=FakeVertexStart) it->second->ToVertex = FakeVertexToReal[it->second->ToVertex];
//	}
//	VertexCount = FakeVertexStart;
	graph.recreateVerticesInfo(VertexCount, longEdges);
}
