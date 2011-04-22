#include "common.hpp"
#include "pairedGraph.hpp"
#include "graphSimplification.hpp"
#include "graphio.hpp"

LOGGER ("p.readTracing");

using namespace paired_assembler;

void traceReads(verticesMap &verts, longEdgesMap &longEdges,
		PairedGraph &graph, int &VertexCount, int &EdgeId) {

	edgePairsMap EdgePairs;
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
		if (fictiveSecondReads) {
			forn(i, readLength) {
				lowerNuclRead[i] = 0;
			}
		}
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
	SplitVertecesByEdgeConnections(graph, EdgePairs, true);
}
