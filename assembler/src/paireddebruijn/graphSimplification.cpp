/*
 * graphSimplification.cpp
 *
 *  Created on: 17.03.2011
 */
#include "graphSimplification.hpp"
#include "common.hpp"
#include "constructHashTable.hpp"

DECL_MODULE_LOGGER("graphSimplification")


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

bool isPath(Edge &e1, Edge &e2) {
	cerr << endl << "(" << e1.FromVertex << ", " << e1.ToVertex << ") " << "("
			<< e2.FromVertex << ", " << e2.ToVertex << ") ";
	if (e1.ToVertex != e2.FromVertex)
		return false;
	int sts = readLength + insertLength;
	if (e1.length + e2.length + k - 1 <= sts)
		return true;
	int left1 = max(0, e1.length - sts);
	int right1 = min(e1.upper->size(), e1.upper->size() - sts + e2.length);
	int left2 = left1 + sts - e1.length;
	int right2 = right1 + sts - e1.length;
	cerr << endl << e1.lower->Subseq(left1, right1).str() << endl
			<< e2.upper->Subseq(left2, right2).str();
	return e1.lower->Subseq(left1, right1) == e2.upper->Subseq(left2, right2);
}

bool isRealPossiblePath(PairedGraph &graph, int &EdgeId, int &FollowerId) {
	int sts = readLength + insertLength;
	bool res =isPath(*(graph.longEdges[EdgeId]), *(graph.longEdges[FollowerId]));
	if (!res) return false;
	if (graph.longEdges[EdgeId]->length>sts) {
		return res;
	}
	else {
		Edge* tmpEdge = new Edge(*(graph.longEdges[EdgeId]));
		tmpEdge->ExpandRight(*(graph.longEdges[FollowerId]));

		PairThreader pg(graph,1);
		vector<pair<int, Edge *> > vp = pg.threadLower(tmpEdge);
		cerr<<"Edge "<<EdgeId<<" + "<<FollowerId<<" extended by "<<vp.size()<<" edges:"<<endl;
		forn(i, vp.size()) {
			cerr<<vp[i].second->EdgeId<<" dist "<<vp[i].first<<endl;
			isPath(*tmpEdge, *(graph.longEdges[vp[i].second->EdgeId]));
		}
		delete tmpEdge;
		return (vp.size()>0);

	}

	return true;

}

pair<bool, int> isPath(Edge *e1, Edge *e2, int shift) {
	int lowerLeft = max(0, shift);
	int lowerRight = min((int) (e1->length), e2->length + shift) + k - 1;
	int upperLeft = max(0, -shift);
	int upperRight = min(e2->length, e1->length - shift) + k - 1;
	Sequence lowerSequence(e1->lower->Subseq(lowerLeft, lowerRight));
	Sequence upperSequence(e2->upper->Subseq(upperLeft, upperRight));
	return make_pair(lowerSequence == upperSequence, lowerRight - lowerLeft);
}

bool addEntry(vector<pair<int, Edge *> > &result, int distance,
		Edge *edgeToAdd) {
	if (distance < 0) {
		return false;
	}
	for (vector<pair<int, Edge*> >::iterator it = result.begin(); it
			!= result.end(); ++it) {
		if(it->first == distance && it->second == edgeToAdd) {
			return false;
		}
	}
	result.push_back(make_pair(distance, edgeToAdd));
	return true;
}

void PairThreader::threadLower(vector<pair<int, Edge *> > &result,
		Edge *currentEdge, int shift, Edge *start) {
	if (shift >= start->length)
		return;
	if (shift + currentEdge->length >= 0) {
		pair<bool, int> intersection = isPath(start, currentEdge, shift);
		if (!intersection.first)
			return;
		if (intersection.second >= minIntersection_) {
			int distance = shift + insertLength + readLength - start->length;
			addEntry(result, distance, currentEdge);
		}
	}
	for (int i = 0; i < g_.degrees[currentEdge->ToVertex][1]; i++) {
		Edge *nextEdge = g_.longEdges[g_.edgeIds[currentEdge->ToVertex][i][1]];
		threadLower(result, nextEdge, shift + currentEdge->length, start);
	}
}

vector<pair<int, Edge *> > PairThreader::threadLower(Edge *start) {
	vector<pair<int, Edge *> > result;
	threadLower(result, start, -insertLength - readLength, start);
	return result;
}

bool processLowerSequence(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount) {
	int countIn[MAX_DEGREE], possibleIn[MAX_DEGREE], countOut[MAX_DEGREE],
			possibleOut[MAX_DEGREE];
	bool res = false;
	forn(curVertId, VertexCount) {
		if ((graph.degrees[curVertId][0] != 0) && (graph.degrees[curVertId][1]
				!= 0)) {
			cerr<<"Process vertex "<<curVertId<<endl;
			forn(i,MAX_DEGREE) {
				countOut[i] = 0;
				countIn[i] = 0;
			}
			forn(i,graph.degrees[curVertId][0]) {
				int curInEdgeId = edgeRealId(
						graph.edgeIds[curVertId][i][IN_EDGE], longEdges);
				graph.edgeIds[curVertId][i][IN_EDGE] = curInEdgeId;
				forn(j,graph.degrees[curVertId][1]) {
					int curOutEdgeId = edgeRealId(
							graph.edgeIds[curVertId][j][OUT_EDGE], longEdges);
					graph.edgeIds[curVertId][j][OUT_EDGE] = curOutEdgeId;
					cerr << "Check isRealPossiblePass for edge " << curInEdgeId << " vs "
							<< curOutEdgeId<<endl;
					if (isRealPossiblePath(graph, curInEdgeId, curOutEdgeId)) {
						countIn[j]++;
						countOut[i]++;
						possibleIn[j] = i;
						possibleOut[i] = j;
						cerr << " POSSIBLE" << endl;
					} else
						cerr << " IMPOSSIBLE" << endl;
				}
			}

			bool AllInHasDefiniteOut = true;
			bool AllOutHasDefiniteIn = true;
			forn(i,graph.degrees[curVertId][0]) {
				if (countOut[i] != 1)
					AllInHasDefiniteOut = false;
			}

			forn(i,graph.degrees[curVertId][1]) {
				if (countIn[i] != 1)
					AllOutHasDefiniteIn = false;
			}

			if (AllInHasDefiniteOut) {
				//				res = true;
				cerr << "Vert " << curVertId
						<< " resolvable: AllInHasDefiniteOut" << endl;
			}
			//			else

			if (AllOutHasDefiniteIn) {
				//			res = true;
				cerr << "Vert " << curVertId
						<< " resolvable: AllOutHasDefiniteIn" << endl;
			}
			//			else
			forn(i,graph.degrees[curVertId][0]) {
				if (countOut[i] == 1) {
					if (countIn[possibleOut[i]] == 1) {
						int InEdge =
								edgeRealId(
										graph.edgeIds[curVertId][i][IN_EDGE],
										longEdges);
						int
								OutEdge =
										edgeRealId(
												graph.edgeIds[curVertId][possibleOut[i]][OUT_EDGE],
												longEdges);
						longEdges[InEdge]->ExpandRight(*longEdges[OutEdge]);
						longEdges[OutEdge] = longEdges[InEdge];
						res = true;
					}
				}
			}
			//REMOVE VERTICES INFO!!! //I did it outside
		}
	}
	return res;
}

/*
 bool commonSequencesExtraction(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount){
 bool res = false;
 forn(curVertexId, VertexCount){


 }
 return res;
 }
 */
/* countEdgesDistintion
 * @params:
 * countEdgesDistinction
 * @param vertexId - index of vertex, distances we are interested in
 * @param graph - pairedGraph
 * @returns pair<int,int> - minimal intervals of nucleotids, such that they uniquely represent Inedges and outedges respectively.
 * If one of results  is more than insertLength + 2 * readLength, returns insertLength + 2 * readLength + 1
 *
 */
pair<int, int> vertexDist(longEdgesMap &longEdges, PairedGraph &graph,
		int vertexId) {
	size_t res1 = k - 1;
	size_t res2 = k - 1;
	size_t max_res = insertLength + 2 * readLength + 1;
	//int colors[MAX_DEGREE];
	forn(i, graph.degrees[vertexId][0]) {
		int e1 = graph.edgeIds[vertexId][i][IN_EDGE];
		Sequence tmp = *longEdges[e1]->upper;
		Sequence tmpl = *longEdges[e1]->lower;

		//		cerr << endl << "1 " << tmp.str();
		forn(j, i) {
			size_t count = 0;

			int e2 = graph.edgeIds[vertexId][j][IN_EDGE];
			while (count < max_res) {
				Sequence tmp2 = *longEdges[e2]->upper;
				Sequence tmp2l = *longEdges[e2]->lower;

				//				cerr << endl <<"2 "<< tmp2.str();
				//				cerr << endl <<count << " ";// << nucl(tmp[tmp.size() - 1 - count]) << " " << nucl (tmp2[tmp2.size() - 1 -count]);
				if (count >= min(tmp.size(), tmp2.size()) || tmp[tmp.size() - 1
						- count] != tmp2[tmp2.size() - 1 - count])
					//				if (*(longEdges[e1]->upper)[count] != *(longEdges[e2]->upper)[count])
					break;
				if (count >= min(tmpl.size(), tmp2l.size()) || tmpl[tmpl.size()
						- 1 - count] != tmp2l[tmp2l.size() - 1 - count])
					//				if (*(longEdges[e1]->upper)[count] != *(longEdges[e2]->upper)[count])
					break;

				count++;
			}
			res1 = max(res1, count);
		}
	}
	forn(i, graph.degrees[vertexId][1]) {
		int e1 = graph.edgeIds[vertexId][i][OUT_EDGE];
		Sequence tmp = *longEdges[e1]->upper;
		Sequence tmpl = *longEdges[e1]->lower;

		//		cerr << endl << "1 " << tmp.str();
		forn(j, i) {
			int e2 = graph.edgeIds[vertexId][j][OUT_EDGE];
			size_t count = 0;
			while (count < max_res) {
				Sequence tmp2 = *longEdges[e2]->upper;
				Sequence tmp2l = *longEdges[e2]->lower;

				//				cerr << endl <<"2 "<< tmp2.str();
				//				cerr << endl <<count << " ";// << nucl(tmp[tmp.size() - 1 - count]) << " " << nucl (tmp2[tmp2.size() - 1 -count]);
				if (count >= min(tmp.size(), tmp2.size()) || tmp[count]
						!= tmp2[count])
					//				if (*(longEdges[e1]->upper)[count] != *(longEdges[e2]->upper)[count])
					break;
				if (count >= min(tmpl.size(), tmp2l.size()) || tmpl[count]
						!= tmp2l[count])
					//				if (*(longEdges[e1]->upper)[count] != *(longEdges[e2]->upper)[count])
					break;

				count++;
			}
			res2 = max(res2, count);
		}
	}

	//	cerr << endl << res1 << " " << res2<<endl;
	//assert(0);
	return make_pair(res1 + 1, res2 + 1);
}


void cutShortTips(PairedGraph &graph, int MaxCutLength){
	cerr<<"Start clip tips: max vert"<<graph.VertexCount<<endl;
	forn(i,graph.VertexCount) {
		if ((graph.rightDegree(i) == 1) && (graph.leftDegree(i) == 0)){
			Edge* curEdge = graph.neighbourEdge(i,0, RIGHT);
			cerr<<"cut tips vert "<<i<<" edge "<< curEdge->EdgeId<<endl;
			if ((curEdge->length<=MaxCutLength)&&((curEdge->coverage<=200))) graph.removeEdge(curEdge);
		} else
		if ((graph.rightDegree(i) == 0) && (graph.leftDegree(i) == 1)){
			Edge* curEdge = graph.neighbourEdge(i,0, LEFT);
			cerr<<"cut tips vert "<<i<<" edge "<< curEdge->EdgeId<<endl;
			if ((curEdge->length<=MaxCutLength)&&(((curEdge->coverage<=200)))) graph.removeEdge(curEdge);
		}
	}
	graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
}


void expandObvious(PairedGraph &graph) {
	int expandEdgeIndex;
	longEdgesMap::iterator it;
	cerr << "expandObviousStart" << endl;
	forn(i,graph.VertexCount) {

		if ((graph.degrees[i][1] == 1) && (graph.degrees[i][0] == 1)) {
			cerr << "expand obvious Vertex " <<i<< endl;
			expandEdgeIndex = edgeRealId(graph.edgeIds[i][0][OUT_EDGE],	graph.longEdges);
			if (expandEdgeIndex != graph.edgeIds[i][0][OUT_EDGE]) {
				cerr<<" BAD edge "<<graph.edgeIds[i][0][OUT_EDGE]<<" rId "<<expandEdgeIndex<<endl;
			}
			int DestVertex = graph.longEdges[expandEdgeIndex]->ToVertex;
			if (DestVertex == i) {
				WARN("Expand obvious has bad loop");
				continue;
			}

			int a = 0;

			while ((edgeRealId(graph.edgeIds[DestVertex][a][IN_EDGE],
					graph.longEdges) != expandEdgeIndex)&&(a < graph.degrees[DestVertex][0])){
				cerr<<"vertex "<<DestVertex<<" in edge "<<edgeRealId(graph.edgeIds[DestVertex][a][IN_EDGE],
						graph.longEdges)<<" on position "<<a<<endl;
				a++;
			}
			cerr<<"Total in edges "<<graph.degrees[DestVertex][0]<<" but we want "<<a<<endl;
			assert(a < graph.degrees[DestVertex][0]);
			//				cerr << a;
			while (a < graph.degrees[DestVertex][0] - 1) {
				graph.edgeIds[DestVertex][a][IN_EDGE]
						= graph.edgeIds[DestVertex][a + 1][IN_EDGE];
				a++;
			}
				//				cerr << "hm"<<" "<< a << endl << graph.degrees[i][0]<< endl;
			graph.degrees[DestVertex][0]--;
			//				assert(graph.degrees[DestVertex][0] > 0);
			forn(j, graph.degrees[i][0]) {
				graph.longEdges[graph.edgeIds[i][j][IN_EDGE]]->ExpandRight(
						*(graph.longEdges[expandEdgeIndex]));
				graph.edgeIds[DestVertex][graph.degrees[DestVertex][0]][IN_EDGE]
						= graph.edgeIds[i][j][IN_EDGE];
				graph.degrees[DestVertex][0]++;
			}
			it = graph.longEdges.find(graph.edgeIds[i][0][OUT_EDGE]);
			graph.longEdges.erase(it);
			graph.degrees[i][0] = 0;
			graph.degrees[i][1] = 0;
		}
	}
	cerr << "expandObviousFinished" << endl;
}
void expandDefinite(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount, bool NotExpandBeyondDefinite) {
	longEdgesMap::iterator it;
	int expandEdgeIndex;
	cerr << "expandDefiniteStart" << endl;
	forn(i,VertexCount) {
		cerr << "expand definite right Vertex " <<i<< endl;

		if ((graph.degrees[i][1] == 1) && (graph.degrees[i][0] > 0)) {
			expandEdgeIndex = edgeRealId(graph.edgeIds[i][0][OUT_EDGE],
					longEdges);
			int DestVertex = longEdges[expandEdgeIndex]->ToVertex;
			if (DestVertex == i) {
				WARN("Expand definite right has bad loop");
				continue;
			}
			pair<int, int> diffDistDest = make_pair(0, 0);
			pair<int, int> diffDistCur = make_pair(0, 0);
			if (NotExpandBeyondDefinite) {
				diffDistDest = vertexDist(longEdges, graph, DestVertex);
				diffDistCur = vertexDist(longEdges, graph, i);
			}
			if ((!NotExpandBeyondDefinite) || (diffDistCur.first
					+ diffDistDest.second + longEdges[expandEdgeIndex]->length
					< readLength + k) ||(graph.degrees[i][0]==1)) {
				if (NotExpandBeyondDefinite)
					cerr << "Check cur vert " << i << " dest vert "
							<< DestVertex << "  " << diffDistCur.first << " + "
							<< diffDistDest.second << " + "
							<< longEdges[expandEdgeIndex]->length << " = "
							<< diffDistCur.first + diffDistDest.second
									+ longEdges[expandEdgeIndex]->length
							<< " < " << readLength + k << endl;
				int a = 0;

				while ((edgeRealId(graph.edgeIds[DestVertex][a][IN_EDGE],
						longEdges) != expandEdgeIndex)&&(a < graph.degrees[DestVertex][0])){
					cerr<<"vertex "<<DestVertex<<" in edge "<<edgeRealId(graph.edgeIds[DestVertex][a][IN_EDGE],
							longEdges)<<" on position "<<a<<endl;
					a++;
				}
				cerr<<"Total in edges "<<graph.degrees[DestVertex][0]<<" but we want "<<a<<endl;
				assert(a < graph.degrees[DestVertex][0]);
				//				cerr << a;
				while (a < graph.degrees[DestVertex][0] - 1) {
					graph.edgeIds[DestVertex][a][IN_EDGE]
							= graph.edgeIds[DestVertex][a + 1][IN_EDGE];
					a++;
				}

				//				cerr << "hm"<<" "<< a << endl << graph.degrees[i][0]<< endl;
				graph.degrees[DestVertex][0]--;
				//				assert(graph.degrees[DestVertex][0] > 0);
				forn(j, graph.degrees[i][0]) {
					longEdges[graph.edgeIds[i][j][IN_EDGE]]->ExpandRight(
							*(longEdges[expandEdgeIndex]));
					graph.edgeIds[DestVertex][graph.degrees[DestVertex][0]][IN_EDGE]
							= graph.edgeIds[i][j][IN_EDGE];
					graph.degrees[DestVertex][0]++;
				}
				it = longEdges.find(expandEdgeIndex);
				longEdges.erase(it);

				graph.degrees[i][0] = 0;
				graph.degrees[i][1] = 0;
				diffDistDest = vertexDist(longEdges, graph, DestVertex);
				if (NotExpandBeyondDefinite)
					cerr << "Must be: dest vert " << DestVertex << "  "
							<< diffDistDest.first << " + "
							<< diffDistDest.second << " = "
							<< diffDistDest.first + diffDistDest.second
							<< " < " << readLength + k << endl;

			}
		}
	}

	cerr << "expandDefinite second attempt\n" << endl;
	forn(i,VertexCount) {
		cerr << "expand definite left Vertex " <<i<< endl;

		if ((graph.degrees[i][0] == 1) && (graph.degrees[i][1] > 0)) {
			cerr << i << endl;
			expandEdgeIndex = edgeRealId(graph.edgeIds[i][0][IN_EDGE],
					longEdges);
			int SourceVertex = longEdges[expandEdgeIndex]->FromVertex;
			if (SourceVertex == i) {
				WARN("Expand definite left has bad loop");
				continue;
			}
			pair<int, int> diffDistSource = make_pair(0, 0);
			pair<int, int> diffDistCur = make_pair(0, 0);
			if (NotExpandBeyondDefinite) {
				diffDistSource = vertexDist(longEdges, graph, SourceVertex);
				diffDistCur = vertexDist(longEdges, graph, i);
			}

			if ((!NotExpandBeyondDefinite) || (diffDistCur.second
					+ diffDistSource.first + longEdges[expandEdgeIndex]->length
					< readLength + k)||(graph.degrees[i][1]==1)) {

				int a = 0;

				while (edgeRealId(graph.edgeIds[SourceVertex][a][OUT_EDGE],
						longEdges) != expandEdgeIndex){
					cerr<<"vertex "<<SourceVertex<<" out edge "<<edgeRealId(graph.edgeIds[SourceVertex][a][OUT_EDGE],
							longEdges)<<" on position "<<a<<endl;
					a++;
					assert(a < graph.degrees[SourceVertex][1]);

				}

				while (a < graph.degrees[SourceVertex][1] - 1) {
					graph.edgeIds[SourceVertex][a][OUT_EDGE]
							= graph.edgeIds[SourceVertex][a + 1][OUT_EDGE];
					a++;
					cerr<<"shift"<<endl;
				}
				graph.degrees[SourceVertex][1]--;

				forn(j,graph.degrees[i][1]) {
//					cerr<<"expand left "<<endl;
					longEdges[graph.edgeIds[i][j][OUT_EDGE]]->ExpandLeft(*(longEdges[expandEdgeIndex]));
					graph.edgeIds[SourceVertex][graph.degrees[SourceVertex][1]][OUT_EDGE]
							= graph.edgeIds[i][j][OUT_EDGE];
					graph.degrees[SourceVertex][1]++;
				}


				it = longEdges.find(expandEdgeIndex);
				longEdges.erase(it);
				graph.degrees[i][1] = 0;
				graph.degrees[i][0] = 0;
			}

		}

	}

	cerr << "expandDefinite finished\n";
}


inline bool equalsAtIndex(longEdgesMap &longEdges, int id1, int id2, int index, int direction) {
//	cerr << "<<";
//	cerr.flush();
	char u1, l1, u2, l2;
	int indu1, indu2, indl1, indl2;
	if (direction == RIGHT){
		indu1 = index;
		indu2 = index;
		indl1 = index;
		indl2 = index;
	//	cerr << index << " "<<id1  << longEdges[id1]->upper->size() << " " << longEdges[id1]->length;
	}
	else
	{
//		cerr << " hm";
		indu1 = longEdges[id1]->upper->size() - index - 1;
		indl1 = longEdges[id1]->lower->size() - index - 1;
		indu2 = longEdges[id2]->upper->size() - index - 1;
		indl2 = longEdges[id2]->lower->size() - index - 1;
//		cerr << indu1<<" "<<indl1<<" "<<indu2<<" "<<indl2 << " "<< index;
//		cerr.flush();
	}
	u1 = (*longEdges[id1]->upper)[indu1];
	l1 = (*longEdges[id1]->lower)[indl1];
	u2 = (*longEdges[id2]->upper)[indu2];
	l2 = (*longEdges[id2]->lower)[indl2];
//	cerr <<">>\n";
//	cerr.flush();
//	return (u1 == u2);
	return ((u1 == u2) && (l1 == l2));
}


inline int fistDifferenceIndex(PairedGraph &graph, Edge* edge1, Edge* edge2, int direction) {
	size_t res = 0;
	size_t maxLength = edge1->upper->size();
	if (maxLength > edge2->upper->size()) maxLength = edge2->upper->size();
	for(res=0; res < maxLength; res++){
		if (!equalsAtIndex(graph.longEdges, edge1->EdgeId, edge2->EdgeId, res, direction))
			break;
	}
	//res++;
	cerr<<" Edge "<< edge1->EdgeId<<" vs "<< edge2->EdgeId<<" for direction "<<direction<<" first diff on "<<res<<endl;
	if ((edge1->length<300)&&((edge2->length<300))){
		cerr<<edge1->upper->str()<<endl;
		cerr<<edge2->upper->str()<<endl;
		cerr<<edge1->lower->str()<<endl;
		cerr<<edge2->lower->str()<<endl;
	}
	return res;
}

void extractDefinite(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount, int dir)
{
//	LOGGER ("p.extractDefinite");
	int edgeIds[MAX_DEGREE];
	//TODO: constant?
	int delta = 20;
	int insertEdgeId = 10000;
	INFO("extractDefinite started");
	int maxV = VertexCount;
	forn(curVert, maxV){
	//	for (int direction = 0; direction < 2; direction ++)
		int direction = dir;
		{
//			INFO("Vertex number" << curVert);
			int index = k;
			int ended = -2;
			forn (edge, graph.degrees[curVert][direction]) {
				edgeIds[edge] = edgeRealId(
						graph.edgeIds[curVert][edge][direction], longEdges);
			}
			while (graph.degrees[curVert][direction] > 1) {
				forn(edge, graph.degrees[curVert][direction]){
					if (index >= longEdges[edgeIds[edge]]->length) {
						ended = edge;
						break;
					}

					//			cerr<< upp.str();
					//assert(0);
					if (!equalsAtIndex(longEdges, edgeIds[0], edgeIds[edge], index, direction))
					{
						ended = -1;
						break;
					}
				}
				if (ended >= -1)
					break;
				if (!(index & ((1 << 16) - 1)))
					INFO("index: "<<curVert << " "<< index);
				index++;
			}
			//index --;
			if (ended != -1) {
				// It seems to be impossible
				DEBUG("Intermediate Vertex " << curVert << "edgeId " << edgeIds[ended]);
			} else {
				if (index > k - 1 + delta) {
					DEBUG("Shortening edge");
					int newVert = VertexCount;
			//		DEBUG(index << " "<< longEdges[edgeIds[0]]->upper->str());
					Sequence *UpperSeq;
					Sequence *LowerSeq;
					int tl;
					if (direction) {
						UpperSeq = new Sequence(longEdges[edgeIds[0]]->upper->Subseq(0, index));
						LowerSeq = new Sequence(longEdges[edgeIds[0]]->lower->Subseq(0, index));
					}else{
						tl = longEdges[edgeIds[0]]->upper->size();
						cerr << "a";
						UpperSeq = new Sequence(longEdges[edgeIds[0]]->upper->Subseq(tl - index, tl));
						LowerSeq = new Sequence(longEdges[edgeIds[0]]->lower->Subseq(tl - index, tl));
						cerr.flush();
					}
					Edge* newEdge;
					if (direction)
						newEdge = new Edge(UpperSeq, LowerSeq, curVert, newVert, index - (k - 1), insertEdgeId);
					else
						newEdge = new Edge(UpperSeq, LowerSeq, newVert, curVert, tl - index - (k - 1), insertEdgeId);
					longEdges.insert(make_pair(insertEdgeId, newEdge));
					DEBUG ("new edge " << newEdge->upper->str());
					forn (tmpId, graph.degrees[curVert][direction]) {
						DEBUG ("before "<< longEdges[edgeIds[tmpId]]->upper->str());
						longEdges[edgeIds[tmpId]]->shortenEdge(index - (k - 1), direction);
						if (direction)
							longEdges[edgeIds[tmpId]]->FromVertex = newVert;
						else
							longEdges[edgeIds[tmpId]]->ToVertex = newVert;
						DEBUG ("after "<< longEdges[edgeIds[tmpId]]->upper->str());
						graph.edgeIds[newVert][tmpId][direction] = graph.edgeIds[curVert][tmpId][direction];
				//		DEBUG (longEdges[graph.edgeIds[newVert][tmpId][direction]]->upper->str());
					}
//					if (direction)
//						assert
					graph.degrees[newVert][direction]
							= graph.degrees[curVert][direction];
					graph.degrees[newVert][1 - direction] = 1;
					graph.edgeIds[newVert][0][1 - direction] = insertEdgeId;
					graph.degrees[curVert][direction] = 1;
					graph.edgeIds[newVert][0][direction] = insertEdgeId;
					insertEdgeId++;
					VertexCount++;

					assert(0);
				}
			}
		}
	}
}

void extractDefinite(PairedGraph &graph, int dir){
	cerr<<"extractDefinite create iterator "<<endl;

/*	int CurVertex = VIter->next();
	VertexIterator *VIter = graph.vertexIterator();
	cerr<<"extractDefinite create iterator OK"<<endl;
	VIter->hasNext();
	cerr<<"Has next OK"<<endl;
	for(;VIter->hasNext();){
	*/
	int CurVertex =0;
	while (CurVertex < graph.VertexCount){
		if (graph.degrees[CurVertex][0]+graph.degrees[CurVertex][1] > 0) break;
		CurVertex++;
	}

	while (CurVertex<graph.VertexCount){
		cerr<<"extractDefinite get next"<<endl;
		cerr<<"CurVertex "<<CurVertex<<endl;
		int FirstEdgeIndex = 0;
		int SecondEdgeIndex = 0;
		Edge* firstEdge;
		Edge* secondEdge;
		pair<Edge*,Edge*> firstEdgePair;
		pair<Edge*,Edge*> secondEdgePair;
		while (1){
			SecondEdgeIndex++;
			if (SecondEdgeIndex >= graph.degree(CurVertex, dir)) {
				FirstEdgeIndex++;
				SecondEdgeIndex = FirstEdgeIndex+1;
			}
			cerr<<"firstIndex "<<FirstEdgeIndex<<" secondIndex "<<SecondEdgeIndex<<" degree "<<graph.degree(CurVertex, dir)<<endl;
			if (FirstEdgeIndex+1 >= graph.degree(CurVertex, dir)) break;
			firstEdge = graph.neighbourEdge(CurVertex, FirstEdgeIndex, dir);
			cerr<<"first edge Id "<<firstEdge->EdgeId<<endl;
			secondEdge = graph.neighbourEdge(CurVertex, SecondEdgeIndex, dir);
			cerr<<"second edge Id "<<secondEdge->EdgeId<<endl;
			int diffIndex = fistDifferenceIndex(graph, firstEdge, secondEdge, dir);
			if (diffIndex <= k-1) continue;
			int splitIndex = diffIndex - k+1;
			cerr<<"splitIndex "<<splitIndex<<endl;
//			if (dir == RIGHT){
				if (diffIndex>firstEdge->length){

					if (diffIndex>secondEdge->length) continue;
					else {
						secondEdgePair = graph.splitEdge(secondEdge, splitIndex, dir);
						graph.glueEdges(secondEdgePair.first, firstEdge);
						SecondEdgeIndex = FirstEdgeIndex;
					}
				}
				else {
					firstEdgePair = graph.splitEdge(firstEdge, splitIndex, dir);
					if (diffIndex>secondEdge->length) {
						graph.glueEdges(firstEdgePair.first, secondEdge);
						SecondEdgeIndex = FirstEdgeIndex;
					}
					else {
						secondEdgePair = graph.splitEdge(secondEdge, splitIndex, dir);
						graph.glueEdges(firstEdgePair.first, secondEdgePair.first);
						SecondEdgeIndex = FirstEdgeIndex;
					}
				}
	//		}
//			else {
//				INFO("extractDefinite for LEFT not implemented yet");
//				assert(0);
//			}
		}
		CurVertex++;
		while (CurVertex < graph.VertexCount){
			if (graph.degrees[CurVertex][0]+graph.degrees[CurVertex][1] > 0) break;
			CurVertex++;
		}

	}
}



void SplitVertecesByEdgeConnections(PairedGraph &graph, edgePairsMap &EdgePairs, bool Strongly1to1){
	//resolve multi case;
	int table[MAX_DEGREE][MAX_DEGREE];
	double fake_cov[MAX_DEGREE][MAX_DEGREE];
	double in_cov[MAX_DEGREE];
	double out_cov[MAX_DEGREE];
	double eps_out_cov[MAX_DEGREE];
	double eps = 10000;


	int FakeVertexCount = graph.VertexCount;
	int FakeVertexStart = graph.VertexCount;
	map<int, int> FakeVertexToReal;
	forn(curVertId,FakeVertexStart) {
		if ((graph.degree(curVertId, LEFT)!=0)&&(graph.degree(curVertId, RIGHT)!=0)) {

			pair<int, int> vDist = vertexDist(graph.longEdges,graph,curVertId);
			//			if (vDist.first+vDist.second-k+1<=readLength)
			{
				cerr<<"vertex "<<curVertId<<" dist <"<<vDist.first<<", "<<vDist.second<<">"<<endl;
				if (graph.degree(curVertId, LEFT)*graph.degree(curVertId, RIGHT) <= EdgePairs[curVertId].size()) continue;

				if (((graph.degree(curVertId, LEFT)<=EdgePairs[curVertId].size())&&(graph.degrees[curVertId][1]<=EdgePairs[curVertId].size()))||(!Strongly1to1)) {
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

					if ((allIns&&allOuts)||(!Strongly1to1)) {

						TRACE ("Faking vert " << curVertId);
						cerr<<"Faking vert " << curVertId<<endl;
						int tmpCurOut = edgeRealId(graph.neighbourEdge(curVertId,0,RIGHT)->EdgeId,graph.longEdges);
						string tmpUpSeq = graph.longEdges[tmpCurOut]->upper->Subseq(0,k-1).str();
//						string tmpLoSeq = graph.longEdges[tmpCurOut]->lower->Subseq(0,l-1).str();
						string tmpLoSeq = "";

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
							int CurIn = edgeRealId(graph.edgeIds[curVertId][i][IN_EDGE],graph.longEdges);
							graph.edgeIds[curVertId][i][IN_EDGE] = CurIn;
							FakeVertexToReal.insert(make_pair(FakeVertexCount,curVertId));

							cerr<<"Edge "<<CurIn<<" ("<<graph.longEdges[CurIn]->FromVertex<<", "<<graph.longEdges[CurIn]->ToVertex<<") maped to ";
							graph.longEdges[CurIn]->ToVertex = FakeVertexCount;
							graph.degrees[FakeVertexCount][0]=1;
							graph.degrees[FakeVertexCount][1]=0;
							in_cov[i] = graph.longEdges[CurIn]->coverage;
							graph.edgeIds[FakeVertexCount][0][IN_EDGE]=CurIn;
							cerr<<"edge "<<CurIn<<" ("<<graph.longEdges[CurIn]->FromVertex<<", "<<graph.longEdges[CurIn]->ToVertex<<")"<<endl;
							FakeVertexCount++;
							graph.VertexCount++;
						}



						//create fake Vertices for out edges;
						int tmpFictStartOut = FakeVertexCount;
						memset(out_cov, 0, sizeof(out_cov));
						forn(i,graph.degrees[curVertId][1]) {
							int CurOut = edgeRealId(graph.edgeIds[curVertId][i][OUT_EDGE],graph.longEdges);
							graph.edgeIds[curVertId][i][OUT_EDGE] = CurOut;
							cerr<<"Edge "<<CurOut<<" ("<<graph.longEdges[CurOut]->FromVertex<<", "<<graph.longEdges[CurOut]->ToVertex<<") maped to ";
							FakeVertexToReal.insert(make_pair(FakeVertexCount,curVertId));
							graph.longEdges[CurOut]->FromVertex = FakeVertexCount;
							graph.degrees[FakeVertexCount][0]=0;
							graph.degrees[FakeVertexCount][1]=1;
							graph.edgeIds[FakeVertexCount][0][OUT_EDGE]=CurOut;

							out_cov[i] = graph.longEdges[CurOut]->coverage;
							cerr<<"edge "<<CurOut<<" ("<<graph.longEdges[CurOut]->FromVertex<<", "<<graph.longEdges[CurOut]->ToVertex<<")"<<endl;
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
							while (edgeRealId(graph.edgeIds[curVertId][tmpFrom][IN_EDGE], graph.longEdges)!=edgeRealId((EdgePairs[curVertId])[tmpEdgePair].first,graph.longEdges)) {
								tmpFrom++;
								assert(tmpFrom<graph.degrees[curVertId][0]);
							}
							while (edgeRealId(graph.edgeIds[curVertId][tmpTo][OUT_EDGE], graph.longEdges)!=edgeRealId((EdgePairs[curVertId])[tmpEdgePair].second, graph.longEdges)) {
								tmpTo++;
								assert(tmpTo<graph.degrees[curVertId][1]);
							}
							Sequence *UpSeq = new Sequence(tmpUpSeq);
							Sequence *LoSeq = new Sequence(tmpLoSeq);

//							Edge *tmpEdge = new Edge(UpSeq, LoSeq, tmpFictStartIn + tmpFrom, tmpFictStartOut + tmpTo,0, EdgeId, fake_cov[edgeIdToLocalId(0, graph, EdgePairs[curVertId][tmpEdgePair].first)][edgeIdToLocalId(1, graph, EdgePairs[curVertId][tmpEdgePair].second)]);
							Edge *tmpEdge = new Edge(UpSeq, LoSeq, tmpFictStartIn + tmpFrom, tmpFictStartOut + tmpTo,0, 0, fake_cov[tmpFrom][tmpTo]);
							int tmpEdgeId = graph.addEdge(tmpEdge,true)->EdgeId;
							//							longEdges.insert(make_pair(EdgeId,tmpEdge));
							//							cerr<<"Virtual edge "<<EdgeId<<
//							cerr<<"Virtual edge "<<graph.EdgeId-1<<" ("<<graph.longEdges[EdgeId]->FromVertex<<", "<<graph.longEdges[EdgeId]->ToVertex<<") ["<<tmpFrom<<", "<<tmpTo<<"] cov "<< fake_cov[tmpFrom][tmpTo]<<endl;
					//		graph.edgeIds[tmpFictStartIn + tmpFrom][graph.degrees[tmpFictStartIn + tmpFrom][1]][OUT_EDGE] = tmpEdgeId;
					//		graph.degrees[tmpFictStartIn + tmpFrom][1]++;
					//		graph.edgeIds[tmpFictStartOut + tmpTo][graph.degrees[tmpFictStartOut + tmpTo][0]][OUT_EDGE] = tmpEdgeId;
					//		graph.degrees[tmpFictStartOut + tmpTo][0]++;
						//	graph.EdgeId++;

						}

						graph.degrees[curVertId][0] = 0;
						graph.degrees[curVertId][1] = 0;
					}
				}
			}
		}
	}

	graph.VertexCount = FakeVertexCount;
	graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
	if (Strongly1to1){
		expandDefinite(graph.longEdges, graph, graph.VertexCount, false);
	}else {
		cerr<<"Start expand obvious"<<endl;
		expandObvious(graph);
		cerr<<"End expand obvious"<<endl;
		cerr<<"Clear fake edges"<<endl;

		for(longEdgesMap::iterator it= graph.longEdges.begin(); it !=graph.longEdges.end();) {
			if (it->second->length == 0) {
				graph.glueVertices(it->second->FromVertex, it->second->ToVertex);
				graph.removeEdge((it++)->second);
//				++it;
			}
			else ++it;
		}
		cerr<<"Clear fake edges END"<<endl;

	}
	//	expandDefinite(longEdges, graph, VertexCount);
//	for(longEdgesMap::iterator it= longEdges.begin(); it !=longEdges.end(); ++it) {
//		if (it->second->FromVertex>=FakeVertexStart) it->second->FromVertex = FakeVertexToReal[it->second->FromVertex];
//		if (it->second->ToVertex>=FakeVertexStart) it->second->ToVertex = FakeVertexToReal[it->second->ToVertex];
//	}
//	VertexCount = FakeVertexStart;
//	graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
}

bool intersectible(Sequence *left, Sequence *right){
	string leftStr;
	string rightStr;
//	cerr<<"Left "<<left->size()<<" VS Right"<<right->size()<<endl;
	if (left->size()<=1) return true;
	if (right->size()<=1) return true;
	if (left->size()>350)
		leftStr = left->Subseq(left->size()-350, left->size()).str();
	else
		leftStr = left->str();
	if (right->size()>350)
		rightStr = right->Subseq(0, 350).str();
	else
		rightStr = right->str();
	pair<int, pair<int,int> > tmp = maxCommonSubstring(leftStr, rightStr);
	int border = min(left->size()/2, right->size()/2);
	if (border > 175) border = 175;
	if (tmp.first > max(l, border)) return true;
	else return false;

}
void dfs (int **table, int color, int * leftcolor, int* rightcolor,  int pos) {
	leftcolor[pos] = color;
	forn(j, MAX_DEGREE) {
		if (table[pos][j] && !(rightcolor[j])) {
			rightcolor[j] = color;
			forn(i, MAX_DEGREE)
				if (table[i][j] && i != pos && !leftcolor[i]) {
					leftcolor[i] = color;
					dfs(table,color,leftcolor,rightcolor,i);
				}
		}
	}
}
void doSplit(PairedGraph &graph, edgePairsMap &EdgePairs) {
	int *table[MAX_DEGREE];
	forn(i, MAX_DEGREE)
		table[i] = new int[MAX_DEGREE];
	int leftcolor [MAX_DEGREE];
	int rightcolor [MAX_DEGREE];
	map<int, int> edgeIds;
	map<int, int> idEdges;
	for(edgePairsMap::iterator iter = EdgePairs.begin(); iter != EdgePairs.end(); iter++) {
		int curVId = iter->first;
		INFO(" splitting vertex" << curVId);
		int len = iter->second.size();
		forn(i, MAX_DEGREE) {
			forn(j, MAX_DEGREE)
				table[i][j] = 0;
			leftcolor[i] = 0;
			rightcolor[i] = 0;
		}
		cerr<<"Vertex "<<curVId<<" in degree "<<graph.degree(curVId,LEFT)<<" out degree "<<graph.degree(curVId,RIGHT)<<" edge pairs:"<<endl;
		forn(i, len){
			cerr<<iter->second[i].first<<" "<<iter->second[i].second<<endl;
		}

		edgeIds.clear();
		idEdges.clear();
		int leftId = 1;
		int rightId = 1;
		forn(i, len) {
			if (edgeIds.find(iter->second[i].first) == edgeIds.end()){
				edgeIds[iter->second[i].first] = leftId;
				idEdges[leftId] = iter->second[i].first;
				leftId ++;
			}
			if (edgeIds.find(-iter->second[i].second) == edgeIds.end()){
				edgeIds[-iter->second[i].second] = rightId;
				idEdges[-rightId] = iter->second[i].second;
				rightId ++;
			}
		//	leftGlobalIds[leftId] = mp(iter->second[i].first,
			table[edgeIds[iter->second[i].first]][edgeIds[-iter->second[i].second]] = 1;
		}
		int in_degree = graph.degree(curVId, LEFT);
		int out_degree = graph.degree(curVId, RIGHT);
		forn(i, in_degree) {
			if (edgeIds.find(graph.leftEdge(curVId, i)->EdgeId) == edgeIds.end()){
				edgeIds[graph.leftEdge(curVId, i)->EdgeId] = leftId;
				idEdges[leftId] = graph.leftEdge(curVId, i)->EdgeId;
				leftId++;

//				EdgePairs.push_back(mp(graph.leftEdge(curVId, i)->EdgeId,0));
			}
		}
		forn(i, out_degree) {
			if (edgeIds.find(-graph.rightEdge(curVId, i)->EdgeId) == edgeIds.end()){
				edgeIds[-graph.rightEdge(curVId, i)->EdgeId] = rightId;
				idEdges[-rightId] = graph.rightEdge(curVId, i)->EdgeId;
				rightId++;
			}
		}
		int color = 0;
		for(int i = 1; i < leftId; i++) {
			if(!leftcolor[i]) {

				color++;
				dfs(table, color, leftcolor, rightcolor, i);
			}
		}
		for(int i =1; i < rightId; i++) {
			if(!rightcolor[i]) {
				color++;
				rightcolor[i] = color;
			}
		}
		cerr<<"left colors: "<<endl;
		forn(i, leftId) {
			cerr<<leftcolor[i]<<" ";
		}
		cerr<<"\n right colors: "<<endl;
		forn(i, rightId) {
			cerr<<rightcolor[i]<<" ";
		}
		cerr<<"table: "<<endl;
		forn(i, leftId) {
			forn(j, rightId) {
			cerr<<table[i][j]<<" ";
			}

			cerr<<endl;
		}



		int tmpVertCount = graph.VertexCount;
		cerr << "idEdges" << endl;
		for(map<int, int>::iterator it = idEdges.begin(); it != idEdges.end(); it++)
			cerr<< it->first << " " << it->second << endl;
		for(int cur_color = 2;cur_color <= color; cur_color++) {

			int tmp = graph.addVertex(graph.VertexCount + cur_color - 2);
			INFO("added vertex" << tmp);
		}
		forn(i, leftId) {
			if (i && leftcolor[i] >= 2)
				graph.addEdgeVertexAdjacency(tmpVertCount + leftcolor[i] - 2, graph.longEdges[idEdges[i]], LEFT);
		}
		forn(i, rightId) {
			if (i && rightcolor[i] >= 2)
				graph.addEdgeVertexAdjacency(tmpVertCount + rightcolor[i] - 2, graph.longEdges[idEdges[-i]], RIGHT);
		}
		//graph.VertexCount += color - 2;
	}
	graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);

}

void SplitByLowers(PairedGraph &graph){
	edgePairsMap EdgePairs;
	forn(CurVert, graph.VertexCount){
		vector<pair<int,int>> tmpVect;
		forn (edgeI, graph.degree(CurVert, LEFT)){
			Edge* leftEdge = graph.neighbourEdge(CurVert,edgeI,LEFT);
			forn (edgeJ, graph.degree(CurVert, RIGHT)){
				Edge* rightEdge = graph.neighbourEdge(CurVert,edgeJ,RIGHT);
	//			cerr<<"Check: "<<CurVert<<" Edge pair "<<leftEdge->EdgeId<<" "<<rightEdge->EdgeId<<endl;

				if (intersectible(leftEdge->lower, rightEdge->lower)){
					tmpVect.push_back(make_pair(leftEdge->EdgeId, rightEdge->EdgeId));
					cerr<<"Vert "<<CurVert<<" Edge pair "<<leftEdge->EdgeId<<" "<<rightEdge->EdgeId<<endl;
				}
			}
		}
		if (tmpVect.size() > 1)
			EdgePairs.insert(make_pair(CurVert,tmpVect));
	}
	cerr<<"Start Spliting"<<endl;
//	SplitVertecesByEdgeConnections(graph, EdgePairs, false);
	doSplit(graph, EdgePairs);
	cerr<<"End Spliting"<<endl;
}

