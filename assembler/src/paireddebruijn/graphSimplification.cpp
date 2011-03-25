/*
 * graphSimplification.cpp
 *
 *  Created on: 17.03.2011
 */
#include "graphSimplification.hpp"
#include "common.hpp"
LOGGER ("p.graphSimplification");
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
					cerr << "Check isPass for edge " << curInEdgeId << " vs "
							<< curOutEdgeId;
					if (isPath(*longEdges[curInEdgeId],
							*longEdges[curOutEdgeId])) {
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
	int res1 = k - 1;
	int res2 = k - 1;
	int max_res = insertLength + 2 * readLength + 1;
	//int colors[MAX_DEGREE];
	forn(i, graph.degrees[vertexId][0]) {
		int e1 = graph.edgeIds[vertexId][i][IN_EDGE];
		Sequence tmp = *longEdges[e1]->upper;
		Sequence tmpl = *longEdges[e1]->lower;

		//		cerr << endl << "1 " << tmp.str();
		forn(j, i) {
			int count = 0;

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
			int count = 0;
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

void expandDefinite(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount, bool NotExpandBeyondDefinite) {
	longEdgesMap::iterator it;
	int expandEdgeIndex;
	cerr << "expandDefiniteStart" << endl;
	forn(i,VertexCount) {
		if ((graph.degrees[i][1] == 1) && (graph.degrees[i][0] > 0)) {
			expandEdgeIndex = edgeRealId(graph.edgeIds[i][0][OUT_EDGE],
					longEdges);
			int DestVertex = longEdges[expandEdgeIndex]->ToVertex;
			pair<int, int> diffDistDest = make_pair(0, 0);
			pair<int, int> diffDistCur = make_pair(0, 0);
			if (NotExpandBeyondDefinite) {
				diffDistDest = vertexDist(longEdges, graph, DestVertex);
				diffDistCur = vertexDist(longEdges, graph, i);
			}
			if ((!NotExpandBeyondDefinite) || (diffDistCur.first
					+ diffDistDest.second + longEdges[expandEdgeIndex]->length
					< readLength + k)) {
				if (NotExpandBeyondDefinite)
					cerr << "Check cur vert " << i << " dest vert "
							<< DestVertex << "  " << diffDistCur.first << " + "
							<< diffDistDest.second << " + "
							<< longEdges[expandEdgeIndex]->length << " = "
							<< diffDistCur.first + diffDistDest.second
									+ longEdges[expandEdgeIndex]->length
							<< " < " << readLength + k << endl;
				int a = 0;
				//				cerr << "trying to expand";
				while ((edgeRealId(graph.edgeIds[DestVertex][a][IN_EDGE],
						longEdges) != expandEdgeIndex))
					a++;
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
		if ((graph.degrees[i][0] == 1) && (graph.degrees[i][1] > 0)) {
			cerr << i << endl;
			expandEdgeIndex = edgeRealId(graph.edgeIds[i][0][IN_EDGE],
					longEdges);
			int SourceVertex = longEdges[expandEdgeIndex]->FromVertex;
			pair<int, int> diffDistSource = make_pair(0, 0);
			pair<int, int> diffDistCur = make_pair(0, 0);
			if (NotExpandBeyondDefinite) {
				diffDistSource = vertexDist(longEdges, graph, SourceVertex);
				diffDistCur = vertexDist(longEdges, graph, i);
			}
			if ((!NotExpandBeyondDefinite) || (diffDistCur.second
					+ diffDistSource.first + longEdges[expandEdgeIndex]->length
					< readLength + k)) {

				int a = 0;
				while (edgeRealId(graph.edgeIds[SourceVertex][a][OUT_EDGE],
						longEdges) != expandEdgeIndex)
					a++;
				assert(a < graph.degrees[SourceVertex][1]);
				while (a < graph.degrees[SourceVertex][1] - 1) {
					graph.edgeIds[SourceVertex][a][OUT_EDGE]
							= graph.edgeIds[SourceVertex][a + 1][OUT_EDGE];
					a++;
				}
				graph.degrees[SourceVertex][1]--;
				forn(j,graph.degrees[i][1]) {
					longEdges[graph.edgeIds[i][j][OUT_EDGE]]->ExpandLeft(
							*(longEdges[expandEdgeIndex]));
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
	cerr << "<<";
	cerr.flush();
	char u1, l1, u2, l2;
	int indu1, indu2, indl1, indl2;
	if (direction){
		indu1 = index;
		indu2 = index;
		indl1 = index;
		indl2 = index;
		cerr << index << " "<<id1  << longEdges[id1]->upper->size() << " " << longEdges[id1]->length;
	}
	else
	{
		cerr << " hm";
		indu1 = longEdges[id1]->upper->size() - index - 1;
		indl1 = longEdges[id1]->lower->size() - index - 1;
		indu2 = longEdges[id2]->upper->size() - index - 1;
		indl2 = longEdges[id2]->lower->size() - index - 1;
		cerr << indu1<<" "<<indl1<<" "<<indu2<<" "<<indl2 << " "<< index;
		cerr.flush();
	}
	u1 = (*longEdges[id1]->upper)[indu1];
	l1 = (*longEdges[id1]->lower)[indl1];
	u2 = (*longEdges[id2]->upper)[indu2];
	l2 = (*longEdges[id2]->lower)[indl2];
	cerr <<">>\n";
	cerr.flush();
	return (u1 == u2 && l1 == l2);
}
void extractDefinite(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount, int dir)
{
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

					forn (tmpId, graph.degrees[curVert][direction]) {
						longEdges[edgeIds[tmpId]]->shortenEdge(index - (k - 1), direction);
						if (direction)
							longEdges[edgeIds[tmpId]]->FromVertex = newVert;
						else
							longEdges[edgeIds[tmpId]]->ToVertex = newVert;
						graph.edgeIds[newVert][tmpId][direction] = graph.edgeIds[curVert][tmpId][direction];
				//		DEBUG (longEdges[graph.edgeIds[newVert][tmpId][direction]]->upper->str());
					}
					graph.degrees[newVert][direction]
							= graph.degrees[curVert][direction];
					graph.degrees[newVert][1 - direction] = 1;
					graph.edgeIds[newVert][0][1 - direction] = insertEdgeId;
					graph.degrees[curVert][direction] = 1;
					graph.edgeIds[newVert][0][direction] = insertEdgeId;
					insertEdgeId++;
					VertexCount++;
					//		assert(0);
				}
			}
		}
	}
}
