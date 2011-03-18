/*
 * graphSimplification.cpp
 *
 *  Created on: 17.03.2011
 */
#include "graphSimplification.hpp"
#include "common.hpp"

bool isPath(Edge &e1, Edge &e2) {
	cerr<<endl<<"("<<e1.FromVertex<<", "<<e1.ToVertex<<") "<<"("<<e2.FromVertex<<", "<<e2.ToVertex<<") ";
	if (e1.ToVertex != e2.FromVertex)
		return false;
	int sts = readLength + insertLength;
	if (e1.length + e2.length + k - 1 <= sts)
		return true;
	int left1 = max(0, e1.length - sts);
	int right1 = min(e1.upper->size(), e1.upper->size() - sts + e2.length);
	int left2 = left1 + sts - e1.length;
	int right2 = right1 + sts - e1.length;
	cerr<<endl<<e1.lower->Subseq(left1, right1).str()<<endl<<e2.upper->Subseq(left2, right2).str();
	return e1.lower->Subseq(left1, right1) == e2.upper->Subseq(left2, right2);
}


void processLowerSequence(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount){
	int countIn[MAX_DEGREE], possibleIn[MAX_DEGREE], countOut[MAX_DEGREE], possibleOut[MAX_DEGREE];
	forn(curVertId, VertexCount){
		if ((graph.inD[curVertId]!=0)&&(graph.outD[curVertId]!=0)) {

			forn(i,MAX_DEGREE){
				countOut [i] =0;
				countIn [i] =0;
			}
			forn(i,graph.inD[curVertId]){
				int curInEdgeId = edgeRealId(graph.inputEdges[curVertId][i], longEdges);
				graph.inputEdges[curVertId][i] = curInEdgeId;
				forn(j,graph.outD[curVertId]){
					int curOutEdgeId = edgeRealId(graph.outputEdges[curVertId][j], longEdges);
					graph.outputEdges[curVertId][j] = curOutEdgeId;
					cerr<<"Check isPass for edge "<<curInEdgeId<<" vs "<< curOutEdgeId;
					if (isPath(*longEdges[curInEdgeId], *longEdges[curOutEdgeId]))
					{
						countIn[j]++;
						countOut[i]++;
						possibleIn[j] = i;
						possibleOut[i] =j;
						cerr<<" POSSIBLE"<<endl;
					}
					else cerr<<" IMPOSSIBLE"<<endl;
				}
			}
			forn(i,graph.inD[curVertId]){
				if (countOut[i]==1){
					if (countIn[possibleOut[i]]==1){
						int InEdge = edgeRealId(graph.inputEdges[curVertId][i], longEdges);
						int OutEdge = edgeRealId(graph.outputEdges[curVertId][possibleOut[i]], longEdges);
						longEdges[InEdge]->ExpandRight(*longEdges[OutEdge]);
						longEdges[OutEdge] = longEdges[InEdge];
					}
				}
			}
			//REMOVE VERTICES INFO!!!
		}
	}
}

