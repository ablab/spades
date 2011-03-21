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


bool processLowerSequence(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount){
	int countIn[MAX_DEGREE], possibleIn[MAX_DEGREE], countOut[MAX_DEGREE], possibleOut[MAX_DEGREE];
	bool res = false;
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
/* countEdgesDistinction
 * @param vertexId - index of vertex, distances we are interested in
 * @param graph - pairedGraph
 * @returns pair<int,int> - minimal intervals of nucleotids, such that they uniquely represent Inedges and outedges respectively.
 * If one of results  is more than insertLength + 2 * readLength, returns insertLength + 2 * readLength + 1
 *
 */
pair<int, int> vertexDist(longEdgesMap &longEdges, PairedGraph &graph, int vertexId){
	int res1 = k-1;
	int res2 = k-1;
	int max_res = insertLength + 2*readLength + 1;
	//int colors[MAX_DEGREE];
	forn(i, graph.inD[vertexId]) {
		int e1 = graph.inputEdges[vertexId][i];
		Sequence tmp = *longEdges[e1]->upper;
		Sequence tmpl = *longEdges[e1]->lower;

//		cerr << endl << "1 " << tmp.str();
		forn(j, i) {
			int count = 0;

			int e2 = graph.inputEdges[vertexId][j];
			while (count < max_res) {
				Sequence tmp2 = *longEdges[e2]->upper;
				Sequence tmp2l = *longEdges[e2]->lower;

//				cerr << endl <<"2 "<< tmp2.str();
//				cerr << endl <<count << " ";// << nucl(tmp[tmp.size() - 1 - count]) << " " << nucl (tmp2[tmp2.size() - 1 -count]);
				if (count >= min(tmp.size(), tmp2.size()) || tmp[tmp.size() - 1 - count] != tmp2[tmp2.size() - 1 -count] )
				//				if (*(longEdges[e1]->upper)[count] != *(longEdges[e2]->upper)[count])
					break;
				if (count >= min(tmpl.size(), tmp2l.size()) || tmpl[tmpl.size() - 1 - count] != tmp2l[tmp2l.size() - 1 -count] )
				//				if (*(longEdges[e1]->upper)[count] != *(longEdges[e2]->upper)[count])
					break;

				count ++;
			}
			res1 = max(res1, count);
		}
	}
	forn(i, graph.outD[vertexId]) {
		int e1 = graph.outputEdges[vertexId][i];
		Sequence tmp = *longEdges[e1]->upper;
		Sequence tmpl = *longEdges[e1]->lower;

//		cerr << endl << "1 " << tmp.str();
		forn(j, i) {
			int e2 = graph.outputEdges[vertexId][j];
			int count = 0;
			while (count < max_res) {
				Sequence tmp2 = *longEdges[e2]->upper;
				Sequence tmp2l = *longEdges[e2]->lower;

//				cerr << endl <<"2 "<< tmp2.str();
//				cerr << endl <<count << " ";// << nucl(tmp[tmp.size() - 1 - count]) << " " << nucl (tmp2[tmp2.size() - 1 -count]);
				if (count >= min(tmp.size(), tmp2.size()) || tmp[tmp.size() - 1 - count] != tmp2[tmp2.size() - 1 -count] )
				//				if (*(longEdges[e1]->upper)[count] != *(longEdges[e2]->upper)[count])
					break;
				if (count >= min(tmpl.size(), tmp2l.size()) || tmpl[tmpl.size() - 1 - count] != tmp2l[tmp2l.size() - 1 -count] )
				//				if (*(longEdges[e1]->upper)[count] != *(longEdges[e2]->upper)[count])
					break;

				count ++;
			}
			res2 = max(res2, count);
		}
	}

	cerr << endl << res1 << " " << res2<<endl;
	//assert(0);
	return make_pair(res1 + 1, res2 + 1);
}



