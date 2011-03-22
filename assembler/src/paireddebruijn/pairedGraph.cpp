#include "pairedGraph.hpp"
//
//int isDeleted(Node* A) {
//	return A->upperSize == 0;
//}
//int deleteNode(Node* A) {
//	return A->upperSize = 0;
//}
namespace paired_assembler {
void PairedGraph::recreateVerticesInfo(int vertCount, longEdgesMap &longEdges)
{
	forn(i, vertCount) {
		forn(j, 2)
		degrees[i][j]=0;
	}
	for(longEdgesMap::iterator it= longEdges.begin(); it !=longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first){
			outputEdges[it->second->FromVertex][degrees[it->second->FromVertex][1]++]=it->first;
			inputEdges[it->second->ToVertex][degrees[it->second->ToVertex][0]++]=it->first;
		}
	}
}
}
