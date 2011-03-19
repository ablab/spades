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
		inD[i]=0;
		outD[i]=0;
	}
	for(longEdgesMap::iterator it= longEdges.begin(); it !=longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first){
			outputEdges[it->second->FromVertex][outD[it->second->FromVertex]++]=it->first;
			inputEdges[it->second->ToVertex][inD[it->second->ToVertex]++]=it->first;
		}
	}
}
}
