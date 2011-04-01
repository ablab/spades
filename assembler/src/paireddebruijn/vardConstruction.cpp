#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"
#include "seq.hpp"
#include "graphVisualizer.hpp"
#include "pairedGraph.hpp"
#include "graphio.hpp"

LOGGER("p.vardConstruction");

char EdgeStr[1000000];
char EdgeStrLo[1000000];

using namespace paired_assembler;
namespace vard {
void createVertices(edgesMap &edges, PairedGraph &graph) {
	int count = 0;
	for (edgesMap::iterator iter = edges.begin(); iter != edges.end();) {
		int size = iter->second.size();
		ll kmer = iter->fi;
		forn(i, size) {
			if ((!(iter->se)[i]->used)) {
				int length = 1;
				int EdgeCoverage;
				count++;
				EdgePrototype* curEdgePrototype = (iter->se)[i];
				curEdgePrototype->used = true;

//				int leftIndex = findPossibleVertex(leftkmer(kmer), curEdgePrototype->lower, graph);
//				assert(leftIndex != -2);
//				expandRight()
			}
		}
	}
}
inline ll leftkmer(ll kmer) {
	return kmer >> 2;
}
inline ll rightkmer(ll kmer) {
	return kmer & (~((ll)3 << (2*k-2)));
}
inline ll subkmer(ll kmer, int direction) {
	if (direction == LEFT)
		return leftkmer(kmer);
	else if (direction == RIGHT)
		return rightkmer(kmer);
	else assert(0);
}

/*
 *
 * @return index of vertex, if this k-1 -seq can be traced to it.
 * If there are multiple, return -2,(we hope there will be no such situation:)
 * if no- returns -1
 */

int findPossibleVertex(ll &kmer, Sequence &down, PairedGraph &graph){
	return 0;
}

int expandDirected(edgesMap &edges, verticesMap &verts, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage, int direction){
	assert(direction == LEFT || direction == RIGHT );

}
int checkUniqueWay(edgesMap &edges, ll finishKmer, Sequence *finishSeq , int direction){
	assert(direction == LEFT || direction == RIGHT );
}
int goUniqueWay(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage, int direction){
	assert(direction == LEFT || direction == RIGHT );
}

}
