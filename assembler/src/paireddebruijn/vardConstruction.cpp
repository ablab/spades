#include "constructHashTable.hpp"
#include "common.hpp"
#include "graphConstruction.hpp"
#include "seq.hpp"
#include "graphVisualizer.hpp"
#include "pairedGraph.hpp"
#include "graphio.hpp"
#include "vardConstruction.hpp"

LOGGER("p.vardConstruction");

using namespace paired_assembler;
namespace vard {
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
 * @param direction  LEFT if we look from leftmost end
 *
 * @return index of vertex, if this k-1 -seq can be traced to it.
 * If there are multiple, return -2,(we hope there will be no such situation:)
 * if no- returns -1
 */

int findPossibleVertex(ll &kmer, Sequence &down, edgesMap &edges, verticesMap &verts){
	return 0;
}
//go left until vertex or fork.
// Go right until vertex or fork, adding nucleotides to curEdge strings.
int expandDirected(edgesMap &edges, curEdgeType &curEdge, verticesMap &verts, ll &curKmer, Sequence* &curSeq, int &EdgeCoverage, int direction){
	assert(direction == LEFT || direction == RIGHT );
	int otherDirection;
	if (direction == LEFT)
		otherDirection = RIGHT;
	else
		otherDirection = LEFT;
	while(!checkUniqueWay(edges, curKmer, curSeq, otherDirection)) {
		if ((findPossibleVertex(curKmer, *curSeq, edges, verts) >= 0) || (!goUniqueWay(edges, curKmer, curSeq, EdgeCoverage, direction))) {
			break;
		}
		if (direction == RIGHT) {
			//TODO: Р’РЎРўРђР’Р�РўР¬ РѕС‚РјР°С‚С‡РµРЅРЅРѕРµ.
		}
	}
}
int checkUniqueWay(edgesMap &edges, ll finishKmer, Sequence *finishSeq , int direction){
	assert(direction == LEFT || direction == RIGHT );
}

//while going left we don't mark anything as used, we just find leftmost possible vert
int goUniqueWay(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq, int &EdgeCoverage, int direction){
	assert(direction == LEFT || direction == RIGHT );
}


void createVertices(edgesMap &edges, PairedGraph &graph) {
	int count = 0;
	for (edgesMap::iterator iter = edges.begin(); iter != edges.end();) {
		int size = iter->second.size();
		ll kmer = iter->fi;
		forn(i, size) {
			if ((!(iter->se)[i]->used)) {
				int length = 0;
				int EdgeCoverage;
				count++;
				curEdgeType curEdge;
				EdgePrototype* curEdgePrototype = (iter->se)[i];
				curEdgePrototype->used = true;
				Sequence * startSeq = curEdgePrototype->lower;
				int curshift = 0;
				ll startKmer = subkmer(kmer, LEFT);


				expandDirected(edges, curEdge, graph.verts, startKmer, startSeq, EdgeCoverage, LEFT);

				curEdge.first = "";
				curEdge.second = "";
				int startVertId = findPossibleVertex(startKmer, *startSeq, edges, graph.verts);

				assert(startVertId != -2);
				ll finKmer = startKmer;
				Sequence *finSeq = new Sequence(&startSeq);
				expandDirected(edges, curEdge, graph.verts, finKmer,finSeq, EdgeCoverage, LEFT);
				while(!checkUniqueWay(edges, finKmer, finSeq, LEFT)) {
					if (! goUniqueWay(edges, finKmer, finSeq, EdgeCoverage, RIGHT));
						break;
				}
				int finVertId = findPossibleVertex(startKmer, *startSeq, edges, graph.verts);
				assert(finVertId != -2);

				if (startVertId < 0  && finVertId < 0) {
					//TODO: in  fact, not 0 but ...
					startVertId = storeVertex(graph, startKmer, startSeq, 0);
					finVertId = storeVertex(graph, finKmer, finSeq, 0);
				}
				//expandDirected(edges, curEdge, graph.verts, startKmer, startSeq, EdgeCoverage, LEFT);
			}
		}
	}
}



}
