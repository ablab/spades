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
/*\
 *\
 *
 * @return coverage of resulting edge when expanding or 0.
 */
int expandDirected(edgesMap &edges, protoEdgeType &curEdge, verticesMap &verts, ll &curKmer, Sequence* &curSeq, int &EdgeCoverage, int direction){
	assert(direction == LEFT || direction == RIGHT );
	while( (findPossibleVertex(curKmer, *curSeq, edges, verts) >= 0) ){
		pair <char, EdgePrototype*> otherdir_res = findUniqueWay(edges, curKmer, curSeq, otherDirection(direction));
		pair <char, EdgePrototype*> dir_res = findUniqueWay(edges, curKmer, curSeq, direction);

		if ((otherdir_res.second == NULL) || (dir_res.second == NULL)) {
			break;
		}
		goUniqueWay(edges, curKmer, curSeq, dir_res, EdgeCoverage, direction);
		if (direction == RIGHT) {
			//TODO:: do it, save nucleo/
		}
	}
	return 0;
}
pair<char, EdgePrototype*> findUniqueWay(edgesMap &edges, ll finishKmer, Sequence *finishSeq , int direction){
	assert(direction == LEFT || direction == RIGHT );
	int count = 0;
//	cerr << "findUniqueWay" << endl;
	pair <char, EdgePrototype*> res = make_pair(0, (EdgePrototype *)NULL);

	for (int Nucl = 0; Nucl < 4; Nucl++) {
		ll tmpKmer = pushNucleotide(finishKmer, k - 1, direction, Nucl);
		edgesMap::iterator iter = edges.find(tmpKmer);
		if (iter != edges.end()) {
			for (vector<EdgePrototype *>::iterator it = iter->second.begin(); it != iter->second.end(); ++it) {
				//TODO: minIntersect?
				if (finishSeq->similar(*((*it)->lower), minIntersect, direction)) {
					count++;
					if (count > 1) {
						return make_pair(0, (EdgePrototype *)NULL);
					} else {
						res = make_pair(Nucl, *it);
					}

				}
			}
		}
	}
	return res;
}

//while going left we don't mark anything as used, we just find leftmost possible vert
int goUniqueWay(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq, pair<char, EdgePrototype*> findResult, int &EdgeCoverage, int direction) {
	assert(direction == LEFT || direction == RIGHT );
	ll tmpKmer = pushNucleotide(finishKmer, k, findResult.first, direction);
	finishKmer = popNucleotide(tmpKmer, k + 1, otherDirection(direction));
	bool sameK = false;
	EdgeCoverage += findResult.second->coverage;
	finishSeq = new Sequence(*findResult.second->lower);//PossibleSequence;
	findResult.second->used = 1;
	return 0;
}

int countWays(vector<EdgePrototype *> &v, Sequence *finishSeq, int direction) {
	int count = 0;
//	cerr <<" countWays started"<< endl;
	for (vector<EdgePrototype *>::iterator it = v.begin(); it != v.end(); ++it) {
//TODO: minIntersect?
		if (finishSeq->similar(*((*it)->lower), minIntersect, direction)) {
			count++;
			if (count > 1) {
				return count;
			}
		}
	}
	return count;
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
				protoEdgeType curEdge;
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
				Sequence *finSeq = new Sequence(*startSeq);
				expandDirected(edges, curEdge, graph.verts, finKmer,finSeq, EdgeCoverage, RIGHT);

				int finVertId = findPossibleVertex(startKmer, *startSeq, edges, graph.verts);
				assert(finVertId != -2);
				//TODO: what about loops?
				if (startVertId < 0) {
					//TODO: in  fact, not 0 but ...
					startVertId = storeVertex(graph, startKmer, startSeq, 0);
				}
				if (finVertId < 0) {
					finVertId = storeVertex(graph, finKmer, finSeq, 0);
				}

				Edge* newEdge = new Edge(curEdge, startVertId, finVertId, graph.EdgeId, EdgeCoverage);
				graph.addEdge(newEdge);
				//expandDirected(edges, curEdge, graph.verts, startKmer, startSeq, EdgeCoverage, LEFT);
			}
		}
	}
}



}
