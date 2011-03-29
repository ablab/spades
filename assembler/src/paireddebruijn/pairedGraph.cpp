#include "pairedGraph.hpp"
#include "common.hpp"
#include "graphio.hpp"
//
//int isDeleted(Node* A) {
//	return A->upperSize == 0;
//}
//int deleteNode(Node* A) {
//	return A->upperSize = 0;
//}
LOGGER("p.pairedGraph");

namespace paired_assembler {

VertexPrototype::VertexPrototype(Sequence *lower_, int id, int coverage_) {
	upper = NULL;
	TRACE("Upper sequence is not defined");
	lower = lower_;
	VertexId = id;
	used = false;
	coverage = coverage_;
}

VertexPrototype::VertexPrototype(ll upper_, Sequence *lower_, int id,
		int coverage_) {
	upper = upper_;
	lower = lower_;
	VertexId = id;
	used = false;
	coverage = coverage_;
}

void PairedGraph::recreateVerticesInfo(int vertCount, longEdgesMap &longEdges) {
	forn(i, vertCount) {
		forn(j, 2)
			degrees[i][j] = 0;
	}
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
			edgeIds[it->second->FromVertex][degrees[it->second->FromVertex][1]++][OUT_EDGE]
					= it->first;
			edgeIds[it->second->ToVertex][degrees[it->second->ToVertex][0]++][IN_EDGE]
					= it->first;
		}
	}
}

//todo: Complete this
void PairedGraph::RebuildVertexMap(void) {
	verts.clear();
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
			ll kmer = 0;
			forn(i,k-1)
				kmer = (kmer << 2) | (it->second->upper->operator[](i));
			Sequence * seq = new Sequence(it->second->lower->Subseq(0, l - 1));
			storeVertex(*this, kmer, seq, it->second->FromVertex);
			kmer = 0;
			forn(i,k-1)
				kmer = (kmer << 2) | (it->second->upper->operator[](
						i + it->second->length));
			seq = new Sequence(
					it->second->lower->Subseq(it->second->length,
							it->second->length + l - 1));
			storeVertex(*this, kmer, seq, it->second->ToVertex);
		}
	}
}

verticesMap::iterator addKmerToMap(verticesMap &verts, ll kmer) {
	verticesMap::iterator position = verts.find(kmer);
	if (position == verts.end()) {
		vector<VertexPrototype *> prototypes;
		return verts.insert(make_pair(kmer, prototypes)).first;
	} else {
		return position;
	}
}

/*First argument of result is id of the vertex. Second argument is true if new entry
 * was created and false otherwise
 *
 */
pair<int, bool> addVertexToMap(PairedGraph &graph, ll newKmer, Sequence* newSeq) {
	verticesMap::iterator position = addKmerToMap(graph.verts, newKmer);
	vector<VertexPrototype *> *sequences = &position->second;
	for (vector<VertexPrototype *>::iterator it = sequences->begin(); it
			!= sequences->end(); ++it) {
		Sequence *otherSequence = (*it)->lower;
		if (newSeq->similar(*otherSequence, minIntersect, 0)) {
			return make_pair((*it)->VertexId, false);
		}
	}
	sequences->push_back(new VertexPrototype(newSeq, graph.VertexCount));
	graph.VertexCount++;
	return make_pair(graph.VertexCount - 1, true);
}
pair<int, bool> addVertexToMap(PairedGraph &graph, ll newKmer,
		Sequence* newSeq, int VertNum) {
	verticesMap::iterator position = addKmerToMap(graph.verts, newKmer);
	vector<VertexPrototype *> *sequences = &position->second;
	for (vector<VertexPrototype *>::iterator it = sequences->begin(); it
			!= sequences->end(); ++it) {
		Sequence *otherSequence = (*it)->lower;
		if (newSeq->similar(*otherSequence, minIntersect, 0)) {
			return make_pair((*it)->VertexId, false);
		}
	}
	sequences->push_back(new VertexPrototype(newSeq, VertNum));
	return make_pair(VertNum, true);
}

int storeVertex(gvis::GraphPrinter<int> &g, PairedGraph &graph, ll newKmer,
		Sequence* newSeq) {
	pair<int, bool> addResult = addVertexToMap(graph, newKmer, newSeq);
	if (addResult.second) {
		g.addVertex(addResult.first, decompress(newKmer, k - 1));
	}
	return addResult.first;
}
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq) {
	pair<int, bool> addResult = addVertexToMap(graph, newKmer, newSeq);
	return addResult.first;
}
int storeVertex(PairedGraph &graph, ll newKmer, Sequence* newSeq, int VertNum) {
	pair<int, bool> addResult = addVertexToMap(graph, newKmer, newSeq, VertNum);
	return addResult.first;
}

void resetVertexCount(PairedGraph &graph) {
	graph.VertexCount = 0;
}

}
