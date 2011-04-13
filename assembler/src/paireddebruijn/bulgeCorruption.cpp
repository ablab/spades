#include "bulgeCorruption.hpp"
#include "common.hpp"
#include "map"

#define SIMILARITY_BOUND 2
#define MAX_GLUE_LENGTH 500

using namespace paired_assembler;

void generateEdgeTypes(longEdgesMap &edges,
		map<pair<int, int> , vector<pair<int, Edge*> > > &edgeTypes) {
	for (longEdgesMap::iterator it = edges.begin(); it != edges.end(); it++) {
		Edge* e = it->second;
		int from = e->FromVertex;
		int to = e->ToVertex;
		edgeTypes[make_pair(from, to)].push_back(make_pair(it->first, e));
	}
}

int hamming(Sequence &s1, Sequence &s2) {
	int result = 0;
	for (size_t i = 0; i < s1.size(); i++)
		if (s1[i] != s2[i])
			result++;
	return result;
}

bool similar(Edge e1, Edge e2, int bound, size_t maxLength) {
	return e1.upper->size() == e2.upper->size() && e1.upper->size()
			<= maxLength && hamming(*(e1.upper), *(e2.upper)) + hamming(
			*(e1.lower), *(e2.lower)) <= bound;
}

void removeSimilarEdges(vector<pair<int, Edge *> > edges){
	for (size_t i = 0; i < edges.size(); i++)
		if (edges[i].first == edges[i].second->EdgeId)
			for (size_t j = i + 1; j < edges.size(); j++)
				if (edges[j].first == edges[j].second->EdgeId)
					if (similar(*edges[i].second, *edges[j].second,
							SIMILARITY_BOUND, MAX_GLUE_LENGTH))
						edges[j].second->EdgeId = edges[i].second->EdgeId;
}

void removeBulges(longEdgesMap edges) {
	map<pair<int, int> , vector<pair<int, Edge*> > > edgeTypes;
	generateEdgeTypes(edges, edgeTypes);
	for (map<pair<int, int> , vector<pair<int, Edge*>> >::iterator it =
			edgeTypes.begin(); it != edgeTypes.end(); it++) {
		removeSimilarEdges(it->second);
	}
}
