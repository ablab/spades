/*
 * main.cpp
 *
 *  Created on: Jul 11, 2011
 *      Author: andrey
 */

#include "long_contigs.hpp"
#include "graph.hpp"

namespace long_contigs {

void BuildSeeds(Graph& g, std::vector<BiderectionalPath>& seeds) {
	BuildDeBruijnGraph(g);
	FindSeeds(g, seeds);
}

void PrintSeedEdgeLengthStats(std::vector<BiderectionalPath>& seeds, std::ostream &os) {
	std::map<size_t, size_t> lengthMap;

	for(auto iter = seeds.begin(); iter != seeds.end(); ++iter) {
		++lengthMap[iter->size()];
	}

	for(auto iter = lengthMap.begin(); iter != lengthMap.end(); ++iter) {
		os << iter->first << " : " << iter->second << std::endl;
	}
}

void PrintSeedLengthStats(Graph& g, std::vector<BiderectionalPath>& seeds, std::ostream &os) {
	std::map<size_t, size_t> lengthMap;

	for(auto iter = seeds.begin(); iter != seeds.end(); ++iter) {
		++lengthMap[PathLength(g, *iter)];
	}

	for(auto iter = lengthMap.begin(); iter != lengthMap.end(); ++iter) {
		os << iter->first << " : " << iter->second << std::endl;
	}
}

//void PrintSeedCoverage(Graph& g, std::vector<BiderectionalPath>& seeds, std::ostream &os) {
//	std::multiset<EdgeId> covered;
//
//	for(auto path = seeds.begin(); path != seeds.end(); ++path) {
//		for(BiderectionalPath::iterator iter = path->back(); iter != path->end(); ++iter) {
//			covered.insert(*iter);
//		}
//	}
//
//	std::map<size_t, size_t> coveredTimes;
//	size_t edgeCount;
//	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
//		++coveredTimes[covered.count(*iter)];
//		++edgeCount;
//	}
//
//	for(auto iter = coveredTimes.begin(); iter != coveredTimes.end(); ++iter) {
//		os << iter->first << " : " << iter->second << std::endl;
//	}
//	os << "Total edge coverage: " << (double) coveredTimes[0]/edgeCount << std::endl;
//}

}

int main() {
	using namespace long_contigs;

	Graph g(K);
	std::vector<BiderectionalPath> seeds;

	BuildSeeds(g, seeds);
	PrintSeedLengthStats(g ,seeds, std::cout);

	return 0;
}

