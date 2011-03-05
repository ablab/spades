#include "constructHashTable.hpp"
#include "graphConstruction.hpp"
#include "sequence.hpp"
#include "common.hpp"



int main() {
	//freopen(error_log.c_str(), "w",stderr);
	freopen("data/graph.dot", "w",stdout);
//	LOG_ASSERT(1 == 0, "Something wrong");
//	readsToPairs(parsed_reads, parsed_k_l_mers);
//	pairsToLmers(parsed_k_l_mers, parsed_l_mers);
//	pairsToSequences(parsed_k_l_mers, parsed_l_mers, parsed_k_sequence);
//	map<>sequencesToMap(parsed_k_sequence);
	constructGraph();
//	testSimilar();
//	testFind();
	return 0;
}
