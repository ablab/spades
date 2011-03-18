#include "common.hpp"

#include "constructHashTable.hpp"
#include "graphConstruction.hpp"
#include "sequence.hpp"


int main() {
	initConstants(ini_file);
//	freopen(error_log.c_str(), "w",stderr);
//	assert( 1== 0);
	cerr << l << " " << k;
	freopen(graph_file.c_str(), "w",stdout);
//	LOG_ASSERT(1 == 0, "Something wrong");
	if (needPairs)
		readsToPairs(parsed_reads, parsed_k_l_mers);
	if (needLmers)
		pairsToLmers(parsed_k_l_mers, parsed_l_mers);
	if (needSequences)
		pairsToSequences(parsed_k_l_mers, parsed_l_mers, parsed_k_sequence);
//	map<>sequencesToMap(parsed_k_sequence);
	if (needGraph)
		constructGraph();
//	testSimilar();
//	testFind();
	return 0;
}
