#include "constructHashTable.hpp"
#include "common.hpp"

void constructGraph() {
	readsToPairs(parsed_reads, parsed_k_l_mers);
	pairsToSequences(parsed_k_l_mers, parsed_k_sequence);
//	map<ll, vector<Sequence *> > sequencesToMap(parsed_k_sequence);
//	Graph * g = new GraphBuilder();
}

