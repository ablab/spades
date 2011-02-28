#include "constructHashTable.hpp"
#include "graphConstruction.hpp"
#include "../seq.hpp"
#include "common.h"
string parsed_reads = "data/reads_var_d.txt";
string parsed_k_l_mers = "data/klmers_var_d.txt";
string parsed_k_sequence = "data/vertices_var_d.txt";
string error_log = "data/error.log";

int main() {
	//	freopen("error_log", "w",stderr);

	readsToPairs(parsed_reads, parsed_k_l_mers);
	pairsToSequences(parsed_k_l_mers, parsed_k_sequence);


}
