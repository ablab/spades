#include "constructHashTable.hpp"
#include "graphConstruction.hpp"
#include "../seq.hpp"
#include "common.hpp"



int testFind() {
	Sequence s("ACATACAGACATACA");
	int t = s.find(Sequence("TACAC"));
	cout << t;
	return 0;
}

int testSimilar() {
	Sequence s("ACATACAGACATACA");
	Sequence t("TACAC");
	forn(i, 10) {
		int ii = i+1;
		//cout << ii;
		int k = s.similar(t, ii);
		cout << k;
	}
	return 0;
}

int main() {
	//	freopen("error_log", "w",stderr);
//	readsToPairs(parsed_reads, parsed_k_l_mers);
//	pairsToSequences(parsed_k_l_mers, parsed_k_sequence);
//	map<>sequencesToMap(parsed_k_sequence);
	constructGraph();
//	testSimilar();
	return 0;
}
