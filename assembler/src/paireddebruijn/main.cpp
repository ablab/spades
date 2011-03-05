#include "constructHashTable.hpp"
#include "graphConstruction.hpp"
#include "../sequence.hpp"
#include "common.hpp"



int testFind() {
	Sequence s("ACATACAGACATACA");
	cerr<<s.str()<<endl;
	Sequence ss = s.Subseq(5,10);
	cerr<<s.str()<<" "<<ss.str()<<endl;
	int t = s.find(ss);
	cout << t;
	return 0;
}

int testSimilar() {
	Sequence s("AAACTCGAAA");
	Sequence t("TTAAACTCGA");
	forn(i, 10) {
		int ii = i+1;
		//cout << ii;
		int k = s.similar(t, ii, -1);
		cerr << "k: "<< ii << " "<<k<<endl;
	}
	return 0;
}

int main() {
	freopen(error_log.c_str(), "w",stderr);
	freopen("data/graph.dot", "w",stdout);

	//readsToPairs(parsed_reads, parsed_k_l_mers);
//	pairsToSequences(parsed_k_l_mers, parsed_k_sequence);
//	map<>sequencesToMap(parsed_k_sequence);
	constructGraph();
//	testSimilar();
//	testFind();
	return 0;
}
