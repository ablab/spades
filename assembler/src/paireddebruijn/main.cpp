#include "constructHashTable.hpp"
#include "graphConstruction.hpp"
#include "../seq.hpp"
#include "common.h"
string parsed_reads = "data/reads_var_d.txt";
string parsed_k_l_mers = "data/klmers_var_d.txt";
string parsed_k_sequence = "data/vertices_var_d.txt";

int main() {
	readsToPairs(parsed_reads, parsed_k_l_mers);
//	pairsToSequences();
	//forn(i, lsize)
	//	cerr << lmers[i] << ":" << decompress(lmers[i]) << " ";
	//cerr << endl << endl;

//	Sequence* tst = new Sequence("01002010301020");
//	string s = tst->str();
//	cout<<s.c_str()<<endl;
}
