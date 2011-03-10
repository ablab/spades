#ifndef CONSTRUCTHASHTABLE_HPP_
#define CONSTRUCTHASHTABLE_HPP_
#include "common.hpp"
#include "string"
using namespace std;

void readsToPairs(string inputFile, string outputFile);
int pairsToSequences(string inputFile, string outputFile);
string decompress(ll a, int l);
void codeRead(char *read, char *code);
ll extractMer(char *read, int shift, int length);
<<<<<<< HEAD
=======
void addPairToTable(myMap& table, ll upper, ll lower);
int pairsToLmers(string inputFile, string outputFile);

inline bool nextReadPair(char * &read1, char * &read2) {
	return (scanf("%s %s", read1, read2)== 2);
}


>>>>>>> c9e9a4c50cceaea6b9b70c96835b90cb2ee4055d
#endif /* CONSTRUCTHASHTABLE_HPP_ */
