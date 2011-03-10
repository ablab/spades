#ifndef CONSTRUCTHASHTABLE_HPP_
#define CONSTRUCTHASHTABLE_HPP_
#include "common.hpp"
#include "string"
using namespace std;
typedef map<ll, vector<ll> > myMap;

void readsToPairs(string inputFile, string outputFile);
int pairsToSequences(string inputFile, string lmerFile, string outputFile);
string decompress(ll a, int l);
void codeRead(char *read, char *code);
ll extractMer(char *read, int shift, int length);
void addPairToTable(myMap& table, ll upper, ll lower);
int pairsToLmers(string inputFile, string outputFile);

inline bool nextReadPair(char * &read1, char * &read2) {
	return (scanf("%s %s", read1, read2)== 2);
}


#endif /* CONSTRUCTHASHTABLE_HPP_ */
