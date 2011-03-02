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
#endif /* CONSTRUCTHASHTABLE_HPP_ */
