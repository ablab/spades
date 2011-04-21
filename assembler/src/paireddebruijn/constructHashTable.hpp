#ifndef CONSTRUCTHASHTABLE_HPP_
#define CONSTRUCTHASHTABLE_HPP_
#include "common.hpp"
#include "string"
using namespace std;
typedef std::tr1::unordered_map<ll, pair<vector<ll>, vector<short> > > myMap;


//reverse- take l-mer from first read and k-mer- from second
pair<int, pair<int,int>> maxCommonSubstring(string &s1,string &s2);
void readsToPairs(string inputFile, string outputFile, bool reverse = false);
int pairsToSequences(string inputFile, string lmerFile, string outputFile);
string decompress(ll a, int l);
//void codeRead(char *read, char *code);
ll extractMer(char *read, int shift, int length);
void addPairToTable(myMap& table, ll upper, ll lower);
int pairsToLmers(string inputFile, string outputFile);
void initGlobal();
void constructReversedReadPairs(string inputFile, string outputFile);
void statsReadPair(map<ll, int>& kmers,map<ll, int>& lmers, char *upperRead, char *lowerRead);
void computeGlobalStats(string kmers_s, string reads, string statsOutputFile);
inline bool ComparePairByFirst(pair<ll,int> i, pair<ll,int> j){
	return i.first < j.first;
}



#endif /* CONSTRUCTHASHTABLE_HPP_ */
