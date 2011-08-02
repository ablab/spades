#ifndef HAMMER_KMERSTAT_HPP_
#define HAMMER_KMERSTAT_HPP_

#include <stdint.h>
#include <vector>

#define KMERSTAT_CHANGE	1e12
#define KMERSTAT_GOOD	1e12 + 1
#define KMERSTAT_BAD	1e12 + 2

struct KMerStat {

	KMerStat (size_t cnt, uint64_t cng) : count(cnt), changeto(cng) { }
	size_t count;
	uint64_t changeto;

	bool isGood() const { return changeto == KMERSTAT_GOOD; }
	bool change() const { return changeto < KMERSTAT_CHANGE; }
};

class PositionKMer;

typedef std::map<PositionKMer, KMerStat> KMerStatMap;
typedef std::pair<PositionKMer, KMerStat> KMerCount;


#endif //  HAMMER_KMERSTAT_HPP_

