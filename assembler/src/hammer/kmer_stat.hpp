#ifndef HAMMER_KMERSTAT_HPP_
#define HAMMER_KMERSTAT_HPP_

#include <stdint.h>
#include <vector>

#define KMERSTAT_CHANGE	18446744073709551610
#define KMERSTAT_GOOD	18446744073709551611
#define KMERSTAT_BAD	18446744073709551612

struct KMerStat {

	KMerStat (size_t cnt, uint64_t cng) : count(cnt), changeto(cng) { }
	size_t count;
	uint64_t changeto;
	// std::vector< std::pair<int64_t, uint32_t> > pos;  // positions in reads

	bool isGood() const { return changeto == KMERSTAT_GOOD; }
	bool change() const { return changeto < KMERSTAT_CHANGE; }
};

class PositionKMer;

typedef std::map<PositionKMer, KMerStat> KMerStatMap;
typedef std::pair<PositionKMer, KMerStat> KMerCount;


#endif //  HAMMER_KMERSTAT_HPP_

