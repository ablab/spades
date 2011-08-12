#ifndef HAMMER_KMERSTAT_HPP_
#define HAMMER_KMERSTAT_HPP_

#include <stdint.h>
#include <vector>

#define BLOBKMER_UNDEFINED	2e12
#define KMERSTAT_CHANGE		2e12
#define KMERSTAT_GOOD		2e12 + 1
#define KMERSTAT_BAD		2e12 + 2

typedef uint64_t hint_t;

struct KMerStat {

	KMerStat (uint32_t cnt, hint_t cng, double qual) : count(cnt), changeto(cng), totalQual(qual) { }
	uint32_t count;
	hint_t changeto;
	double totalQual;

	bool isGood() const { return changeto == KMERSTAT_GOOD; }
	bool change() const { return changeto < KMERSTAT_CHANGE; }
};

class PositionKMer;

typedef std::map<PositionKMer, KMerStat> KMerStatMap;
typedef std::pair<PositionKMer, KMerStat> KMerCount;


#endif //  HAMMER_KMERSTAT_HPP_

