#ifndef HAMMER_KMERSTAT_HPP_
#define HAMMER_KMERSTAT_HPP_

#include <stdint.h>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <bitset>

const uint32_t K = 55;
typedef uint64_t hint_t;

#define BLOBKMER_UNDEFINED		2e12
#define KMERSTAT_CHANGE			2e12
#define KMERSTAT_BAD			2e12 + 1
#define KMERSTAT_GOOD			2e12 + 2
#define KMERSTAT_GOODITER		2e12 + 3
#define KMERSTAT_GOODITER_BAD		2e12 + 4
#define KMERSTAT_MARKED_FOR_GOODITER	2e12 + 5

#define MAX_SHORT 1000

class Read;
class PositionRead;
class PositionKMer;
class KMerStat;

typedef std::map<PositionKMer, KMerStat> KMerStatMap;
typedef std::pair<PositionKMer, KMerStat> KMerCount;
typedef std::pair<std::string, std::pair<uint32_t, double> > StringCount;

struct QualBitSet {
	std::bitset<K*10> q;
	QualBitSet() { q.reset(); }
	unsigned short operator[](size_t n) const {
		unsigned short x = 0;
		for (size_t i=0; i<10; ++i) {
			x = x << 1;
			if (q[n*10+9-i]) ++x;
		}
		return x;
	}
	void set(size_t n, unsigned short value) {
		unsigned short x = value;
		for (size_t j = n*10; j < n*10+10; ++j) {
			q[j] = x % 2;
			x = x >> 1;
		}
	}
};

struct KMerStat {

	KMerStat (uint32_t cnt, hint_t cng, double quality) : count(cnt), changeto(cng), totalQual(quality), qual() { }
	KMerStat () : count(0), changeto(KMERSTAT_BAD), totalQual(1), qual() { }

	uint32_t count;
	hint_t changeto;
	double totalQual;
	QualBitSet qual;

	bool isGood() const { return changeto >= KMERSTAT_GOOD; }
	bool isGoodForIterative() const { return (changeto == KMERSTAT_GOODITER); }
	void makeGoodForIterative() { changeto = KMERSTAT_GOODITER; }
	void markGoodForIterative() { changeto = KMERSTAT_MARKED_FOR_GOODITER; }
	bool isMarkedGoodForIterative() { return (changeto == KMERSTAT_MARKED_FOR_GOODITER); }
	bool change() const { return changeto < KMERSTAT_CHANGE; }
};

#endif //  HAMMER_KMERSTAT_HPP_

