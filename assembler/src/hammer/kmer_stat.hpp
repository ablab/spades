#ifndef HAMMER_KMERSTAT_HPP_
#define HAMMER_KMERSTAT_HPP_

#include <stdint.h>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <bitset>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/bitset.hpp>

const uint32_t K = 55;
typedef uint64_t hint_t;

#define BLOBKMER_UNDEFINED		2e12
#define KMERSTAT_CHANGE			2e12
#define KMERSTAT_BAD			2e12 + 1
#define KMERSTAT_GOOD			2e12 + 2
#define KMERSTAT_GOODITER		2e12 + 3
#define KMERSTAT_GOODITER_BAD		2e12 + 4
#define KMERSTAT_MARKED_FOR_GOODITER	2e12 + 5

#define MAX_SHORT 254

class Read;
class PositionRead;
class PositionKMer;
class KMerStat;

typedef std::map<PositionKMer, KMerStat> KMerStatMap;
typedef std::pair<PositionKMer, KMerStat> KMerCount;
typedef std::pair<std::string, std::pair<uint32_t, double> > StringCount;

struct QualBitSet {
	std::vector<unsigned char> q;
	QualBitSet() : q(K) {
		// q.reset();
	}
	QualBitSet(size_t n) : q(n) {
		// q.reset();
	}
	unsigned short operator[](size_t n) const {
		return (unsigned short)q[n];
	}

	unsigned short at(size_t n) const {
			return (unsigned short)q.at(n);
	}


	void set(size_t n, unsigned short value) {
		q[n] = (unsigned char)value;
	}

	friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & q;
    }
};

struct KMerStat {

	KMerStat (bool first, uint32_t cnt, hint_t cng, double quality) : count(cnt), changeto(cng), totalQual(quality), qual(first ? 0 : K) { }
	// KMerStat (uint32_t cnt, hint_t cng, double quality) : count(cnt), changeto(cng), totalQual(quality), qual() { }
	KMerStat () : count(0), changeto(KMERSTAT_BAD), totalQual(1), qual(0) { }

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

	friend class boost::serialization::access;
	// When the class Archive corresponds to an output archive, the
	// & operator is defined similar to <<.  Likewise, when the class Archive
	// is a type of input archive the & operator is defined similar to >>.
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
	    ar & count;
		ar & changeto;
		ar & totalQual;
		ar & qual;
	}
};

char getQual(const KMerCount & kmc, int i);

#endif //  HAMMER_KMERSTAT_HPP_

