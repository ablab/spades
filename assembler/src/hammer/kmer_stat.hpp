#ifndef HAMMER_KMERSTAT_HPP_
#define HAMMER_KMERSTAT_HPP_
#include <vector>
struct KMerStat {
	size_t count;
	bool change;
	bool good;
	uint64_t changeto;
	vector< pair<int64_t, uint32_t> > pos;  // positions in reads
};

#endif //  HAMMER_KMERSTAT_HPP_
