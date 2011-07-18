#ifndef HAMMER_KMERSTAT_HPP_
#define HAMMER_KMERSTAT_HPP_
#include <vector>
struct KMerStat {
	size_t count;
	float freq;
	bool change;
	uint64_t changeto;
	vector< pair<uint64_t, uint32_t> > pos;  // positions in reads
};
#endif //  HAMMER_KMERSTAT_HPP_
