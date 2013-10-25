#ifndef __SUBCLUSTER_HPP__
#define __SUBCLUSTER_HPP__

#include <vector>

class KMerData;

size_t subcluster(KMerData &kmer_data, std::vector<size_t> &cluster);

#endif // __SUBCLUSTER_HPP__
