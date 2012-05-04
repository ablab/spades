#pragma once
#include "omni/parallel_unordered_map.hpp"
#include <omp.h>

template <size_t size_, typename Value> 
class ParallelSeqMap {
public:
    typedef Seq<size_> Kmer;
    typedef parallel_unordered_map<Kmer, pair<Value, size_t>, typename Kmer::hash, typename Kmer::equal_to> par_map_t;
    typedef std::tr1::unordered_map<Kmer, pair<Value, size_t>, typename Kmer::hash, typename Kmer::equal_to> map_t; // size_t is offset

private:
	
    size_t nthreads_;
	std::vector<par_map_t> nodes;

public:

    ParallelSeqMap(size_t nthreads) : nthreads_(nthreads), nodes(nthreads, par_map_t(nthreads)) 
    {
    }
	
    void AddEdge(const Kmer &k, size_t i) {
        //TRACE("Adding ''edge'' to map");
		nodes[i].insert(make_pair(k, make_pair(Value(), -1)), k.GetHash() & (nthreads_ - 1));
    }

	void CountSequence(const Sequence& s, size_t thread_number) {
        if (s.size() < size_)
			return;
        //TRACE("Threading sequence " << s.str());
        //TRACE("Thread number " << thread_number << " or " << omp_get_thread_num()); 
		Kmer kmer = s.start<size_>();
		AddEdge(kmer, thread_number);
        //TRACE("start kmer " << kmer.str());
		for (size_t j = size_; j < s.size(); ++j) {
			kmer = kmer << s[j];
            //TRACE("shifted left " << kmer.str() << " thread_number " << thread_number << " " << char('a'+s[j]));
			AddEdge(kmer, thread_number);
            //TRACE("edge has been added");
		}
        //TRACE("Threading finished");
	}

    void MergeMaps(map_t & temp_map, size_t i) {
        if (temp_map.bucket_count() < nodes[0][i].bucket_count()) {
            temp_map.rehash(nodes[0][i].bucket_count() * nthreads_);
        }

        for (size_t j = 0; j < nthreads_; ++j) {
            temp_map.insert(nodes[j][i].begin(), nodes[j][i].end());
        }
    }

    bool Contains(const Sequence& s) {
        for (size_t i = 0; i < nthreads_; ++i) 
            if (nodes[i].find(s) != nodes[i].end()) return true;
        return false;
    }


    size_t SingleBucketCount() const {
        return nodes[0][0].bucket_count();
    }

    void clear() {
        for (size_t i = 0; i < nthreads_; ++i) {
            nodes[i].clear();
        }
    }
};
