#pragma once
#include "omni/parallel_unordered_map.hpp"
#include <omp.h>

#define DESTINATION_MAP 2

//0 -- std::tr1::unordered_map
//2 -- mct::closed_hash_map

#if DESTINATION_MAP == 0
    #include <tr1/unordered_map>
    #include <tr1/unordered_set>
#elif DESTINATION_MAP == 2
    #include "mct/hash-set.hpp"
    #include "mct/hash-map.hpp"
#endif



using std::make_pair;
using std::pair;

template <size_t size_, typename Value> 
class ParallelSeqMap {

public:
    typedef Seq<size_> Kmer;

    typedef parallel_unordered_map<Kmer, pair<Value, size_t>, typename Kmer::hash, typename Kmer::equal_to> par_container_t;

#if DESTINATION_MAP == 0
    typedef std::tr1::unordered_map<Kmer, pair<Value, size_t>, typename Kmer::hash, typename Kmer::equal_to> destination_container_t; // size_t is offset
#elif DESTINATION_MAP == 2
    typedef mct::closed_hash_map<Kmer, pair<Value, size_t>, typename Kmer::hash, typename Kmer::equal_to> destination_container_t;
#endif

private:
	
    size_t nthreads_;
	std::vector<par_container_t> nodes_;

public:

    ParallelSeqMap(size_t nthreads, size_t cell_size = 100000) :
        nthreads_(nthreads),
        nodes_(nthreads, par_container_t(nthreads, cell_size))
    {
    }
	
    void AddEdge(const Kmer &k, size_t i) {
		nodes_[i].insert(make_pair(k, make_pair(Value(), -1)), k.GetHash() % nthreads_);
    }

	void CountSequence(const Sequence& s, size_t thread_number) {
        if (s.size() < size_)
			return;

		Kmer kmer = s.start<size_>();
		AddEdge(kmer, thread_number);
		for (size_t j = size_; j < s.size(); ++j) {
			kmer = kmer << s[j];
			AddEdge(kmer, thread_number);
		}
	}

    void MergeMaps(destination_container_t & temp_map, size_t i) {
        if (temp_map.bucket_count() < nodes_[0][i].bucket_count()) {
            temp_map.rehash(nodes_[0][i].bucket_count() * nthreads_);
        }

        for (size_t j = 0; j < nthreads_; ++j) {
            temp_map.insert(nodes_[j][i].begin(), nodes_[j][i].end());
        }
    }

    bool Contains(const Sequence& s) {
        for (size_t i = 0; i < nthreads_; ++i) 
            if (nodes_[i].find(s) != nodes_[i].end()) return true;
        return false;
    }


    size_t SingleBucketCount() const {
        return nodes_[0][0].bucket_count();
    }

    size_t GetSize(size_t i, size_t j) const {
        return nodes_[i][j].size();
    }

    void clear() {
        for (size_t i = 0; i < nthreads_; ++i) {
            nodes_[i].clear();
        }
    }
};


template <size_t size_>
class ParallelSeqSet {

public:
    typedef Seq<size_> Kmer;

    typedef parallel_unordered_set<Kmer, typename Kmer::hash, typename Kmer::equal_to> par_container_t;

#if DESTINATION_MAP == 0
    typedef std::tr1::unordered_set<Kmer, typename Kmer::hash, typename Kmer::equal_to> destination_container_t;
#elif DESTINATION_MAP == 2
    typedef mct::closed_hash_set<Kmer, typename Kmer::hash, typename Kmer::equal_to> destination_container_t;
#endif


private:

    size_t nthreads_;
    std::vector<par_container_t> nodes_;

public:

    ParallelSeqSet(size_t nthreads, size_t cell_size = 100000) :
        nthreads_(nthreads),
        nodes_(nthreads, par_container_t(nthreads, cell_size))
    {
    }

    void AddEdge(const Kmer &k, size_t i) {
        nodes_[i].insert(k, k.GetHash() % nthreads_);
    }

    void CountSequence(const Sequence& s, size_t thread_number) {
        if (s.size() < size_)
            return;

        Kmer kmer = s.start<size_>();
        AddEdge(kmer, thread_number);
        for (size_t j = size_; j < s.size(); ++j) {
            kmer = kmer << s[j];
            AddEdge(kmer, thread_number);
        }
    }

    void MergeMaps(destination_container_t & temp_map, size_t i) {
        if (temp_map.bucket_count() < nodes_[0][i].bucket_count()) {
            temp_map.rehash(nodes_[0][i].bucket_count() * nthreads_);
        }

        for (size_t j = 0; j < nthreads_; ++j) {
            temp_map.insert(nodes_[j][i].begin(), nodes_[j][i].end());
        }
    }

    bool Contains(const Sequence& s) {
        for (size_t i = 0; i < nthreads_; ++i)
            if (nodes_[i].find(s) != nodes_[i].end()) return true;
        return false;
    }


    size_t SingleBucketCount() const {
        return nodes_[0][0].bucket_count();
    }

    void clear() {
        for (size_t i = 0; i < nthreads_; ++i) {
            nodes_[i].clear();
        }
    }
};


template <size_t size_>
class ParallelSeqVector {

public:
    typedef Seq<size_> Kmer;

    typedef parallel_vector<Kmer> par_container_t;

#if DESTINATION_MAP == 0
    typedef std::tr1::unordered_set<Kmer, typename Kmer::hash, typename Kmer::equal_to> destination_container_t;
#elif DESTINATION_MAP == 2
    typedef mct::closed_hash_set<Kmer, typename Kmer::hash, typename Kmer::equal_to> destination_container_t;
#endif


private:

    size_t nthreads_;
    std::vector<par_container_t> nodes_;

public:

    ParallelSeqVector(size_t nthreads, size_t cell_size = 100000) :
        nthreads_(nthreads),
        nodes_(nthreads, par_container_t(nthreads, cell_size))
    {
    }

    void AddEdge(const Kmer &k, size_t i) {
        nodes_[i].insert(k, k.GetHash() % nthreads_);
    }

    void CountSequence(const Sequence& s, size_t thread_number) {
        if (s.size() < size_)
            return;

        Kmer kmer = s.start<size_>();
        AddEdge(kmer, thread_number);
        for (size_t j = size_; j < s.size(); ++j) {
            kmer = kmer << s[j];
            AddEdge(kmer, thread_number);
        }
    }

    void MergeMaps(destination_container_t & temp_map, size_t i) {
        for (size_t j = 0; j < nthreads_; ++j) {
            temp_map.insert(nodes_[j][i].begin(), nodes_[j][i].end());
        }
    }

    void Dump(destination_container_t & temp_map, size_t i) {
        for (size_t j = 0; j < nthreads_; ++j) {
            temp_map.insert(nodes_[j][i].begin(), nodes_[j][i].end());
            nodes_[j][i].clear();
        }
    }

    bool Contains(const Sequence& s) {
        for (size_t i = 0; i < nthreads_; ++i)
            if (nodes_[i].find(s) != nodes_[i].end()) return true;
        return false;
    }


    size_t SingleBucketCount() const {
        return nodes_[0][0].capacity();
    }

    bool IsFull(size_t i) const {
        return nodes_[i].is_full();
    }

    void Clear(size_t i) {
        nodes_[i].clear();
    }

    void Clear() {
        for (size_t i = 0; i < nthreads_; ++i) {
            nodes_[i].clear();
        }
    }

};
