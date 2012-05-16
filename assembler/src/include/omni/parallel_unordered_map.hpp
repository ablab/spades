#pragma once 


//0 -- std::tr1::unordered_map
//1 -- google::dense_hash_map
//2 -- mct::closed_hash_map
#define MAP_IN_USE 2

#if MAP_IN_USE == 0
    #include <tr1/unordered_map>
    #include <tr1/unordered_set>
#elif MAP_IN_USE == 1
    #include "google/dense_hash_map"
    #include "google/dense_hash_set"
#elif MAP_IN_USE == 2
    #include "mct/hash-set.hpp"
    #include "mct/hash-map.hpp"
#endif

template<class T, class Value, class Hash, class KeyEqual>
struct parallel_unordered_map
{
    private:

#if MAP_IN_USE == 0
    typedef std::tr1::unordered_map<T, Value, Hash, KeyEqual>               origin_container_t;
#elif MAP_IN_USE == 1
    typedef google::dense_hash_map<T, Value, Hash, KeyEqual>                origin_container_t;
#elif MAP_IN_USE == 2
    typedef mct::closed_hash_map<T, Value, Hash, KeyEqual>               origin_container_t;
#endif

        typedef	std::vector<origin_container_t>							        container_arr_t;
        typedef typename origin_container_t::value_type                         value_type;

    public:    
        parallel_unordered_map(size_t nthreads, size_t cell_size = 100000)
            : nthreads_		(nthreads)
            , buckets_		(nthreads, origin_container_t(cell_size)) {

#if MAP_IN_USE == 1
            for (size_t i = 0; i < nthreads_; ++i) {
                buckets_[i].set_empty_key(T::GetZero());
            }
#endif
        }

        void insert(const value_type& value, size_t bucket_num)
        {
            buckets_[bucket_num].insert(value);
        }

        const origin_container_t & operator[](size_t i) const
        {
            return buckets_[i];
        }

        size_t get_threads_num() const 
        {
            return nthreads_;
        }

        const container_arr_t & get_buckets() const
        {
            return buckets_;   
        }

        void clear() {
            for (size_t i = 0; i < nthreads_; ++i) {
                buckets_[i].clear();
            }
        }

    private:
        parallel_unordered_map& operator=(const parallel_unordered_map&);

    private:
        size_t 		nthreads_;
        container_arr_t 	buckets_;
};



template<class T, class Hash, class KeyEqual>
struct parallel_unordered_set
{
    private:
#if MAP_IN_USE == 0
    typedef std::tr1::unordered_set<T, Hash, KeyEqual>               origin_container_t;
#elif MAP_IN_USE == 1
    typedef google::dense_hash_set<T, Hash, KeyEqual>                origin_container_t;
#elif MAP_IN_USE == 2
    typedef mct::closed_hash_set<T, Hash, KeyEqual>                  origin_container_t;
#endif
        typedef std::vector<origin_container_t>                                 container_arr_t;
        typedef typename origin_container_t::value_type                         value_type;

    public:
        parallel_unordered_set(size_t nthreads, size_t cell_size = 100000)
            : nthreads_     (nthreads)
            , buckets_      (nthreads, origin_container_t(cell_size)) {

#if MAP_IN_USE == 1
            for (size_t i = 0; i < nthreads_; ++i) {
                buckets_[i].set_empty_key(T::GetZero());
            }
#endif
        }

        void insert(const value_type& value, size_t bucket_num)
        {
            buckets_[bucket_num].insert(value);
        }

        const origin_container_t & operator[](size_t i) const
        {
            return buckets_[i];
        }

        size_t get_threads_num() const
        {
            return nthreads_;
        }

        const container_arr_t & get_buckets() const
        {
            return buckets_;
        }

        void clear() {
            for (size_t i = 0; i < nthreads_; ++i) {
                buckets_[i].clear();
            }
        }

    private:
        parallel_unordered_set& operator=(const parallel_unordered_set&);

    private:
        size_t      nthreads_;
        container_arr_t     buckets_;
};



template<class T>
struct parallel_vector
{
    private:
        typedef std::vector<T>                                                  origin_container_t;
        typedef std::vector<origin_container_t>                                 container_arr_t;
        typedef typename origin_container_t::value_type                         value_type;

    public:
        parallel_vector(size_t nthreads, size_t cell_size = 100000)
            : nthreads_     (nthreads)
            , buckets_      (nthreads) {

            for (size_t i = 0; i < nthreads_; ++i) {
                buckets_[i].reserve(cell_size);
            }
        }

        void insert(const value_type& value, size_t bucket_num)
        {
            buckets_[bucket_num].push_back(value);
        }

        const origin_container_t & operator[](size_t i) const
        {
            return buckets_[i];
        }

        size_t get_threads_num() const
        {
            return nthreads_;
        }

        const container_arr_t & get_buckets() const
        {
            return buckets_;
        }

        void clear() {
            for (size_t i = 0; i < nthreads_; ++i) {
                buckets_[i].clear();
            }
        }

    private:
        parallel_vector& operator=(const parallel_vector&);

    private:
        size_t      nthreads_;
        container_arr_t     buckets_;
};

