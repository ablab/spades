#pragma once 

#define TR1_UNORDERED 0
#define GOOGLE_DENSE 1
#define MCT_CLOSED_HASH 2

#define PARALLEL_MAP_IN_USE TR1_UNORDERED

#if PARALLEL_MAP_IN_USE == TR1_UNORDERED
    #include <tr1/unordered_set>
#elif PARALLEL_MAP_IN_USE == GOOGLE_DENSE
    #include "google/dense_hash_set"
#elif PARALLEL_MAP_IN_USE == MCT_CLOSED_HASH
    #include "mct/hash-set.hpp"
#endif


template<class T, class Hash, class KeyEqual>
struct parallel_unordered_set
{
    private:
#if PARALLEL_MAP_IN_USE == TR1_UNORDERED
    typedef std::tr1::unordered_set<T, Hash, KeyEqual>               origin_container_t;
#elif PARALLEL_MAP_IN_USE == GOOGLE_DENSE
    typedef google::dense_hash_set<T, Hash, KeyEqual>                origin_container_t;
#elif PARALLEL_MAP_IN_USE == MCT_CLOSED_HASH
    typedef mct::closed_hash_set<T, Hash, KeyEqual>                  origin_container_t;
#endif
        typedef std::vector<origin_container_t>                                 container_arr_t;
        typedef typename origin_container_t::value_type                         value_type;

    public:
        parallel_unordered_set(size_t nthreads, size_t cell_size = 100000)
            : nthreads_     (nthreads)
            , buckets_      (nthreads, origin_container_t(cell_size)) {

#if PARALLEL_MAP_IN_USE == GOOGLE_DENSE
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
        static const size_t LOAD_OVERHEAD = 1000;

        typedef std::vector<T>                                                  origin_container_t;
        typedef std::vector<origin_container_t>                                 container_arr_t;
        typedef typename origin_container_t::value_type                         value_type;

    public:
        parallel_vector(size_t nthreads, size_t cell_size = 100000)
            : nthreads_     (nthreads)
            , cell_size_    (cell_size)
            , buckets_      (nthreads) {

            for (size_t i = 0; i < nthreads_; ++i) {
                buckets_[i].reserve(cell_size + LOAD_OVERHEAD);
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

        origin_container_t & operator[](size_t i)  {
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

        bool is_full() const {
            return buckets_[0].size() >= cell_size_;
        }

        bool is_presisely_full() const {
            for (size_t i = 0; i < nthreads_; ++i) {
                if (buckets_[i].size() >= cell_size_)
                    return true;
            }
            return false;
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
        size_t      cell_size_;
        container_arr_t     buckets_;
};
