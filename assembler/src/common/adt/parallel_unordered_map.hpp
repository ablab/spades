//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once 

#include <unordered_set>

template<class T, class Hash, class KeyEqual>
struct parallel_unordered_set
{
private:

    typedef std::unordered_set<T, Hash, KeyEqual>                      origin_container_t;

    typedef std::vector<origin_container_t>                                 container_arr_t;

    typedef typename origin_container_t::value_type                         value_type;

    public:
        parallel_unordered_set(size_t nthreads, size_t cell_size = 100000)
            : nthreads_     (nthreads)
            , buckets_      (nthreads, origin_container_t(cell_size)) {

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
