/*
    pdqsort.h - Pattern-defeating quicksort.

    Copyright (c) 2015 Orson Peters

    This software is provided 'as-is', without any express or implied warranty. In no event will the
    authors be held liable for any damages arising from the use of this software.

    Permission is granted to anyone to use this software for any purpose, including commercial
    applications, and to alter it and redistribute it freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not claim that you wrote the
       original software. If you use this software in a product, an acknowledgment in the product
       documentation would be appreciated but is not required.

    2. Altered source versions must be plainly marked as such, and must not be misrepresented as
       being the original software.

    3. This notice may not be removed or altered from any source distribution.
*/


#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <utility>
#include <iterator>

namespace pdqsort_pod_detail {
    enum {
        // Partitions below this size are sorted using insertion sort.
        insertion_sort_threshold = 24,

        // Partitions above this size use Tukey's ninther to select the pivot.
        ninther_threshold = 128,

        // When we detect an already sorted partition, attempt an insertion sort that allows this
        // amount of element moves before giving up.
        partial_insertion_sort_limit = 8,

        // Must be multiple of 8 due to loop unrolling, and < 256 to fit in unsigned char.
        block_size = 64,

        // Cacheline size, assumes power of two.
        cacheline_size = 64

    };

    // Returns floor(log2(n)), assumes n > 0.
    template<class T>
    inline int log2(T n) {
        int log = 0;
        while (n >>= 1) ++log;
        return log;
    }

    template<class WrappedIter>
    void array_swap(WrappedIter lhs, WrappedIter rhs,
                    size_t array_n) {
        // Note that iterator here are wrapped, so we need to "unwrap" them
        typename std::iterator_traits<WrappedIter>::pointer lhs_ptr = &*lhs, rhs_ptr = &*rhs;
        for (size_t i = 0; i < array_n; ++i) {
            using namespace std;
            swap(lhs_ptr[i], rhs_ptr[i]);
        }
    }

    template<class Ptr>
    void array_move(Ptr lhs, Ptr rhs,
                    size_t array_n) {
        for (size_t i = 0; i < array_n; ++i)
            lhs[i] = rhs[i];
    }

    // Sorts [begin, end) using insertion sort with the given comparison function.
    template<class Iter, class Compare>
    inline void insertion_sort(Iter begin, Iter end, size_t array_n, Compare comp) {
        typedef typename std::iterator_traits<Iter>::value_type T;
        if (begin == end) return;

        T tmp[array_n];
        for (Iter cur = begin + 1; cur != end; ++cur) {
            Iter sift = cur;
            Iter sift_1 = cur - 1;

            // Compare first so we can avoid 2 moves for an element already positioned correctly.
            if (comp(&*sift, &*sift_1, array_n)) {
                array_move(tmp, &*sift, array_n);

                do { array_move(&*sift--, &*sift_1, array_n); }
                while (sift != begin && comp(tmp, &*--sift_1, array_n));

                array_move(&*sift, tmp, array_n);
            }
        }
    }

    // Sorts [begin, end) using insertion sort with the given comparison function. Assumes
    // *(begin - 1) is an element smaller than or equal to any element in [begin, end).
    template<class Iter, class Compare>
    inline void unguarded_insertion_sort(Iter begin, Iter end, size_t array_n, Compare comp) {
        typedef typename std::iterator_traits<Iter>::value_type T;
        if (begin == end) return;

        T tmp[array_n];
        for (Iter cur = begin + 1; cur != end; ++cur) {
            Iter sift = cur;
            Iter sift_1 = cur - 1;

            // Compare first so we can avoid 2 moves for an element already positioned correctly.
            if (comp(&*sift, &*sift_1, array_n)) {
                array_move(tmp, &*sift, array_n);

                do { array_move(&*sift--, &*sift_1, array_n); }
                while (comp(tmp, &*--sift_1, array_n));

                array_move(&*sift, tmp, array_n);
            }
        }
    }

    // Attempts to use insertion sort on [begin, end). Will return false if more than
    // partial_insertion_sort_limit elements were moved, and abort sorting. Otherwise it will
    // successfully sort and return true.
    template<class Iter, class Compare>
    inline bool partial_insertion_sort(Iter begin, Iter end, size_t array_n, Compare comp) {
        typedef typename std::iterator_traits<Iter>::value_type T;
        if (begin == end) return true;

        std::size_t limit = 0;
        T tmp[array_n];
        for (Iter cur = begin + 1; cur != end; ++cur) {
            Iter sift = cur;
            Iter sift_1 = cur - 1;

            // Compare first so we can avoid 2 moves for an element already positioned correctly.
            if (comp(&*sift, &*sift_1, array_n)) {
                array_move(tmp, &*sift, array_n);

                do { array_move(&*sift--, &*sift_1, array_n); }
                while (sift != begin && comp(tmp, &*--sift_1, array_n));

                array_move(&*sift, tmp, array_n);
                limit += cur - sift;
            }

            if (limit > partial_insertion_sort_limit) return false;
        }

        return true;
    }

    template<class Iter, class Compare>
    inline void sort2(Iter a, Iter b, size_t array_n, Compare comp) {
        if (comp(&*b, &*a, array_n)) array_swap(a, b, array_n);
    }

    // Sorts the elements *a, *b and *c using comparison function comp.
    template<class Iter, class Compare>
    inline void sort3(Iter a, Iter b, Iter c, size_t array_n, Compare comp) {
        sort2(a, b, array_n, comp);
        sort2(b, c, array_n, comp);
        sort2(a, b, array_n, comp);
    }

    template<class T>
    inline T* align_cacheline(T* p) {
#if defined(UINTPTR_MAX) && __cplusplus >= 201103L
        std::uintptr_t ip = reinterpret_cast<std::uintptr_t>(p);
#else
        std::size_t ip = reinterpret_cast<std::size_t>(p);
#endif
        ip = (ip + cacheline_size - 1) & -cacheline_size;
        return reinterpret_cast<T*>(ip);
    }

    template<class Iter>
    inline void swap_offsets(Iter first, Iter last, size_t array_n,
                             unsigned char* offsets_l, unsigned char* offsets_r,
                             int num, bool use_swaps) {
        typedef typename std::iterator_traits<Iter>::value_type T;
        if (use_swaps) {
            // This case is needed for the descending distribution, where we need
            // to have proper swapping for pdqsort to remain O(n).
            for (size_t i = 0; i < num; ++i) {
                array_swap(first + offsets_l[i], last - offsets_r[i], array_n);
            }
        } else if (num > 0) {
            Iter l = first + offsets_l[0]; Iter r = last - offsets_r[0];
            T tmp[array_n];
            array_move(tmp, &*l, array_n);
            array_move(&*l, &*r, array_n);
            for (size_t i = 1; i < num; ++i) {
                l = first + offsets_l[i]; array_move(&*r, &*l, array_n);
                r = last - offsets_r[i]; array_move(&*l, &*r, array_n);
            }
            array_move(&*r, tmp, array_n);
        }
    }

    // Partitions [begin, end) around pivot *begin using comparison function comp. Elements equal
    // to the pivot are put in the right-hand partition. Returns the position of the pivot after
    // partitioning and whether the passed sequence already was correctly partitioned. Assumes the
    // pivot is a median of at least 3 elements and that [begin, end) is at least
    // insertion_sort_threshold long. Uses branchless partitioning.
    template<class Iter, class Compare>
    inline std::pair<Iter, bool> partition_right_branchless(Iter begin, Iter end, size_t array_n, Compare comp) {
        typedef typename std::iterator_traits<Iter>::value_type T;

        // Move pivot into local for speed.
        T pivot[array_n];
        array_move(pivot, &*begin, array_n);
        Iter first = begin;
        Iter last = end;

        // Find the first element greater than or equal than the pivot (the median of 3 guarantees
        // this exists).
        while (comp(&*++first, pivot, array_n));

        // Find the first element strictly smaller than the pivot. We have to guard this search if
        // there was no element before *first.
        if (first - 1 == begin) while (first < last && !comp(&*--last, pivot, array_n));
        else                    while (                !comp(&*--last, pivot, array_n));

        // If the first pair of elements that should be swapped to partition are the same element,
        // the passed in sequence already was correctly partitioned.
        bool already_partitioned = first >= last;
        if (!already_partitioned) {
            array_swap(first, last, array_n);
            ++first;

            // The following branchless partitioning is derived from "BlockQuicksort: How Branch
            // Mispredictions don’t affect Quicksort" by Stefan Edelkamp and Armin Weiss, but
            // heavily micro-optimized.
            unsigned char offsets_l_storage[block_size + cacheline_size];
            unsigned char offsets_r_storage[block_size + cacheline_size];
            unsigned char* offsets_l = align_cacheline(offsets_l_storage);
            unsigned char* offsets_r = align_cacheline(offsets_r_storage);

            Iter offsets_l_base = first;
            Iter offsets_r_base = last;
            size_t num_l, num_r, start_l, start_r;
            num_l = num_r = start_l = start_r = 0;
            
            while (first < last) {
                // Fill up offset blocks with elements that are on the wrong side.
                // First we determine how much elements are considered for each offset block.
                size_t num_unknown = last - first;
                size_t left_split = num_l == 0 ? (num_r == 0 ? num_unknown / 2 : num_unknown) : 0;
                size_t right_split = num_r == 0 ? (num_unknown - left_split) : 0;

                // Fill the offset blocks.
                if (left_split >= block_size) {
                    for (size_t i = 0; i < block_size;) {
                        offsets_l[num_l] = i++; num_l += !comp(&*first, pivot, array_n); ++first;
                        offsets_l[num_l] = i++; num_l += !comp(&*first, pivot, array_n); ++first;
                        offsets_l[num_l] = i++; num_l += !comp(&*first, pivot, array_n); ++first;
                        offsets_l[num_l] = i++; num_l += !comp(&*first, pivot, array_n); ++first;
                        offsets_l[num_l] = i++; num_l += !comp(&*first, pivot, array_n); ++first;
                        offsets_l[num_l] = i++; num_l += !comp(&*first, pivot, array_n); ++first;
                        offsets_l[num_l] = i++; num_l += !comp(&*first, pivot, array_n); ++first;
                        offsets_l[num_l] = i++; num_l += !comp(&*first, pivot, array_n); ++first;
                    }
                } else {
                    for (size_t i = 0; i < left_split;) {
                        offsets_l[num_l] = i++; num_l += !comp(&*first, pivot, array_n); ++first;
                    }
                }

                if (right_split >= block_size) {
                    for (size_t i = 0; i < block_size;) {
                        offsets_r[num_r] = ++i; num_r += comp(&*--last, pivot, array_n);
                        offsets_r[num_r] = ++i; num_r += comp(&*--last, pivot, array_n);
                        offsets_r[num_r] = ++i; num_r += comp(&*--last, pivot, array_n);
                        offsets_r[num_r] = ++i; num_r += comp(&*--last, pivot, array_n);
                        offsets_r[num_r] = ++i; num_r += comp(&*--last, pivot, array_n);
                        offsets_r[num_r] = ++i; num_r += comp(&*--last, pivot, array_n);
                        offsets_r[num_r] = ++i; num_r += comp(&*--last, pivot, array_n);
                        offsets_r[num_r] = ++i; num_r += comp(&*--last, pivot, array_n);
                    }
                } else {
                    for (size_t i = 0; i < right_split;) {
                        offsets_r[num_r] = ++i; num_r += comp(&*--last, pivot, array_n);
                    }
                }

                // Swap elements and update block sizes and first/last boundaries.
                size_t num = std::min(num_l, num_r);
                swap_offsets(offsets_l_base, offsets_r_base, array_n,
                             offsets_l + start_l, offsets_r + start_r,
                             num, num_l == num_r);
                num_l -= num; num_r -= num;
                start_l += num; start_r += num;

                if (num_l == 0) {
                    start_l = 0;
                    offsets_l_base = first;
                }
                
                if (num_r == 0) {
                    start_r = 0;
                    offsets_r_base = last;
                }
            }

            // We have now fully identified [first, last)'s proper position. Swap the last elements.
            if (num_l) {
                offsets_l += start_l;
                while (num_l--) array_swap(offsets_l_base + offsets_l[num_l], --last, array_n);
                first = last;
            }
            if (num_r) {
                offsets_r += start_r;
                while (num_r--) array_swap(offsets_r_base - offsets_r[num_r], first, array_n), ++first;
                last = first;
            }
        }

        // Put the pivot in the right place.
        Iter pivot_pos = first - 1;
        array_move(&*begin, &*pivot_pos, array_n);
        array_move(&*pivot_pos, &*pivot, array_n);

        return std::make_pair(pivot_pos, already_partitioned);
    }

    // Similar function to the one above, except elements equal to the pivot are put to the left of
    // the pivot and it doesn't check or return if the passed sequence already was partitioned.
    // Since this is rarely used (the many equal case), and in that case pdqsort already has O(n)
    // performance, no block quicksort is applied here for simplicity.
    template<class Iter, class Compare>
    inline Iter partition_left(Iter begin, Iter end, size_t array_n, Compare comp) {
        typedef typename std::iterator_traits<Iter>::value_type T;

        T pivot[array_n];
        array_move(pivot, &*begin, array_n);
        Iter first = begin;
        Iter last = end;

        while (comp(pivot, &*--last, array_n));

        if (last + 1 == end) while (first < last && !comp(pivot, &*++first, array_n));
        else                 while (                !comp(pivot, &*++first, array_n));

        while (first < last) {
            array_swap(first, last, array_n);
            while (comp(pivot, &*--last, array_n));
            while (!comp(pivot, &*++first, array_n));
        }

        Iter pivot_pos = last;
        array_move(&*begin, &*pivot_pos, array_n);
        array_move(&*pivot_pos, pivot, array_n);

        return pivot_pos;
    }


    template<class Iter, class Compare>
    inline void pdqsort_loop(Iter begin, Iter end, size_t array_n, Compare comp,
                             int bad_allowed, bool leftmost = true) {
        typedef typename std::iterator_traits<Iter>::difference_type diff_t;

        // Use a while loop for tail recursion elimination.
        while (true) {
            diff_t size = end - begin;

            // Insertion sort is faster for small arrays.
            if (size < insertion_sort_threshold) {
                if (leftmost) insertion_sort(begin, end, array_n, comp);
                else unguarded_insertion_sort(begin, end, array_n, comp);
                return;
            }

            // Choose pivot as median of 3 or pseudomedian of 9.
            diff_t s2 = size / 2;
            if (size > ninther_threshold) {
                sort3(begin, begin + s2, end - 1, array_n, comp);
                sort3(begin + 1, begin + (s2 - 1), end - 2, array_n, comp);
                sort3(begin + 2, begin + (s2 + 1), end - 3, array_n, comp);
                sort3(begin + (s2 - 1), begin + s2, begin + (s2 + 1), array_n, comp);
                array_swap(begin, begin + s2, array_n);
            } else sort3(begin + s2, begin, end - 1, array_n, comp);

            // If *(begin - 1) is the end of the right partition of a previous partition operation
            // there is no element in [begin, end) that is smaller than *(begin - 1). Then if our
            // pivot compares equal to *(begin - 1) we change strategy, putting equal elements in
            // the left partition, greater elements in the right partition. We do not have to
            // recurse on the left partition, since it's sorted (all equal).
            if (!leftmost && !comp(&*(begin - 1), &*begin, array_n)) {
                begin = partition_left(begin, end, array_n, comp) + 1;
                continue;
            }

            // Partition and get results.
            std::pair<Iter, bool> part_result = partition_right_branchless(begin, end, array_n, comp);
            Iter pivot_pos = part_result.first;
            bool already_partitioned = part_result.second;

            // Check for a highly unbalanced partition.
            diff_t l_size = pivot_pos - begin;
            diff_t r_size = end - (pivot_pos + 1);
            bool highly_unbalanced = l_size < size / 8 || r_size < size / 8;

            // If we got a highly unbalanced partition we shuffle elements to break many patterns.
            if (highly_unbalanced) {
                if (l_size >= insertion_sort_threshold) {
                    array_swap(begin,             begin + l_size / 4, array_n);
                    array_swap(pivot_pos - 1, pivot_pos - l_size / 4, array_n);

                    if (l_size > ninther_threshold) {
                        array_swap(begin + 1,         begin + (l_size / 4 + 1), array_n);
                        array_swap(begin + 2,         begin + (l_size / 4 + 2), array_n);
                        array_swap(pivot_pos - 2, pivot_pos - (l_size / 4 + 1), array_n);
                        array_swap(pivot_pos - 3, pivot_pos - (l_size / 4 + 2), array_n);
                    }
                }

                if (r_size >= insertion_sort_threshold) {
                    array_swap(pivot_pos + 1, pivot_pos + (1 + r_size / 4), array_n);
                    array_swap(end - 1,                   end - r_size / 4, array_n);

                    if (r_size > ninther_threshold) {
                        array_swap(pivot_pos + 2, pivot_pos + (2 + r_size / 4), array_n);
                        array_swap(pivot_pos + 3, pivot_pos + (3 + r_size / 4), array_n);
                        array_swap(end - 2,             end - (1 + r_size / 4), array_n);
                        array_swap(end - 3,             end - (2 + r_size / 4), array_n);
                    }
                }
            } else {
                // If we were decently balanced and we tried to sort an already partitioned
                // sequence try to use insertion sort.
                if (already_partitioned &&
                    partial_insertion_sort(begin, pivot_pos, array_n, comp) &&
                    partial_insertion_sort(pivot_pos + 1, end, array_n, comp))
                    return;
            }

            // Sort the left partition first using recursion and do tail recursion elimination for
            // the right-hand partition.
            pdqsort_loop<Iter, Compare>(begin, pivot_pos, array_n, comp, bad_allowed, leftmost);
            begin = pivot_pos + 1;
            leftmost = false;
        }
    }

    template<class Iter, class Compare, size_t array_n>
    inline void pdqsort_loop_constant_size(Iter begin, Iter end, Compare comp,
                                           int bad_allowed, bool leftmost = true) {
        typedef typename std::iterator_traits<Iter>::difference_type diff_t;
        while (true) {
            diff_t size = end - begin;

            // Insertion sort is faster for small arrays.
            if (size < insertion_sort_threshold) {
                if (leftmost) insertion_sort(begin, end, array_n, comp);
                else unguarded_insertion_sort(begin, end, array_n, comp);
                return;
            }

            // Choose pivot as median of 3 or pseudomedian of 9.
            diff_t s2 = size / 2;
            if (size > ninther_threshold) {
                sort3(begin, begin + s2, end - 1, array_n, comp);
                sort3(begin + 1, begin + (s2 - 1), end - 2, array_n, comp);
                sort3(begin + 2, begin + (s2 + 1), end - 3, array_n, comp);
                sort3(begin + (s2 - 1), begin + s2, begin + (s2 + 1), array_n, comp);
                array_swap(begin, begin + s2, array_n);
            } else sort3(begin + s2, begin, end - 1, array_n, comp);

            // If *(begin - 1) is the end of the right partition of a previous partition operation
            // there is no element in [begin, end) that is smaller than *(begin - 1). Then if our
            // pivot compares equal to *(begin - 1) we change strategy, putting equal elements in
            // the left partition, greater elements in the right partition. We do not have to
            // recurse on the left partition, since it's sorted (all equal).
            if (!leftmost && !comp(&*(begin - 1), &*begin, array_n)) {
                begin = partition_left(begin, end, array_n, comp) + 1;
                continue;
            }

            // Partition and get results.
            std::pair<Iter, bool> part_result = partition_right_branchless(begin, end, array_n, comp);
            Iter pivot_pos = part_result.first;
            bool already_partitioned = part_result.second;

            // Check for a highly unbalanced partition.
            diff_t l_size = pivot_pos - begin;
            diff_t r_size = end - (pivot_pos + 1);
            bool highly_unbalanced = l_size < size / 8 || r_size < size / 8;

            // If we got a highly unbalanced partition we shuffle elements to break many patterns.
            if (highly_unbalanced) {
                if (l_size >= insertion_sort_threshold) {
                    array_swap(begin,             begin + l_size / 4, array_n);
                    array_swap(pivot_pos - 1, pivot_pos - l_size / 4, array_n);

                    if (l_size > ninther_threshold) {
                        array_swap(begin + 1,         begin + (l_size / 4 + 1), array_n);
                        array_swap(begin + 2,         begin + (l_size / 4 + 2), array_n);
                        array_swap(pivot_pos - 2, pivot_pos - (l_size / 4 + 1), array_n);
                        array_swap(pivot_pos - 3, pivot_pos - (l_size / 4 + 2), array_n);
                    }
                }

                if (r_size >= insertion_sort_threshold) {
                    array_swap(pivot_pos + 1, pivot_pos + (1 + r_size / 4), array_n);
                    array_swap(end - 1,                   end - r_size / 4, array_n);

                    if (r_size > ninther_threshold) {
                        array_swap(pivot_pos + 2, pivot_pos + (2 + r_size / 4), array_n);
                        array_swap(pivot_pos + 3, pivot_pos + (3 + r_size / 4), array_n);
                        array_swap(end - 2,             end - (1 + r_size / 4), array_n);
                        array_swap(end - 3,             end - (2 + r_size / 4), array_n);
                    }
                }
            } else {
                // If we were decently balanced and we tried to sort an already partitioned
                // sequence try to use insertion sort.
                if (already_partitioned &&
                    partial_insertion_sort(begin, pivot_pos, array_n, comp) &&
                    partial_insertion_sort(pivot_pos + 1, end, array_n, comp))
                    return;
            }

            // Sort the left partition first using recursion and do tail recursion elimination for
            // the right-hand partition.
            pdqsort_loop_constant_size<Iter, Compare, array_n>(begin, pivot_pos, comp, bad_allowed, leftmost);
            begin = pivot_pos + 1;
            leftmost = false;
        }
    }

    template<typename Iter>
    class array_iterator {
      public:
        typedef typename std::iterator_traits<Iter>::difference_type difference_type;
        typedef typename std::iterator_traits<Iter>::pointer pointer;
        typedef typename std::iterator_traits<Iter>::reference reference;
        typedef typename std::iterator_traits<Iter>::value_type value_type;
        typedef std::random_access_iterator_tag iterator_category;

      private:
        Iter wrapped_;
        size_t el_sz_;

      public:
        array_iterator(Iter wrapped, size_t el_sz)
                : wrapped_(std::move(wrapped)), el_sz_(el_sz) { }

        size_t size() const { return el_sz_; }

        reference operator*() const { return *wrapped_; }
        reference operator[](difference_type n) const { return *(*this + n); }

        array_iterator &operator++() {
            wrapped_ += el_sz_;
            return *this;
        }
        array_iterator &operator--() {
            wrapped_ -= el_sz_;
            return *this;
        }
        array_iterator operator++(int) {
            array_iterator res = *this;
            wrapped_ += el_sz_;
            return res;
        }
        array_iterator operator--(int) {
            array_iterator res = *this;
            wrapped_ -= el_sz_;
            return res;
        }

        array_iterator operator+(const difference_type &n) const {
            return array_iterator(wrapped_ + n * el_sz_, el_sz_);
        }
        array_iterator operator-(const difference_type &n) const {
            return array_iterator(wrapped_ - n * el_sz_, el_sz_);
        }

        array_iterator &operator+=(const difference_type &n) {
            wrapped_ += n * el_sz_;
            return *this;
        }
        array_iterator &operator-=(const difference_type &n) {
            wrapped_ -= n * el_sz_;
            return *this;
        }

        friend bool operator==(const array_iterator &r1, const array_iterator &r2) {
            return r1.wrapped_ == r2.wrapped_;
        }

        friend bool operator!=(const array_iterator &r1, const array_iterator &r2) {
            return r1.wrapped_ != r2.wrapped_;
        }

        friend bool operator<(const array_iterator &r1, const array_iterator &r2) {
            return r1.wrapped_ < r2.wrapped_;
        }

        friend bool operator<=(const array_iterator &r1, const array_iterator &r2) {
            return r1.wrapped_ <= r2.wrapped_;
        }

        friend bool operator>(const array_iterator &r1, const array_iterator &r2) {
            return r1.wrapped_ > r2.wrapped_;
        }

        friend bool operator>=(const array_iterator &r1, const array_iterator &r2) {
            return r1.wrapped_ >= r2.wrapped_;
        }


        friend array_iterator
        operator+(difference_type n,
                  const array_iterator &r2) {
            return r2 + n;
        }

        friend difference_type
        operator-(const array_iterator &r1,
                  const array_iterator &r2) {
            return (r1.wrapped_ - r2.wrapped_) / r1.el_sz_;
        }
    };

    template<typename Iter, size_t el_sz>
    class array_iterator_constant_size {
      public:
        typedef typename std::iterator_traits<Iter>::difference_type difference_type;
        typedef typename std::iterator_traits<Iter>::pointer pointer;
        typedef typename std::iterator_traits<Iter>::reference reference;
        typedef typename std::iterator_traits<Iter>::value_type value_type;
        typedef std::random_access_iterator_tag iterator_category;

      private:
        Iter wrapped_;

      public:
        array_iterator_constant_size(Iter wrapped)
                : wrapped_(std::move(wrapped)) { }

        size_t size() const { return el_sz; }

        reference operator*() const { return *wrapped_; }
        reference operator[](difference_type n) const { return *(*this + n); }

        array_iterator_constant_size &operator++() {
            wrapped_ += el_sz;
            return *this;
        }
        array_iterator_constant_size &operator--() {
            wrapped_ -= el_sz;
            return *this;
        }
        array_iterator_constant_size operator++(int) {
            array_iterator_constant_size res = *this;
            wrapped_ += el_sz;
            return res;
        }
        array_iterator_constant_size operator--(int) {
            array_iterator_constant_size res = *this;
            wrapped_ -= el_sz;
            return res;
        }

        array_iterator_constant_size operator+(const difference_type &n) const {
            return array_iterator_constant_size(wrapped_ + n * el_sz);
        }
        array_iterator_constant_size operator-(const difference_type &n) const {
            return array_iterator_constant_size(wrapped_ - n * el_sz);
        }

        array_iterator_constant_size &operator+=(const difference_type &n) {
            wrapped_ += n * el_sz;
            return *this;
        }
        array_iterator_constant_size &operator-=(const difference_type &n) {
            wrapped_ -= n * el_sz;
            return *this;
        }

        friend bool operator==(const array_iterator_constant_size &r1, const array_iterator_constant_size &r2) {
            return r1.wrapped_ == r2.wrapped_;
        }

        friend bool operator!=(const array_iterator_constant_size &r1, const array_iterator_constant_size &r2) {
            return r1.wrapped_ != r2.wrapped_;
        }

        friend bool operator<(const array_iterator_constant_size &r1, const array_iterator_constant_size &r2) {
            return r1.wrapped_ < r2.wrapped_;
        }

        friend bool operator<=(const array_iterator_constant_size &r1, const array_iterator_constant_size &r2) {
            return r1.wrapped_ <= r2.wrapped_;
        }

        friend bool operator>(const array_iterator_constant_size &r1, const array_iterator_constant_size &r2) {
            return r1.wrapped_ > r2.wrapped_;
        }

        friend bool operator>=(const array_iterator_constant_size &r1, const array_iterator_constant_size &r2) {
            return r1.wrapped_ >= r2.wrapped_;
        }


        friend array_iterator_constant_size
        operator+(difference_type n,
                  const array_iterator_constant_size &r2) {
            return r2 + n;
        }

        friend difference_type
        operator-(const array_iterator_constant_size &r1,
                  const array_iterator_constant_size &r2) {
            return (r1.wrapped_ - r2.wrapped_) / el_sz;
        }
    };

    template<class Ptr>
    struct array_compare {
        bool operator()(const Ptr this_ptr, const Ptr that_ptr, size_t size) const {
            for (size_t i = 0; i < size; ++i) {
                if (this_ptr[i] != that_ptr[i])
                    return this_ptr[i] < that_ptr[i];
            }

            return false;
        }
    };

    template<class Iter, size_t array_n>
    inline void pdqsort_const_array(Iter begin, Iter end) {
        using WrappedIter = array_iterator_constant_size<Iter, array_n>;
        using WrappedCompare = array_compare<typename std::iterator_traits<Iter>::pointer>;

        return pdqsort_loop_constant_size<WrappedIter, WrappedCompare, array_n> (
            WrappedIter(begin), WrappedIter(end),
            WrappedCompare(),
            pdqsort_pod_detail::log2(end - begin));
    }
}


template<class Iter>
inline void pdqsort_pod(Iter begin, Iter end, size_t array_n) {
    if (begin == end) return;

    using namespace pdqsort_pod_detail;
    using WrappedIter = array_iterator<Iter>;
    using WrappedCompare = array_compare<typename std::iterator_traits<Iter>::pointer>;

    switch (array_n) {
        case 1:
            return pdqsort_const_array<Iter, 1>(begin, end);
        case 2:
            return pdqsort_const_array<Iter, 2>(begin, end);
        case 3:
            return pdqsort_const_array<Iter, 3>(begin, end);
        case 4:
            return pdqsort_const_array<Iter, 4>(begin, end);
        default:
            break;
    }

    return pdqsort_loop<WrappedIter, WrappedCompare> (
        WrappedIter(begin, array_n), WrappedIter(end, array_n), array_n,
        WrappedCompare(),
        pdqsort_pod_detail::log2(end - begin));
}
