#ifndef __LIBCXX_SORT_HPP__
#define __LIBCXX_SORT_HPP__

#include <type_traits>
#include <cstring>
#include <utility>
#include <memory>
#include <iterator>
#include <cstddef>

namespace libcxx {

#ifndef __has_feature
#define __has_feature(x) 0
#endif

#if __has_feature(has_trivial_destructor) || (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3)

template <class _Tp> struct is_trivially_destructible
    : public std::integral_constant<bool, __has_trivial_destructor(_Tp)> {};

#else

template <class _Tp> struct __libcpp_trivial_destructor
    : public std::integral_constant<bool, std::is_scalar<_Tp>::value ||
                                          std::is_reference<_Tp>::value> {};

template <class _Tp> struct is_trivially_destructible
    : public __libcpp_trivial_destructor<typename std::remove_all_extents<_Tp>::type> {};

#endif

template <class _T1, class _T2 = _T1>
struct __less
{
    __attribute__ ((__always_inline__)) bool operator()(const _T1& __x, const _T1& __y) const {return __x < __y;}
    __attribute__ ((__always_inline__)) bool operator()(const _T1& __x, const _T2& __y) const {return __x < __y;}
    __attribute__ ((__always_inline__)) bool operator()(const _T2& __x, const _T1& __y) const {return __x < __y;}
    __attribute__ ((__always_inline__)) bool operator()(const _T2& __x, const _T2& __y) const {return __x < __y;}
};

template <class _T1>
struct __less<_T1, _T1>
{
    __attribute__ ((__always_inline__)) bool operator()(const _T1& __x, const _T1& __y) const {return __x < __y;}
};

template <class _T1>
struct __less<const _T1, _T1>
{
    __attribute__ ((__always_inline__)) bool operator()(const _T1& __x, const _T1& __y) const {return __x < __y;}
};

template <class _T1>
struct __less<_T1, const _T1>
{
    __attribute__ ((__always_inline__)) bool operator()(const _T1& __x, const _T1& __y) const {return __x < __y;}
};

struct __destruct_n
{
private:
    size_t size;

    template <class _Tp>
    __attribute__ ((__always_inline__)) void __process(_Tp* __p, std::false_type) throw()
        {for (size_t __i = 0; __i < size; ++__i, ++__p) __p->~_Tp();}

    template <class _Tp>
    __attribute__ ((__always_inline__)) void __process(_Tp*, std::true_type) throw()
        {}

    __attribute__ ((__always_inline__)) void __incr(std::false_type) throw()
        {++size;}
    __attribute__ ((__always_inline__)) void __incr(std::true_type) throw()
        {}

    __attribute__ ((__always_inline__)) void __set(size_t __s, std::false_type) throw()
        {size = __s;}
    __attribute__ ((__always_inline__)) void __set(size_t, std::true_type) throw()
        {}
public:
    __attribute__ ((__always_inline__)) explicit __destruct_n(size_t __s) throw()
        : size(__s) {}

    template <class _Tp>
    __attribute__ ((__always_inline__)) void __incr(_Tp*) throw()
        {__incr(std::integral_constant<bool, is_trivially_destructible<_Tp>::value>());}

    template <class _Tp>
    __attribute__ ((__always_inline__)) void __set(size_t __s, _Tp*) throw()
        {__set(__s, std::integral_constant<bool, is_trivially_destructible<_Tp>::value>());}

    template <class _Tp>
    __attribute__ ((__always_inline__)) void operator()(_Tp* __p) throw()
        {__process(__p, std::integral_constant<bool, is_trivially_destructible<_Tp>::value>());}
};

// Make sure we can find custom swap via ADL.
using std::swap;

// stable, 2-3 compares, 0-2 swaps
template <class _Compare, class _ForwardIterator>
unsigned
__sort3(_ForwardIterator __x, _ForwardIterator __y, _ForwardIterator __z, _Compare __c)
{
    unsigned __r = 0;
    if (!__c(*__y, *__x))          // if x <= y
    {
        if (!__c(*__z, *__y))      // if y <= z
            return __r;            // x <= y && y <= z
                                   // x <= y && y > z
        swap(*__y, *__z);          // x <= z && y < z
        __r = 1;
        if (__c(*__y, *__x))       // if x > y
        {
            swap(*__x, *__y);      // x < y && y <= z
            __r = 2;
        }
        return __r;                // x <= y && y < z
    }
    if (__c(*__z, *__y))           // x > y, if y > z
    {
        swap(*__x, *__z);          // x < y && y < z
        __r = 1;
        return __r;
    }
    swap(*__x, *__y);              // x > y && y <= z
    __r = 1;                       // x < y && x <= z
    if (__c(*__z, *__y))           // if y > z
    {
        swap(*__y, *__z);          // x <= y && y < z
        __r = 2;
    }
    return __r;
}                                  // x <= y && y <= z

// stable, 3-6 compares, 0-5 swaps

template <class _Compare, class _ForwardIterator>
unsigned
__sort4(_ForwardIterator __x1, _ForwardIterator __x2, _ForwardIterator __x3,
            _ForwardIterator __x4, _Compare __c)
{
    unsigned __r = libcxx::__sort3<_Compare>(__x1, __x2, __x3, __c);
    if (__c(*__x4, *__x3))
    {
        swap(*__x3, *__x4);
        ++__r;
        if (__c(*__x3, *__x2))
        {
            swap(*__x2, *__x3);
            ++__r;
            if (__c(*__x2, *__x1))
            {
                swap(*__x1, *__x2);
                ++__r;
            }
        }
    }
    return __r;
}

// stable, 4-10 compares, 0-9 swaps

template <class _Compare, class _ForwardIterator>
unsigned
__sort5(_ForwardIterator __x1, _ForwardIterator __x2, _ForwardIterator __x3,
            _ForwardIterator __x4, _ForwardIterator __x5, _Compare __c)
{
    unsigned __r = libcxx::__sort4<_Compare>(__x1, __x2, __x3, __x4, __c);
    if (__c(*__x5, *__x4))
    {
        swap(*__x4, *__x5);
        ++__r;
        if (__c(*__x4, *__x3))
        {
            swap(*__x3, *__x4);
            ++__r;
            if (__c(*__x3, *__x2))
            {
                swap(*__x2, *__x3);
                ++__r;
                if (__c(*__x2, *__x1))
                {
                    swap(*__x1, *__x2);
                    ++__r;
                }
            }
        }
    }
    return __r;
}

// Assumes size > 0
template <class _Compare, class _BirdirectionalIterator>
void
__selection_sort(_BirdirectionalIterator __first, _BirdirectionalIterator __last, _Compare __comp)
{
    _BirdirectionalIterator __lm1 = __last;
    for (--__lm1; __first != __lm1; ++__first)
    {
        _BirdirectionalIterator __i = std::min_element<_BirdirectionalIterator,
                                                       typename std::add_lvalue_reference<_Compare>::type>
                                                       (__first, __last, __comp);
        if (__i != __first)
            swap(*__first, *__i);
    }
}

template <class _Compare, class _BirdirectionalIterator>
void
__insertion_sort(_BirdirectionalIterator __first, _BirdirectionalIterator __last, _Compare __comp)
{
    typedef typename std::iterator_traits<_BirdirectionalIterator>::value_type value_type;
    if (__first != __last)
    {
        _BirdirectionalIterator __i = __first;
        for (++__i; __i != __last; ++__i)
        {
            _BirdirectionalIterator __j = __i;
            value_type __t(std::move(*__j));
            for (_BirdirectionalIterator __k = __i; __k != __first && __comp(__t,  *--__k); --__j)
                *__j = std::move(*__k);
            *__j = std::move(__t);
        }
    }
}

template <class _Compare, class _RandomAccessIterator>
void
__insertion_sort_3(_RandomAccessIterator __first, _RandomAccessIterator __last, _Compare __comp)
{
    typedef typename std::iterator_traits<_RandomAccessIterator>::value_type value_type;
    _RandomAccessIterator __j = __first+2;
    libcxx::__sort3<_Compare>(__first, __first+1, __j, __comp);
    for (_RandomAccessIterator __i = __j+1; __i != __last; ++__i)
    {
        if (__comp(*__i, *__j))
        {
            value_type __t(std::move(*__i));
            _RandomAccessIterator __k = __j;
            __j = __i;
            do
            {
                *__j = std::move(*__k);
                __j = __k;
            } while (__j != __first && __comp(__t, *--__k));
            *__j = std::move(__t);
        }
        __j = __i;
    }
}

template <class _Compare, class _RandomAccessIterator>
bool
__insertion_sort_incomplete(_RandomAccessIterator __first, _RandomAccessIterator __last, _Compare __comp)
{
    switch (__last - __first)
    {
    case 0:
    case 1:
        return true;
    case 2:
        if (__comp(*--__last, *__first))
            swap(*__first, *__last);
        return true;
    case 3:
        libcxx::__sort3<_Compare>(__first, __first+1, --__last, __comp);
        return true;
    case 4:
        libcxx::__sort4<_Compare>(__first, __first+1, __first+2, --__last, __comp);
        return true;
    case 5:
        libcxx::__sort5<_Compare>(__first, __first+1, __first+2, __first+3, --__last, __comp);
        return true;
    }
    typedef typename std::iterator_traits<_RandomAccessIterator>::value_type value_type;
    _RandomAccessIterator __j = __first+2;
    libcxx::__sort3<_Compare>(__first, __first+1, __j, __comp);
    const unsigned __limit = 8;
    unsigned __count = 0;
    for (_RandomAccessIterator __i = __j+1; __i != __last; ++__i)
    {
        if (__comp(*__i, *__j))
        {
            value_type __t(std::move(*__i));
            _RandomAccessIterator __k = __j;
            __j = __i;
            do
            {
                *__j = std::move(*__k);
                __j = __k;
            } while (__j != __first && __comp(__t, *--__k));
            *__j = std::move(__t);
            if (++__count == __limit)
                return ++__i == __last;
        }
        __j = __i;
    }
    return true;
}

template <class _Compare, class _BirdirectionalIterator>
void
__insertion_sort_move(_BirdirectionalIterator __first1, _BirdirectionalIterator __last1,
                      typename std::iterator_traits<_BirdirectionalIterator>::value_type* __first2, _Compare __comp)
{
    typedef typename std::iterator_traits<_BirdirectionalIterator>::value_type value_type;
    if (__first1 != __last1)
    {
        __destruct_n __d(0);
        std::unique_ptr<value_type, __destruct_n&> __h(__first2, __d);
        value_type* __last2 = __first2;
        ::new(__last2) value_type(std::move(*__first1));
        __d.__incr((value_type*)0);
        for (++__last2; ++__first1 != __last1; ++__last2)
        {
            value_type* __j2 = __last2;
            value_type* __i2 = __j2;
            if (__comp(*__first1, *--__i2))
            {
                ::new(__j2) value_type(std::move(*__i2));
                __d.__incr((value_type*)0);
                for (--__j2; __i2 != __first2 && __comp(*__first1,  *--__i2); --__j2)
                    *__j2 = std::move(*__i2);
                *__j2 = std::move(*__first1);
            }
            else
            {
                ::new(__j2) value_type(std::move(*__first1));
                __d.__incr((value_type*)0);
            }
        }
        __h.release();
    }
}

template <class _Compare, class _RandomAccessIterator>
void
__sort(_RandomAccessIterator __first, _RandomAccessIterator __last, _Compare __comp)
{
    // _Compare is known to be a reference type
    typedef typename std::iterator_traits<_RandomAccessIterator>::difference_type difference_type;
    typedef typename std::iterator_traits<_RandomAccessIterator>::value_type value_type;
#if 0
    const difference_type __limit = std::is_trivially_copy_constructible<value_type>::value &&
                                    std::is_trivially_copy_assignable<value_type>::value ? 30 : 6;
#else
    const difference_type __limit = 10;
#endif

    while (true)
    {
    __restart:
        difference_type __len = __last - __first;
        switch (__len)
        {
        case 0:
        case 1:
            return;
        case 2:
            if (__comp(*--__last, *__first))
                swap(*__first, *__last);
            return;
        case 3:
            libcxx::__sort3<_Compare>(__first, __first+1, --__last, __comp);
            return;
        case 4:
            libcxx::__sort4<_Compare>(__first, __first+1, __first+2, --__last, __comp);
            return;
        case 5:
            libcxx::__sort5<_Compare>(__first, __first+1, __first+2, __first+3, --__last, __comp);
            return;
        }
        if (__len <= __limit)
        {
            libcxx::__insertion_sort_3<_Compare>(__first, __last, __comp);
            return;
        }
        // __len > 5
        _RandomAccessIterator __m = __first;
        _RandomAccessIterator __lm1 = __last;
        --__lm1;
        unsigned __n_swaps;
        {
        difference_type __delta;
        if (__len >= 1000)
        {
            __delta = __len/2;
            __m += __delta;
            __delta /= 2;
            __n_swaps = libcxx::__sort5<_Compare>(__first, __first + __delta, __m, __m+__delta, __lm1, __comp);
        }
        else
        {
            __delta = __len/2;
            __m += __delta;
            __n_swaps = libcxx::__sort3<_Compare>(__first, __m, __lm1, __comp);
        }
        }
        // *__m is median
        // partition [__first, __m) < *__m and *__m <= [__m, __last)
        // (this inhibits tossing elements equivalent to __m around unnecessarily)
        _RandomAccessIterator __i = __first;
        _RandomAccessIterator __j = __lm1;
        // j points beyond range to be tested, *__m is known to be <= *__lm1
        // The search going up is known to be guarded but the search coming down isn't.
        // Prime the downward search with a guard.
        if (!__comp(*__i, *__m))  // if *__first == *__m
        {
            // *__first == *__m, *__first doesn't go in first part
            // manually guard downward moving __j against __i
            while (true)
            {
                if (__i == --__j)
                {
                    // *__first == *__m, *__m <= all other elements
                    // Parition instead into [__first, __i) == *__first and *__first < [__i, __last)
                    ++__i;  // __first + 1
                    __j = __last;
                    if (!__comp(*__first, *--__j))  // we need a guard if *__first == *(__last-1)
                    {
                        while (true)
                        {
                            if (__i == __j)
                                return;  // [__first, __last) all equivalent elements
                            if (__comp(*__first, *__i))
                            {
                                swap(*__i, *__j);
                                ++__n_swaps;
                                ++__i;
                                break;
                            }
                            ++__i;
                        }
                    }
                    // [__first, __i) == *__first and *__first < [__j, __last) and __j == __last - 1
                    if (__i == __j)
                        return;
                    while (true)
                    {
                        while (!__comp(*__first, *__i))
                            ++__i;
                        while (__comp(*__first, *--__j))
                            ;
                        if (__i >= __j)
                            break;
                        swap(*__i, *__j);
                        ++__n_swaps;
                        ++__i;
                    }
                    // [__first, __i) == *__first and *__first < [__i, __last)
                    // The first part is sorted, sort the secod part
                    // libcxx::__sort<_Compare>(__i, __last, __comp);
                    __first = __i;
                    goto __restart;
                }
                if (__comp(*__j, *__m))
                {
                    swap(*__i, *__j);
                    ++__n_swaps;
                    break;  // found guard for downward moving __j, now use unguarded partition
                }
            }
        }
        // It is known that *__i < *__m
        ++__i;
        // j points beyond range to be tested, *__m is known to be <= *__lm1
        // if not yet partitioned...
        if (__i < __j)
        {
            // known that *(__i - 1) < *__m
            // known that __i <= __m
            while (true)
            {
                // __m still guards upward moving __i
                while (__comp(*__i, *__m))
                    ++__i;
                // It is now known that a guard exists for downward moving __j
                while (!__comp(*--__j, *__m))
                    ;
                if (__i > __j)
                    break;
                swap(*__i, *__j);
                ++__n_swaps;
                // It is known that __m != __j
                // If __m just moved, follow it
                if (__m == __i)
                    __m = __j;
                ++__i;
            }
        }
        // [__first, __i) < *__m and *__m <= [__i, __last)
        if (__i != __m && __comp(*__m, *__i))
        {
            swap(*__i, *__m);
            ++__n_swaps;
        }
        // [__first, __i) < *__i and *__i <= [__i+1, __last)
        // If we were given a perfect partition, see if insertion sort is quick...
        if (__n_swaps == 0)
        {
            bool __fs = libcxx::__insertion_sort_incomplete<_Compare>(__first, __i, __comp);
            if (libcxx::__insertion_sort_incomplete<_Compare>(__i+1, __last, __comp))
            {
                if (__fs)
                    return;
                __last = __i;
                continue;
            }
            else
            {
                if (__fs)
                {
                    __first = ++__i;
                    continue;
                }
            }
        }
        // sort smaller range with recursive call and larger with tail recursion elimination
        if (__i - __first < __last - __i)
        {
            libcxx::__sort<_Compare>(__first, __i, __comp);
            // libcxx::__sort<_Compare>(__i+1, __last, __comp);
            __first = ++__i;
        }
        else
        {
            libcxx::__sort<_Compare>(__i+1, __last, __comp);
            // libcxx::__sort<_Compare>(__first, __i, __comp);
            __last = __i;
        }
    }
}

// This forwarder keeps the top call and the recursive calls using the same instantiation, forcing a reference _Compare
template <class _RandomAccessIterator, class _Compare>
inline __attribute__ ((__always_inline__))
void
sort(_RandomAccessIterator __first, _RandomAccessIterator __last, _Compare __comp)
{
    typedef typename std::add_lvalue_reference<_Compare>::type _Comp_ref;
    libcxx::__sort<_Comp_ref>(__first, __last, __comp);
}

template <class _RandomAccessIterator>
inline __attribute__ ((__always_inline__))
void
sort(_RandomAccessIterator __first, _RandomAccessIterator __last)
{
    libcxx::sort(__first, __last, __less<typename std::iterator_traits<_RandomAccessIterator>::value_type>());
}
};

#endif // __LIBCXX_SORT_HPP__
