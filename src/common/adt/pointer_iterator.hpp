//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_POINTER_ITERATOR_HPP__
#define __HAMMER_POINTER_ITERATOR_HPP__

#include <iterator>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace adt {

template<typename T>
class pointer_iterator : public std::iterator<std::random_access_iterator_tag, T> {
protected:
    T *data_;

public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef typename std::iterator<std::random_access_iterator_tag, T>::value_type value_type;
    typedef typename std::iterator<std::random_access_iterator_tag, T>::difference_type difference_type;
    typedef typename std::iterator<std::random_access_iterator_tag, T>::reference reference;
    typedef typename std::iterator<std::random_access_iterator_tag, T>::pointer pointer;

    pointer_iterator() : data_(NULL) { }

    template<typename T2>
    pointer_iterator(const pointer_iterator<T2> &r) : data_(&(*r)) { }

    pointer_iterator(pointer data) : data_(data) { }

    template<typename T2>
    pointer_iterator &operator=(const pointer_iterator<T2> &r) {
        data_ = &(*r);
        return *this;
    }

    pointer_iterator &operator++() {
        data_ += 1;
        return *this;
    }

    pointer_iterator &operator--() {
        data_ -= 1;
        return *this;
    }

    pointer_iterator operator++(int) {
        pointer_iterator res = *this;
        data_ += 1;

        return res;
    }

    pointer_iterator operator--(int) {
        pointer_iterator res = *this;
        data_ -= 1;

        return res;
    }

    pointer_iterator operator+(const difference_type &n) const {
        return pointer_iterator(data_ + n);
    }

    pointer_iterator &operator+=(const difference_type &n) {
        data_ += n;
        return *this;
    }

    pointer_iterator operator-(const difference_type &n) const {
        return pointer_iterator(pointer(data_ - n));
    }

    pointer_iterator &operator-=(const difference_type &n) {
        data_ -= n;
        return *this;
    }

    reference operator*() const {
        return *data_;
    }

    pointer operator->() const {
        return data_;
    }

    reference operator[](const difference_type &n) const {
        return data_[n];
    }

    template<typename T2>
    friend bool operator==(const pointer_iterator<T2> &r1,
                           const pointer_iterator<T2> &r2);

    template<typename T2>
    friend bool operator!=(const pointer_iterator<T2> &r1,
                           const pointer_iterator<T2> &r2);

    template<typename T2>
    friend bool operator<(const pointer_iterator<T2> &r1,
                          const pointer_iterator<T2> &r2);

    template<typename T2>
    friend bool operator>(const pointer_iterator<T2> &r1,
                          const pointer_iterator<T2> &r2);

    template<typename T2>
    friend bool operator<=(const pointer_iterator<T2> &r1,
                           const pointer_iterator<T2> &r2);

    template<typename T2>
    friend bool operator>=(const pointer_iterator<T2> &r1,
                           const pointer_iterator<T2> &r2);

    template<typename T2>
    friend typename pointer_iterator<T2>::difference_type
            operator+(const pointer_iterator<T2> &r1,
                      const pointer_iterator<T2> &r2);

    template<typename T2>
    friend typename pointer_iterator<T2>::difference_type
            operator-(const pointer_iterator<T2> &r1,
                      const pointer_iterator<T2> &r2);
};

template<typename T>
inline bool operator==(const pointer_iterator<T> &r1,
                       const pointer_iterator<T> &r2) {
    return (r1.data_ == r2.data_);
}

template<typename T>
inline bool operator!=(const pointer_iterator<T> &r1,
                       const pointer_iterator<T> &r2) {
    return (r1.data_ != r2.data_);
}

template<typename T>
inline bool operator<(const pointer_iterator<T> &r1,
                      const pointer_iterator<T> &r2) {
    return (r1.data_ < r2.data_);
}

template<typename T>
inline bool operator>(const pointer_iterator<T> &r1,
                      const pointer_iterator<T> &r2) {
    return (r1.data_ > r2.data_);
}

template<typename T>
inline bool operator<=(const pointer_iterator<T> &r1,
                       const pointer_iterator<T> &r2) {
    return (r1.data_ <= r2.data_);
}

template<typename T>
inline bool operator>=(const pointer_iterator<T> &r1,
                       const pointer_iterator<T> &r2) {
    return (r1.data_ >= r2.data_);
}

template<typename T>
inline typename pointer_iterator<T>::difference_type
operator-(const pointer_iterator<T> &r1,
          const pointer_iterator<T> &r2) {
    return (r1.data_ - r2.data_);
}

} //adt

#endif // __HAMMER_POINTER_ITERATOR_HPP__
