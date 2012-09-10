#ifndef __HAMMER_POINTER_ITERATOR_HPP__
#define __HAMMER_POINTER_ITERATOR_HPP__

#include <iterator>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

template<typename T>
class pointer_iterator : public std::iterator<std::random_access_iterator_tag, T> {
protected:
  T *data_;
  size_t stride_;

public:
  typedef std::random_access_iterator_tag iterator_category;
  typedef typename std::iterator<std::random_access_iterator_tag, T>::value_type value_type;
  typedef typename std::iterator<std::random_access_iterator_tag, T>::difference_type difference_type;
  typedef typename std::iterator<std::random_access_iterator_tag, T>::reference reference;
  typedef typename std::iterator<std::random_access_iterator_tag, T>::pointer   pointer;

  pointer_iterator() : data_(NULL), stride_(1) {}

  template<typename T2>
  pointer_iterator(const pointer_iterator<T2> &r) : data_(&(*r)), stride_(r.stride_) {}

  pointer_iterator(pointer data, size_t stride = 1) : data_(data), stride_(stride) {}

  template<typename T2>
  pointer_iterator& operator=(const pointer_iterator<T2> &r) {
    data_ = &(*r);
    stride_ = r.stride_;
    return *this;
  }

  pointer_iterator& operator++() {
    data_ += stride_;
    return *this;
  }

  pointer_iterator& operator--() {
    data_ -= stride_;
    return *this;
  }

  pointer_iterator operator++(int) {
    pointer_iterator res = *this;
    data_ += stride_;

    return res;
  }

  pointer_iterator operator--(int) {
    pointer_iterator res = *this;
    data_ -= stride_;

    return res;
  }

  pointer_iterator operator+(const difference_type &n) const {
    return pointer_iterator(data_ + n*stride_);
  }

  pointer_iterator& operator+=(const difference_type &n) {
    data_ += n*stride_; return *this;
  }

  pointer_iterator operator-(const difference_type &n) const {
    return pointer_iterator(pointer(data_ - n*stride_));
  }

  pointer_iterator& operator-=(const difference_type &n) {
    data_ -= n*stride_; return *this;
  }

  reference operator*() const {
    return *data_;
  }

  pointer operator->() const {
    return data_;
  }

  reference operator[](const difference_type &n) const {
    return data_[n*stride_];
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
  return (r1.data_ == r2.data_ &&
          r1.stride_ == r2.stride_);
}

template<typename T>
inline bool operator!=(const pointer_iterator<T> &r1,
                       const pointer_iterator<T> &r2) {
  return (r1.data_ != r2.data_ ||
          r1.stride_ != r2.stride_);
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
  return (r1.data_ - r2.data_) / r1.stride_;
}

#endif // __HAMMER_POINTER_ITERATOR_HPP__
