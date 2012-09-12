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

public:
  typedef std::random_access_iterator_tag iterator_category;
  typedef typename std::iterator<std::random_access_iterator_tag, T>::value_type value_type;
  typedef typename std::iterator<std::random_access_iterator_tag, T>::difference_type difference_type;
  typedef typename std::iterator<std::random_access_iterator_tag, T>::reference reference;
  typedef typename std::iterator<std::random_access_iterator_tag, T>::pointer   pointer;

  pointer_iterator() : data_(NULL) {}

  template<typename T2>
  pointer_iterator(const pointer_iterator<T2> &r) : data_(&(*r)) {}

  pointer_iterator(pointer data, size_t stride = 1) : data_(data) {}

  template<typename T2>
  pointer_iterator& operator=(const pointer_iterator<T2> &r) {
    data_ = &(*r);
    return *this;
  }

  pointer_iterator& operator++() {
    data_ += 1;
    return *this;
  }

  pointer_iterator& operator--() {
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

  pointer_iterator& operator+=(const difference_type &n) {
    data_ += n; return *this;
  }

  pointer_iterator operator-(const difference_type &n) const {
    return pointer_iterator(pointer(data_ - n));
  }

  pointer_iterator& operator-=(const difference_type &n) {
    data_ -= n; return *this;
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


template<typename T>
struct array_ref {
  T *ptr;
  size_t size;

  array_ref() = delete;
  array_ref(T *p, size_t sz) : ptr(p), size(sz) {  }

  operator T*() {
    return ptr;
  }

  array_ref& operator=(const array_ref &that) {
    if (this != &that)
      memcpy(ptr, that.ptr, size*sizeof(T));

    return *this;
  }

  bool operator==(const array_ref &that) const {
    return (ptr == that.ptr && size == that.size);
  }
  bool operator!=(const array_ref &that) const {
    return !operator==(that);
  }

  bool operator<(const array_ref &that) const {
    return 0 > memcmp(ptr, that.ptr, that.size*sizeof(T));
  }

  struct equal_to : std::binary_function<array_ref, array_ref, bool> {
    bool operator()(const array_ref &x, const array_ref &y) {
      return (0 == memcmp(x.ptr, y.ptr, x.size*sizeof(T)));
    }
  };
};

template<typename T>
void swap(array_ref<T> &lhs, array_ref<T> &rhs) {
  std::swap_ranges(lhs.ptr, lhs.ptr + lhs.size, rhs.ptr);
}

template<typename T>
class pointer_array_iterator : public std::iterator<std::random_access_iterator_tag, T> {
 private:
  array_ref<T> data_;

 public:
  typedef std::random_access_iterator_tag iterator_category;
  typedef typename std::iterator<std::random_access_iterator_tag, T>::difference_type difference_type;

  typedef array_ref<T>  value_type;
  typedef array_ref<T>& reference;
  typedef array_ref<T>* pointer;

  pointer_array_iterator() : data_(NULL, 0) {}

  template<typename T2>
  pointer_array_iterator(const pointer_array_iterator<T2> &r) : data_(r.data_.ptr, r.data_.size) {}

  pointer_array_iterator(T *data, size_t size = 1) : data_(data, size) { }

  pointer_array_iterator& operator=(const pointer_array_iterator<T> &r) {
    data_.ptr = r.data_.ptr;
    data_.size = r.data_.size;

    return *this;
  }

  pointer_array_iterator& operator++() {
    data_.ptr += data_.size;
    return *this;
  }

  pointer_array_iterator& operator--() {
    data_.ptr -= data_.size;
    return *this;
  }

  pointer_array_iterator operator++(int) {
    pointer_array_iterator res = *this;
    data_.ptr += data_.size;
    return res;
  }

  pointer_array_iterator operator--(int) {
    pointer_array_iterator res = *this;
    data_.ptr -= data_.size_;
    return res;
  }

  pointer_array_iterator operator+(const difference_type &n) const {
    return pointer_array_iterator(data_.ptr + n*data_.size, data_.size);
  }

  pointer_array_iterator& operator+=(const difference_type &n) {
    data_.ptr += n*data_.size;
    return *this;
  }

  pointer_array_iterator operator-(const difference_type &n) const {
    return pointer_array_iterator(data_.ptr - n*data_.size, data_.size);
  }

  pointer_array_iterator& operator-=(const difference_type &n) {
    data_.ptr -= n*data_.size;
    return *this;
  }

  reference operator*() {
    return data_;
  }

  const reference operator*() const {
    return data_;
  }

  template<typename T2>
  friend bool operator==(const pointer_array_iterator<T2> &r1,
                         const pointer_array_iterator<T2> &r2);

  template<typename T2>
  friend bool operator!=(const pointer_array_iterator<T2> &r1,
                         const pointer_array_iterator<T2> &r2);

  template<typename T2>
  friend bool operator<(const pointer_array_iterator<T2> &r1,
                        const pointer_array_iterator<T2> &r2);

  template<typename T2>
  friend bool operator>(const pointer_array_iterator<T2> &r1,
                        const pointer_array_iterator<T2> &r2);

  template<typename T2>
  friend bool operator<=(const pointer_array_iterator<T2> &r1,
                         const pointer_array_iterator<T2> &r2);

  template<typename T2>
  friend bool operator>=(const pointer_array_iterator<T2> &r1,
                         const pointer_array_iterator<T2> &r2);

  template<typename T2>
  friend typename pointer_array_iterator<T2>::difference_type
  operator+(const pointer_array_iterator<T2> &r1,
            const pointer_array_iterator<T2> &r2);

  template<typename T2>
  friend typename pointer_array_iterator<T2>::difference_type
  operator-(const pointer_array_iterator<T2> &r1,
            const pointer_array_iterator<T2> &r2);
};

template<typename T>
inline bool operator==(const pointer_array_iterator<T> &r1,
                       const pointer_array_iterator<T> &r2) {
  return (r1.data_.ptr == r2.data_.ptr);
}

template<typename T>
inline bool operator!=(const pointer_array_iterator<T> &r1,
                       const pointer_array_iterator<T> &r2) {
  return (r1.data_.ptr != r2.data_.ptr);
}

template<typename T>
inline bool operator<(const pointer_array_iterator<T> &r1,
                      const pointer_array_iterator<T> &r2) {
  return (r1.data_.ptr < r2.data_.ptr);
}

template<typename T>
inline bool operator>(const pointer_array_iterator<T> &r1,
                      const pointer_array_iterator<T> &r2) {
  return (r1.data_.ptr > r2.data_.ptr);
}

template<typename T>
inline bool operator<=(const pointer_array_iterator<T> &r1,
                       const pointer_array_iterator<T> &r2) {
  return (r1.data_.ptr <= r2.data_.ptr);
}

template<typename T>
inline bool operator>=(const pointer_array_iterator<T> &r1,
                       const pointer_array_iterator<T> &r2) {
  return (r1.data_.ptr >= r2.data_.ptr);
}

template<typename T>
inline typename pointer_array_iterator<T>::difference_type
operator-(const pointer_array_iterator<T> &r1,
          const pointer_array_iterator<T> &r2) {
  size_t sz = (r1.data_.ptr - r2.data_.ptr) / r1.data_.size;
  return sz;
}

#endif // __HAMMER_POINTER_ITERATOR_HPP__
