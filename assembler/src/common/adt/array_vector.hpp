//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __ARRAY_VECTOR_HPP__
#define __ARRAY_VECTOR_HPP__

#include <algorithm>
#include <memory>

#include <cstdlib>
#include <cstring>
#include <cstddef>

namespace adt {

template<class _Cp, bool _IsConst>
class array_vector_iterator;

template<class _Cp>
class array_reference;

template<class _Cp>
class array_const_reference;

template<typename ElTy>
struct array_equal_to;


template<class _Cp>
class array {
    typedef typename _Cp::storage_type storage_type;
    typedef typename _Cp::storage_pointer storage_pointer;
    typedef typename _Cp::const_storage_pointer const_storage_pointer;
    typedef typename _Cp::size_type size_type;

#if defined(__clang__)
    friend typename _Cp::self;
#else

    friend class _Cp::self;

#endif

    friend class array_vector_iterator<_Cp, false>;

    friend class array_reference<_Cp>;

    friend class array_const_reference<_Cp>;

    storage_pointer ptr_;
    size_type size_;
    bool allocated;

public:
    ~array() {
        if (allocated)
            delete[] ptr_;
    }

    size_t size() const {
        return size_;
    }

    size_t data_size() const {
        return size_ * sizeof(storage_type);
    }

    storage_pointer data() const {
        return ptr_;
    }

    array(const array &that) {
        size_ = that.size_;
        ptr_ = new storage_type[size_];
        allocated = true;
        memcpy(ptr_, that.ptr_, data_size());
    }

    array(const array_reference<_Cp> that) {
        size_ = that.size();
        ptr_ = new storage_type[size_];
        allocated = true;
        memcpy(ptr_, that.data(), data_size());
    }

    array &operator=(const array &that) {
        storage_pointer this_ptr = data(), that_ptr = that.data();
        if (this_ptr != that_ptr)
            memcpy(this_ptr, that_ptr, data_size());

        return *this;
    }

    array &operator=(const array_reference<_Cp> that) {
        storage_pointer this_ptr = data(), that_ptr = that.data();
        if (this_ptr != that_ptr)
            memcpy(this_ptr, that_ptr, data_size());

        return *this;
    }

    array &operator=(const_storage_pointer that_ptr) {
        storage_pointer this_ptr = data();
        if (this_ptr != that_ptr)
            memcpy(this_ptr, that_ptr, data_size());

        return *this;
    }

    bool operator<(const array &that) const {
        storage_pointer this_ptr = data(), that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return this_ptr[i] < that_ptr[i];
        }

        return false;
    }

    bool operator<(const array_reference<_Cp> that) const {
        storage_pointer this_ptr = data(), that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return this_ptr[i] < that_ptr[i];
        }

        return false;
    }

    bool operator==(const array &that) const {
        storage_pointer this_ptr = data(), that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return false;
        }

        return true;
    }

    bool operator==(const array_reference<_Cp> that) const {
        storage_pointer this_ptr = data(), that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return false;
        }

        return true;
    }

    bool operator!=(const array &that) const {
        return !operator==(that);
    }

    bool operator!=(const array_reference<_Cp> that) const {
        return !operator==(that);
    }

private:
    array(storage_pointer p, size_type sz) :
            ptr_(p), size_(sz), allocated(false) { }
};


template<class _Cp>
class array_reference {
    typedef typename _Cp::storage_type storage_type;
    typedef typename _Cp::storage_pointer storage_pointer;
    typedef typename _Cp::const_storage_pointer const_storage_pointer;
    typedef typename _Cp::size_type size_type;

#if defined(__clang__)
    friend typename _Cp::self;
#else

    friend class _Cp::self;

#endif

    friend class array_vector_iterator<_Cp, false>;

    friend class array<_Cp>;

    friend struct adt::array_equal_to<storage_type>;

    storage_pointer ptr_;
    size_type size_;

public:
    size_t size() const {
        return size_;
    }

    size_t data_size() const {
        return size() * sizeof(storage_type);
    }

    storage_pointer data() const {
        return ptr_;
    }

    array_reference &operator=(const array<_Cp> &that) {
        storage_pointer this_ptr = data(), that_ptr = that.data();
        if (this_ptr != that_ptr)
            memcpy(this_ptr, that_ptr, data_size());

        return *this;
    }

    array_reference &operator=(const_storage_pointer that_ptr) {
        storage_pointer this_ptr = data();
        if (this_ptr != that_ptr)
            memcpy(this_ptr, that_ptr, data_size());

        return *this;
    }

    array_reference &operator=(const array_reference that) {
        storage_pointer this_ptr = data(), that_ptr = that.data();
        if (this_ptr != that_ptr)
            memcpy(this_ptr, that_ptr, data_size());

        return *this;
    }

    bool operator<(const array<_Cp> &that) const {
        storage_pointer this_ptr = data(), that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return this_ptr[i] < that_ptr[i];
        }

        return false;
    }

    bool operator<(const array_reference that) const {
        storage_pointer this_ptr = data(), that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return this_ptr[i] < that_ptr[i];
        }

        return false;
    }

    bool operator==(const array<_Cp> &that) const {
        storage_pointer this_ptr = data(), that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return false;
        }

        return true;
    }

    bool operator==(const array_reference that) const {
        storage_pointer this_ptr = data(), that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return false;
        }

        return true;
    }

    bool operator!=(const array_reference that) const {
        return !operator==(that);
    }

    bool operator!=(const array<_Cp> &that) const {
        return !operator==(that);
    }

private:
    array_reference(storage_pointer p, size_type sz) :
            ptr_(p), size_(sz) { }
};

template<class _Cp>
class array_const_reference {
    typedef typename _Cp::storage_type storage_type;
    typedef typename _Cp::storage_pointer storage_pointer;
    typedef typename _Cp::const_storage_pointer const_storage_pointer;
    typedef typename _Cp::size_type size_type;

#if defined(__clang__)
    friend typename _Cp::self;
#else

    friend class _Cp::self;

#endif

    friend class array_vector_iterator<_Cp, true>;

    friend struct adt::array_equal_to<storage_type>;

    const_storage_pointer ptr_;
    size_type size_;

public:
    size_t size() const {
        return size_;
    }

    size_t data_size() const {
        return size() * sizeof(storage_type);
    }

    const_storage_pointer data() const {
        return ptr_;
    }

    array_const_reference(const array_const_reference &that)
            : ptr_(that.ptr_), size_(that.size_) { }

    bool operator<(array_const_reference that) const {
        const storage_pointer this_ptr = data(), that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return this_ptr[i] < that_ptr[i];
        }

        return false;
    }

    bool operator==(array_const_reference that) const {
        const storage_pointer this_ptr = data(), that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return false;
        }

        return true;
    }

    bool operator==(const array_reference<_Cp> that) const {
        const_storage_pointer this_ptr = data();
        const storage_pointer that_ptr = that.data();

        for (size_t i = 0; i < size(); ++i) {
            if (this_ptr[i] != that_ptr[i])
                return false;
        }

        return true;
    }

    bool operator!=(const array_const_reference that) const {
        return !operator==(that);
    }

    bool operator!=(const array_reference<_Cp> that) const {
        return !operator==(that);
    }

private:
    array_const_reference(const_storage_pointer p, size_type sz) :
            ptr_(p), size_(sz) { }

    array_const_reference &operator=(const array_const_reference &that);
};

}
// This is hack. Never do this again!
#ifdef __GLIBCXX__
namespace std {
    template<typename _Cp>
    struct __are_same<adt::array_reference<_Cp>, adt::array<_Cp> &> {
        enum {
            __value = 1
        };
        typedef __true_type __type;
    };

    template<typename _Cp>
    struct __are_same<adt::array<_Cp> &, adt::array_reference<_Cp> > {
        enum {
            __value = 1
        };
        typedef __true_type __type;
    };
}
#endif
namespace adt {

template<typename _Cp>
void swap(array_reference<_Cp> lhs, array_reference<_Cp> rhs) {
    std::swap_ranges(lhs.data(), lhs.data() + lhs.size(), rhs.data());
}

template<typename _Cp>
void swap(array<_Cp> &lhs, array_reference<_Cp> rhs) {
    std::swap_ranges(lhs.data(), lhs.data() + lhs.size(), rhs.data());
}

template<typename _Cp>
void swap(array_reference<_Cp> lhs, array<_Cp> &rhs) {
    std::swap_ranges(lhs.data(), lhs.data() + lhs.size(), rhs.data());
}

template<typename _Cp, bool _IsConst>
class array_vector_iterator {
public:
    typedef typename _Cp::difference_type difference_type;
    typedef array_vector_iterator pointer;
    typedef typename std::conditional<_IsConst, array_const_reference<_Cp>, array_reference<_Cp> >::type reference;
    typedef array<_Cp> value_type;

    typedef std::random_access_iterator_tag iterator_category;

private:
    typedef typename _Cp::storage_type storage_type;
    typedef typename _Cp::storage_pointer storage_pointer;
    typedef typename _Cp::size_type size_type;

#if defined(__clang__)
  friend typename _Cp::self;
#else

    friend class _Cp::self;

#endif

    storage_pointer data_;
    size_type el_sz_;

public:
    array_vector_iterator(storage_pointer data, size_type el_sz)
            : data_(data), el_sz_(el_sz) { }

    size_t size() const {
        return el_sz_;
    }

    size_t data_size() const {
        return el_sz_ * sizeof(storage_type);
    }

    storage_pointer data() const {
        return data_;
    }

    reference operator*() const {
        return reference(data_, el_sz_);
    }

    reference operator[](difference_type n) const {
        return *(*this + n);
    }

    array_vector_iterator &operator++() {
        data_ += el_sz_;
        return *this;
    }

    array_vector_iterator &operator--() {
        data_ -= el_sz_;
        return *this;
    }

    array_vector_iterator operator++(int) {
        array_vector_iterator res = *this;
        data_ += el_sz_;
        return res;
    }

    array_vector_iterator operator--(int) {
        array_vector_iterator res = *this;
        data_ -= el_sz_;
        return res;
    }

    array_vector_iterator operator+(const difference_type &n) const {
        return array_vector_iterator(data_ + n * el_sz_, el_sz_);
    }

    array_vector_iterator &operator+=(const difference_type &n) {
        data_ += n * el_sz_;
        return *this;
    }

    array_vector_iterator operator-(const difference_type &n) const {
        return array_vector_iterator(data_ - n * el_sz_, el_sz_);
    }

    array_vector_iterator &operator-=(const difference_type &n) {
        data_ -= n * el_sz_;
        return *this;
    }

    friend bool operator==(const array_vector_iterator &r1,
                           const array_vector_iterator &r2) {
        return r1.data_ == r2.data_;
    }

    friend bool operator!=(const array_vector_iterator &r1,
                           const array_vector_iterator &r2) {
        return r1.data_ != r2.data_;
    }

    friend bool operator<(const array_vector_iterator &r1,
                          const array_vector_iterator &r2) {
        return r1.data_ < r2.data_;
    }

    friend bool operator<=(const array_vector_iterator &r1,
                           const array_vector_iterator &r2) {
        return r1.data_ <= r2.data_;
    }

    friend bool operator>(const array_vector_iterator &r1,
                          const array_vector_iterator &r2) {
        return r1.data_ > r2.data_;
    }

    friend bool operator>=(const array_vector_iterator &r1,
                           const array_vector_iterator &r2) {
        return r1.data_ >= r2.data_;
    }


    friend array_vector_iterator
    operator+(difference_type n,
              const array_vector_iterator &r2) {
        return r2 + n;
    }

    friend difference_type
    operator-(const array_vector_iterator &r1,
              const array_vector_iterator &r2) {
        return (r1.data_ - r2.data_) / r1.el_sz_;
    }
};

template<typename ElTy>
class array_vector {
public:
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    typedef array_reference<array_vector> reference;
    typedef array_const_reference<array_vector> const_reference;
    typedef array<array_vector> value_type;
    typedef array_vector_iterator<array_vector, false> iterator;
    typedef array_vector_iterator<array_vector, true> const_iterator;

private:
    typedef ElTy storage_type;
    typedef array_vector self;
    typedef storage_type *storage_pointer;
    typedef const storage_type *const_storage_pointer;

    friend class array<self>;

    friend class array_reference<self>;

    friend class array_const_reference<self>;

    friend class array_vector_iterator<self, true>;

    friend class array_vector_iterator<self, false>;

    storage_pointer data_;
    size_type size_;
    size_type el_sz_;

public:
    array_vector(storage_pointer data, size_type sz, size_type el_sz)
            : data_(data), size_(sz), el_sz_(el_sz) { }

    reference operator[](size_t pos) {
        return reference(data_ + pos * el_sz_, el_sz_);
    }

    const ElTy *operator[](size_t pos) const {
        return data_ + pos * el_sz_;
    }

    iterator begin() {
        return iterator(data_, el_sz_);
    }

    iterator end() {
        return iterator(data_ + size_ * el_sz_, el_sz_);
    }

    const_iterator begin() const {
        return const_iterator(data_, el_sz_);
    }

    const_iterator end() const {
        return const_iterator(data_ + size_ * el_sz_, el_sz_);
    }

    const_iterator cbegin() const {
        return const_iterator(data_, el_sz_);
    }

    const_iterator cend() const {
        return const_iterator(data_ + size_ * el_sz_, el_sz_);
    }

    size_t size() const { return size_; }

    storage_pointer data() const { return data_; }

    void set_size(size_t size) {
        size_ = size;
    }

    void set_data(storage_pointer data) {
        data_ = data;
    }
};

template<typename ElTy>
struct array_less {
    typedef typename array_vector<ElTy>::value_type value;
    typedef typename array_vector<ElTy>::reference reference;

    bool operator()(const value &lhs, const value &rhs) const {
        return lhs < rhs;
    }

    bool operator()(const value &lhs, const reference rhs) const {
        return lhs < rhs;
    }

    bool operator()(const reference lhs, const value &rhs) const {
        return lhs < rhs;
    }

    bool operator()(const reference lhs, const reference rhs) const {
        return lhs < rhs;
    }
};

template<typename ElTy>
struct array_equal_to {
    typedef typename array_vector<ElTy>::value_type value;
    typedef typename array_vector<ElTy>::reference reference;
    typedef typename array_vector<ElTy>::const_reference const_reference;

    bool operator()(const value &lhs, const value &rhs) const {
        return lhs == rhs;
    }

    bool operator()(const value &lhs, const reference rhs) const {
        return lhs == rhs;
    }

    bool operator()(const reference lhs, const value &rhs) const {
        return lhs == rhs;
    }

    bool operator()(const reference lhs, const ElTy *rhs, size_t sz) const {
        return lhs == reference(rhs, sz);
    }

    bool operator()(const reference lhs, const reference rhs) const {
        return lhs == rhs;
    }

    bool operator()(const ElTy *lhs, size_t sz, const reference rhs) const {
        return const_reference(lhs, sz) == rhs;
    }
};

} //adt
#endif
