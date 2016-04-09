#ifndef __ADT_SMALL_POD_VECTOR__
#define __ADT_SMALL_POD_VECTOR__

#pragma once

#include <llvm/PointerIntPair.h>

#include <vector>
#include <type_traits>

namespace adt {

#define LIKELY(EXPR) __builtin_expect((bool)(EXPR), true)
#define UNLIKELY(EXPR) __builtin_expect((bool)(EXPR), false)

template<class T>
class SmallPODVector {
    template<typename PT1, typename PT2>
    class PointerUnionTraits {
    public:
        static inline void *getAsVoidPointer(void *P) { return P; }

        static inline void *getFromVoidPointer(void *P) { return P; }

        enum {
            PT1BitsAv = (int) (llvm::PointerLikeTypeTraits<PT1>::NumLowBitsAvailable),
            PT2BitsAv = (int) (llvm::PointerLikeTypeTraits<PT2>::NumLowBitsAvailable),
            NumLowBitsAvailable = PT1BitsAv < PT2BitsAv ? PT1BitsAv : PT2BitsAv
        };
    };

    static const unsigned SmallSizeIntBits = 3;
    static const unsigned MaxSmall = (1 << SmallSizeIntBits) - 1;

    typedef typename std::vector<T> vector_type;

    typedef llvm::PointerIntPair<void *, SmallSizeIntBits, size_t,
            PointerUnionTraits<T *, vector_type *> > container_type;

    typedef SmallPODVector<T> self;
    container_type data_;

public:
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T value_type;
    typedef T *iterator;
    typedef const T *const_iterator;

    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;

    typedef T &reference;
    typedef const T &const_reference;
    typedef T *pointer;
    typedef const T *const_pointer;

// workaround missing "is_trivially_copyable" in g++ < 5.0
#if __GNUG__ && __GNUC__ < 5
#define IS_TRIVIALLY_COPYABLE(T) __has_trivial_copy(T)
#else
#define IS_TRIVIALLY_COPYABLE(T) std::is_trivially_copyable<T>::value
#endif

    static_assert(IS_TRIVIALLY_COPYABLE(value_type), "Value type for SmallPODVector should be trivially copyable");

#undef IS_TRIVIALLY_COPYABLE

private:
    vector_type *vector() const {
        return (data_.getInt() == 0 ? static_cast<vector_type *>(data_.getPointer()) : nullptr);
    }

    void impl_resize(size_type N) {
        void *data = data_.getPointer(), *new_data = data;
        size_t sz = data_.getInt(), new_sz = N;

        if (UNLIKELY(sz == 0 && data != nullptr)) { // vector case
            vector_type *v = static_cast<vector_type *>(data);
            if (N > MaxSmall) {
                v->resize(N);
                new_data = v;
                new_sz = 0;
            } else { // We have to turn vector into array
                if (N) {
                    new_data = malloc(N * sizeof(T));
                    new_sz = N;
                    memcpy(new_data, v->data(), N * sizeof(T));
                } else {
                    new_data = nullptr;
                    new_sz = 0;
                }
                delete v;
            }
        } else if (UNLIKELY(N > MaxSmall)) {
            // Ok, we have to grow too much - allocate new vector
            vector_type *new_vector = new vector_type((T *) data, (T *) data + sz);
            new_vector->resize(N);
            if (data)
                free(data);
            new_data = new_vector;
            new_sz = 0;
        } else {
            // Otherwise, simply change the size of the allocated space
            if (N) {
                new_data = realloc(data, N * sizeof(T));
                new_sz = N;
            } else {
                free(data);
                new_data = nullptr;
                new_sz = 0;
            }
        }

        data_.setPointer(new_data);
        data_.setInt(new_sz);
    }

public:
    SmallPODVector<T>() = default;

    SmallPODVector<T>(size_type size, const T &value = T()) {
        this->assign(size, value);
    }

    SmallPODVector<T>(const self &that) {
        assign(that.begin(), that.end());
    }

    const self &operator=(const self &that) {
        // Avoid self-assignment.
        if (this == &that) return *this;
        assign(that.begin(), that.end());
        return *this;
    }

    SmallPODVector<T>(self &&that) {
        data_ = that.data_;
        that.data_.setPointer(nullptr);
        that.data_.setInt(0);
    }

    const self &operator=(const self &&that) {
        // Avoid self-assignment.
        if (this == &that) return *this;

        this->impl_resize(0);
        data_ = that.data_;
        that.data_.setPointer(nullptr);
        that.data_.setInt(0);

        return *this;
    }

    ~SmallPODVector<T>() {
        this->impl_resize(0);
    }

    __attribute__((always_inline))
    bool empty() const {
        return data_.getInt() == 0 && data_.getPointer() == nullptr;
    }

    __attribute__((always_inline))
    size_type size() const {
        const auto v = vector();
        if (UNLIKELY(v != nullptr))
            return v->size();

        return data_.getInt();
    }

    __attribute__((always_inline))
    pointer data() {
        const auto v = vector();
        if (UNLIKELY(v != nullptr))
            return v->data();

        return pointer(data_.getPointer());
    }

    __attribute__((always_inline))
    const_pointer cdata() const {
        const auto v = vector();
        if (UNLIKELY(v != nullptr))
            return v->data();

        return const_pointer(data_.getPointer());
    }

    size_type max_size() const { return size_type(-1) / sizeof(T); }

    size_t capacity() const {
        const auto v = vector();
        if (UNLIKELY(v != nullptr))
            return v->capacity();

        return data_.getInt();
    }

    // forward iterator creation methods.
    __attribute__((always_inline))
    iterator begin() {
        return (iterator)(data());
    }

    __attribute__((always_inline))
    const_iterator begin() const {
        return (const_iterator)(cdata());
    }

    __attribute__((always_inline))
    const_iterator cbegin() const {
        return (const_iterator)(cdata());
    }

    __attribute__((always_inline))
    iterator end() {
        return (iterator)(data() + size());
    }

    __attribute__((always_inline))
    const_iterator end() const {
        return (const_iterator)(cdata() + size());
    }

    __attribute__((always_inline))
    const_iterator cend() const {
        return (const_iterator)(cdata() + size());
    }

    // reverse iterator creation methods.
    reverse_iterator rbegin() { return reverse_iterator(end()); }

    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }

    reverse_iterator rend() { return reverse_iterator(begin()); }

    const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }

    __attribute__((always_inline))
    reference operator[](size_type idx) {
        assert(idx < size());
        return begin()[idx];
    }

    __attribute__((always_inline))
    const_reference operator[](size_type idx) const {
        assert(idx < size());
        return begin()[idx];
    }

    reference front() {
        assert(!empty());
        return begin()[0];
    }

    const_reference front() const {
        assert(!empty());
        return begin()[0];
    }

    reference back() {
        assert(!empty());
        return end()[-1];
    }

    const_reference back() const {
        assert(!empty());
        return end()[-1];
    }

    void push_back(const T &value) {
        const auto v = vector();
        if (UNLIKELY(v != nullptr)) {
            v->push_back(value);
            return;
        }

        this->impl_resize(this->size() + 1);
        memcpy(this->end() - 1, &value, sizeof(T));
    }

    void pop_back() {
        // This will reallocate to array, if necessary.
        this->impl_resize(this->size() - 1);
    }

    T pop_back_val() {
        T res = ::std::move(this->back());
        this->pop_back();
        return res;
    }

    void clear() {
        this->impl_resize(0);
    }

    void resize(size_type count) {
        this->impl_resize(count);
        std::uninitialized_fill(this->begin() + count, this->end(), T());
    }

    void resize(size_type count, const T &value) {
        this->impl_resize(count);
        std::uninitialized_fill(this->begin() + count, this->end(), value);
    }

    void reserve(size_type count) {
        if (auto v = vector()) {
            v->reserve(count);
        }
    }

    void assign(size_type count, const T &value) {
        this->impl_resize(count);
        std::uninitialized_fill(this->begin(), this->end(), value);
    }

    template<class InputIt>
    void assign(InputIt first, InputIt last) {
        this->impl_resize(last - first);
        std::uninitialized_copy(first, last, this->begin());
    }

    iterator erase(const_iterator pos) {
        size_type idx = pos - this->begin();
        std::copy(iterator(pos + 1), this->end(), iterator(pos));
        this->impl_resize(this->size() - 1); // This might invalidate iterators

        return this->begin() + idx;
    }

    iterator erase(const_iterator first, const_iterator last) {
        difference_type idx = first - this->begin();
        std::copy(iterator(last), this->end(), iterator(first));
        this->impl_resize(this->size() - (last - first)); // This might invalidate iterators

        return this->begin() + idx;
    }

    iterator insert(const_iterator pos, const T &value) {
        if (pos == this->end()) {
            this->push_back(value);
            return this->end() - 1;
        }

        difference_type idx = pos - this->begin();
        size_type sz = this->size();

        this->impl_resize(sz + 1); // This might invalidate iterators

        iterator it = this->begin() + idx;
        std::copy_backward(it, this->end() - 1, this->end());

        // If we just moved the element we're inserting, be sure to update the
        // reference.
        const T *vptr = &value;
        if (it <= vptr && vptr < this->end())
            ++vptr;

        *it = *vptr;

        return it;
    }

    template<typename... ArgTypes>
    void emplace_back(ArgTypes &&... args) {
        value_type tmp(std::forward<ArgTypes>(args)...);
        push_back(std::move(tmp));
    }

    template<typename... ArgTypes>
    iterator emplace(const_iterator pos, ArgTypes &&... args) {
        value_type tmp(std::forward<ArgTypes>(args)...);
        return insert(pos, std::move(tmp));
    }

    bool operator==(const self &rhs) const {
        if (this->size() != rhs.size()) return false;
        return std::equal(this->begin(), this->end(), rhs.begin());
    }

    bool operator!=(const self &rhs) const {
        return !(*this == rhs);
    }

    bool operator<(const self &rhs) const {
        return std::lexicographical_compare(this->begin(), this->end(),
                                            rhs.begin(), rhs.end());
    }
};

#undef LIKELY
#undef UNLIKELY

}

#endif // __ADT_SMALL_POD_VECTOR__
