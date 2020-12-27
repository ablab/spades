#ifndef __ADT_SMALL_POD_VECTOR__
#define __ADT_SMALL_POD_VECTOR__

#pragma once

#include <llvm/ADT/PointerIntPair.h>

#include <cstring>
#include <vector>
#include <type_traits>

namespace adt {

#define LIKELY(EXPR) __builtin_expect((bool)(EXPR), true)
#define UNLIKELY(EXPR) __builtin_expect((bool)(EXPR), false)

namespace impl {
// Class holding the data for SmallPOD vector. It is expected that descendants
// could implement different memory allocation policies.
template<class T>
struct SmallPODVectorData {
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T value_type;
    typedef T *iterator;

    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;

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

    static constexpr unsigned SmallSizeIntBits = 3;
    static constexpr unsigned MaxSmall = (1 << SmallSizeIntBits) - 1;

    typedef typename std::vector<T> vector_type;

    typedef llvm::PointerIntPair<void *, SmallSizeIntBits, size_t,
                                 PointerUnionTraits<T *, vector_type *> > container_type;

    container_type data_;

    SmallPODVectorData() {
        reset();
    }

    void reset() {
        data_.setPointer(nullptr);
        data_.setInt(0);
    }

    __attribute__((always_inline))
    bool empty() const {
        return data_.getInt() == 0 && data_.getPointer() == nullptr;
    }

    __attribute__((always_inline))
    vector_type *vector() const {
        return (data_.getInt() == 0 ?
                static_cast<vector_type *>(data_.getPointer()) : nullptr);
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

    size_t capacity() const {
        const auto v = vector();
        if (UNLIKELY(v != nullptr))
            return v->capacity();

        return data_.getInt();
    }
};

// Allocate all SmallPODvector memory on heap
template<class T>
struct HeapAllocatedStorage : public SmallPODVectorData<T> {
    using typename SmallPODVectorData<T>::vector_type;
    using typename SmallPODVectorData<T>::size_type;
    using SmallPODVectorData<T>::MaxSmall;

    HeapAllocatedStorage() {
        this->reset();
    }

    void grow(size_type N) {
        void *data = this->data_.getPointer(), *new_data = data;
        size_t sz = this->data_.getInt(), new_sz = N;

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
            // Otherwise, simply change the size of the allocated space if we have to grow
            if (N > sz) {
                new_data = realloc(data, N * sizeof(T));
                new_sz = N;
            } else if (N == 0) {
                free(data);
                new_data = nullptr;
                new_sz = 0;
            }
        }

        this->data_.setPointer(new_data);
        this->data_.setInt(new_sz);
    }
};

// Store all data into pre-allocated storage. On overflow bail down to heap
// allocation
template<class T, unsigned PreAllocated = 3>
struct PreAllocatedStorage : public SmallPODVectorData<T> {
    using typename SmallPODVectorData<T>::vector_type;
    using typename SmallPODVectorData<T>::size_type;

    T buffer[PreAllocated];
    static_assert(PreAllocated <= SmallPODVectorData<T>::MaxSmall,
                  "Cannot fit so much pre-allocated data");

    using SmallPODVectorData<T>::MaxSmall;

    PreAllocatedStorage() {
        this->reset();
    }

    PreAllocatedStorage(PreAllocatedStorage &&that) {
        void *data = that.data_.getPointer(), *new_data = data;
        size_t sz = that.data_.getInt(), new_sz = sz;

        if (sz == 0 && data != nullptr) { // vector case
            // Do nothing special;
        } else if (sz) {
            new_data = buffer;
            memcpy(new_data, data, sz * sizeof(T));
        }


        this->data_.setPointer(new_data);
        this->data_.setInt(new_sz);

        that.reset();
    }

    void grow(size_type N) {
        void *data = this->data_.getPointer(), *new_data = data;
        size_t sz = this->data_.getInt(), new_sz = N;

        if (UNLIKELY(sz == 0 && data != nullptr)) { // vector case
            vector_type *v = static_cast<vector_type *>(data);
            if (N > PreAllocated) {
                v->resize(N);
                new_data = v;
                new_sz = 0;
            } else { // We have to turn vector into array
                if (N) {
                    new_data = buffer;
                    new_sz = N;
                    memcpy(new_data, v->data(), N * sizeof(T));
                } else {
                    new_data = nullptr;
                    new_sz = 0;
                }
                delete v;
            }
        } else if (UNLIKELY(N > PreAllocated)) {
            // Ok, we have to grow too much - allocate new vector
            vector_type *new_vector = new vector_type((T *) data, (T *) data + sz);
            new_vector->resize(N);
            new_data = new_vector;
            new_sz = 0;
        } else {
            // Otherwise, simply change the size of the allocated space if we have to grow
            if (N > sz) {
                new_data = buffer;
                new_sz = N;
            } else if (N == 0) {
                new_data = nullptr;
                new_sz = 0;
            }
        }

        this->data_.setPointer(new_data);
        this->data_.setInt(new_sz);
    }
};

// Hybrid allocation strategy between pre-allocated and heap allocated
template<class T, unsigned PreAllocated = 3>
struct HybridAllocatedStorage : public SmallPODVectorData<T> {
    using typename SmallPODVectorData<T>::vector_type;
    using typename SmallPODVectorData<T>::size_type;

    T buffer[PreAllocated];
    static_assert(PreAllocated <= SmallPODVectorData<T>::MaxSmall,
                  "Cannot fit so much pre-allocated data");

    using SmallPODVectorData<T>::MaxSmall;

    HybridAllocatedStorage() {
        this->reset();
    }

    HybridAllocatedStorage(const HybridAllocatedStorage &that) {
        // This will allocate memory and setup buffers as necessary
        grow(that.size());

        void *data = that.data_.getPointer(), *new_data = data;
        memcpy(new_data, data, this->size() * sizeof(T));
    }
    
    HybridAllocatedStorage(HybridAllocatedStorage &&that) {
        void *data = that.data_.getPointer(), *new_data = data;
        size_t sz = that.data_.getInt(), new_sz = sz;

        if (sz == 0 && data != nullptr) { // vector case
            // Do nothing special;
        } else if (sz > PreAllocated) { // heap allocated case
            // Do nothing special;
        } else if (sz) {
            new_data = buffer;
            memcpy(new_data, data, sz * sizeof(T));
        }

        this->data_.setPointer(new_data);
        this->data_.setInt(new_sz);

        that.reset();
    }

    void grow(size_type N) {
        void *data = this->data_.getPointer(), *new_data = data;
        size_t sz = this->data_.getInt(), new_sz = N;

        // Destructor
        if (sz == 0 && data == nullptr && N == 0)
            return;

        if (UNLIKELY(sz == 0 && data != nullptr)) { // vector case
            vector_type *v = static_cast<vector_type *>(data);
            if (N > MaxSmall) {
                v->resize(N);
                new_data = v;
                new_sz = 0;
            } else { // We have to turn vector into array
                if (N > PreAllocated) {
                    new_data = malloc(N * sizeof(T));
                    new_sz = N;
                    memcpy(new_data, v->data(), N * sizeof(T));
                } else if (N) {
                    new_data = buffer;
                    new_sz = N;
                    memcpy(new_data, v->data(), N * sizeof(T));
                } else {
                    new_data = nullptr;
                    new_sz = 0;
                }
                delete v;
            }
        } else if (UNLIKELY(N > MaxSmall)) { // Ok, we have to grow too much - allocate new vector
            vector_type *new_vector = new vector_type((T *) data, (T *) data + sz);
            new_vector->resize(N);
            new_data = new_vector;
            new_sz = 0;
            if (sz > PreAllocated)
                free(data);
        } else if (UNLIKELY(sz > PreAllocated)) { // Heap-allocated storage
            if (N > sz) { // Realloc if necessary
                new_data = realloc(data, N * sizeof(T));
                new_sz = N;
            } else if (N > 0 && N <= PreAllocated) { // Go to static buffer
                new_data = buffer;
                new_sz = N;
                memcpy(new_data, data, N * sizeof(T));
                free(data);
            } else if (N == 0) {
                free(data);
                new_data = nullptr;
                new_sz = 0;
            }
        } else { // Pre-allocated storage
            if (UNLIKELY(N > PreAllocated)) { // Turn to heap-allocated storage
                new_data = malloc(N * sizeof(T));
                new_sz = N;
                memcpy(new_data, data, sz * sizeof(T));
            } else if (N > sz) { // Otherwise, simply change the size of the allocated space if we have to grow
                new_data = buffer;
                new_sz = N;
            } else if (N == 0) {
                new_data = nullptr;
                new_sz = 0;
            }
        }

        this->data_.setPointer(new_data);
        this->data_.setInt(new_sz);
    }
};

}

template<class T, class Container = impl::HybridAllocatedStorage<T>>
class SmallPODVector {
    Container data_;
    typedef SmallPODVector<T> self;

public:
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T value_type;
    typedef T *iterator;
    typedef const T *const_iterator;

    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;

    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;

    static_assert(std::is_trivially_copyable<value_type>::value,
                  "Value type for SmallPODVector should be trivially copyable");

private:


public:
    SmallPODVector() = default;

    SmallPODVector(size_type size, const T &value = T()) {
        this->assign(size, value);
    }

    SmallPODVector(const self &that) {
        assign(that.begin(), that.end());
    }

    const self &operator=(const self &that) {
        // Avoid self-assignment.
        if (this == &that) return *this;
        assign(that.begin(), that.end());
        return *this;
    }

    SmallPODVector(self &&that)
            : data_(std::move(that.data_)) {
        that.data_.reset();
    }

    const self &operator=(const self &&that) {
        // Avoid self-assignment.
        if (this == &that) return *this;

        data_.grow(0);
        data_ = std::move(that.data_);
        that.data_.reset();

        return *this;
    }

    ~SmallPODVector() {
        data_.grow(0);
    }

    size_type max_size() const { return size_type(-1) / sizeof(T); }

    __attribute__((always_inline))
    bool empty() const { return data_.empty(); }

    __attribute__((always_inline))
    size_type size() const { return data_.size(); }

    __attribute__((always_inline))
    pointer data() { return data_.data(); }
    __attribute__((always_inline))
    const_pointer cdata() const { return data_.cdata(); }
    size_t capacity() const { return data_.capacity(); }

    // forward iterator creation methods.
    iterator begin() { return (iterator)(data()); }
    const_iterator begin() const { return (const_iterator)(cdata()); }
    const_iterator cbegin() const { return (const_iterator)(cdata()); }
    iterator end() { return (iterator)(data() + size()); }
    const_iterator end() const { return (const_iterator)(cdata() + size()); }
    const_iterator cend() const { return (const_iterator)(cdata() + size()); }

    // reverse iterator creation methods.
    reverse_iterator rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }

    __attribute__((always_inline))
    reference operator[](size_type idx) {
        return begin()[idx];
    }

    __attribute__((always_inline))
    const_reference operator[](size_type idx) const {
        return begin()[idx];
    }

    reference front() {
        return begin()[0];
    }

    const_reference front() const {
        return begin()[0];
    }

    reference back() {
        return end()[-1];
    }

    const_reference back() const {
        return end()[-1];
    }

    void push_back(const T &value) {
        const auto v = data_.vector();
        if (UNLIKELY(v != nullptr)) {
            v->push_back(value);
            return;
        }

        data_.grow(this->size() + 1);
        memcpy(this->end() - 1, &value, sizeof(T));
    }

    void pop_back() {
        // This will reallocate to array, if necessary.
        data_.grow(this->size() - 1);
    }

    T pop_back_val() {
        T res = ::std::move(this->back());
        this->pop_back();
        return res;
    }

    void clear() {
        data_.grow(0);
    }

    void resize(size_type count) {
        data_.grow(count);
        std::uninitialized_fill(this->begin() + count, this->end(), T());
    }

    void resize(size_type count, const T &value) {
        data_.grow(count);
        std::uninitialized_fill(this->begin() + count, this->end(), value);
    }

    void reserve(size_type count) {
        if (auto v = data_.vector()) {
            v->reserve(count);
        }
    }

    void assign(size_type count, const T &value) {
        data_.grow(count);
        std::uninitialized_fill(this->begin(), this->end(), value);
    }

    template<class InputIt>
    void assign(InputIt first, InputIt last) {
        data_.grow(last - first);
        std::uninitialized_copy(first, last, this->begin());
    }

    iterator erase(const_iterator pos) {
        size_type idx = pos - this->begin();
        std::copy(iterator(pos + 1), this->end(), iterator(pos));
        data_.grow(this->size() - 1); // This might invalidate iterators

        return this->begin() + idx;
    }

    iterator erase(const_iterator first, const_iterator last) {
        difference_type idx = first - this->begin();
        std::copy(iterator(last), this->end(), iterator(first));
        data_.grow(this->size() - (last - first)); // This might invalidate iterators

        return this->begin() + idx;
    }

    iterator insert(const_iterator pos, const T &value) {
        if (pos == this->end()) {
            this->push_back(value);
            return this->end() - 1;
        }

        difference_type idx = pos - this->begin();
        size_type sz = this->size();

        data_.grow(sz + 1); // This might invalidate iterators

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

} //adt
#endif // __ADT_SMALL_POD_VECTOR__
