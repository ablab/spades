//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace omnigraph {
namespace de {

template<class T>
class StrongWeakPtr {
  public:
    typedef T  element_type;
    typedef T* pointer;

    StrongWeakPtr() noexcept
            : ptr_(pointer(), false) {}

    StrongWeakPtr(std::nullptr_t) noexcept
            : ptr_(pointer(), false) {}

    StrongWeakPtr(pointer p, bool owning = true) noexcept
            : ptr_(std::move(p), owning) { }

    StrongWeakPtr(StrongWeakPtr &&p) noexcept
            : ptr_(p.release(), p.owning()) {}

    StrongWeakPtr& operator=(StrongWeakPtr &&p) noexcept {
        reset(p.release(), p.owning());
        return *this;
    }

    ~StrongWeakPtr() {
        reset();
    }

    StrongWeakPtr &operator=(std::nullptr_t) noexcept {
        reset();
        return *this;
    }

    typename std::add_lvalue_reference<T>::type operator*() const {
        return *ptr_.getPointer();
    }

    pointer operator->() const noexcept {
        return ptr_.getPointer();
    }

    pointer get() const noexcept {
        return ptr_.getPointer();
    }
    
    explicit operator bool() const noexcept {
        return ptr_.getPointer() != nullptr;
    }

    bool owning() const noexcept {
        return ptr_.getInt();
    }
    
    pointer release() noexcept {
        pointer p = ptr_.getPointer();
        ptr_ = raw_type();
        return p;
    }
    
    void reset(pointer p = pointer(), bool own = true) {
        pointer tmp = ptr_.getPointer(); bool is_owning = ptr_.getInt();
        ptr_ = raw_type(p, own);
        if (is_owning)
            delete tmp;
    }

    void swap(StrongWeakPtr &p) noexcept {
        std::swap(p.ptr_, ptr_);
    }
    
  private:
    llvm::PointerIntPair<pointer, 1, bool> ptr_;
  public:
    typedef decltype(ptr_) raw_type;
};


template<class T>
inline void swap(StrongWeakPtr<T> &x, StrongWeakPtr<T> &y) noexcept {
    x.swap(y);
}

template<class T>
inline bool operator==(const StrongWeakPtr<T> &x, const StrongWeakPtr<T> &y) noexcept {
    return x.get() == y.get();
}

template<class T>
inline bool operator!=(const StrongWeakPtr<T> &x, const StrongWeakPtr<T> &y) noexcept {
    return !(x == y);
}

template<class T1, class T2>
inline bool operator<(const StrongWeakPtr<T1> &x, const StrongWeakPtr<T2> &y) noexcept {
    typedef typename StrongWeakPtr<T1>::pointer P1;
    typedef typename StrongWeakPtr<T2>::pointer P2;
    typedef typename std::common_type<P1, P2>::type Common;

    using namespace std;
    return less<Common>()(x.get(), y.get());
}

template<class T1, class T2>
inline bool operator>(const StrongWeakPtr<T1> &x, const StrongWeakPtr<T2> &y) noexcept {
    return y < x;
}

template<class T1, class T2>
inline bool operator<=(const StrongWeakPtr<T1> &x, const StrongWeakPtr<T2> &y) noexcept {
    return !(y < x);
}

template<class T1, class T2>
inline bool operator>=(const StrongWeakPtr<T1> &x, const StrongWeakPtr<T2> &y) noexcept {
    return !(x < y);
}

template<class T>
inline bool operator==(const StrongWeakPtr<T> &x, std::nullptr_t) noexcept {
    return !x;
}

template<class T>
inline bool operator==(std::nullptr_t, const StrongWeakPtr<T> &x) noexcept {
    return !x;
}

template<class T>
inline bool operator!=(const StrongWeakPtr<T> &x, std::nullptr_t) noexcept {
    return static_cast<bool>(x);
}

template<class T>
inline bool operator!=(std::nullptr_t, const StrongWeakPtr<T> &x) noexcept {
    return static_cast<bool>(x);
}

template<class T, class... Args>
StrongWeakPtr<T>
make_sw(Args&&... args) {
    return StrongWeakPtr<T>(new T(std::forward<Args>(args)...));
}

}
}

