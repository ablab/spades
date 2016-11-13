//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <functional>

namespace func {

//to use with std::function-s
template<class T>
void Compose(T t, std::function<void(T)> f1,
        std::function<void(T)> f2) {
    if (f1)
        f1(t);
    if (f2)
        f2(t);
}

template<class T>
std::function<void(T)> Composition(std::function<void(T)> f1,
                                     std::function<void(T)> f2) {
    return std::bind(func::Compose<T>, std::placeholders::_1, f1, f2);
}

template<class T>
class Predicate {
public:
    typedef T checked_type;

    virtual bool Check(T t) const = 0;

    bool operator()(T t) const { return Check(t); }

    virtual ~Predicate() {}
};


}
