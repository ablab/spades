//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <functional>

namespace func {

template<class T>
std::function<void(T)> CombineCallbacks(const std::function<void(T)>& f1,
                                        const std::function<void(T)>& f2) {
    return [=] (T t) {
        if (f1)
            f1(t);
        if (f2)
            f2(t);
    };
}

}
