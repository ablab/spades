//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/library.hpp"

#include <functional>

namespace io {

static inline std::pair<bool, bool> GetRCFlags(LibraryOrientation orientation) {
    switch (orientation) {
        case LibraryOrientation::RR:
            return {true, true};
        case LibraryOrientation::FR:
            return {false, true};
        case LibraryOrientation::RF:
            return {true, false};
        default:
            return {false, false};
    }
}

template<typename ReadType>
ReadType GetRCRead(const ReadType &read, bool rc) {
    return rc ? !read : read;
}

template<class ReadType>
using OrientationF = std::function<ReadType (const ReadType&)>;

template<typename ReadType>
OrientationF<ReadType> GetOrientationChanger(LibraryOrientation orientation) {
    auto rc_flags = GetRCFlags(orientation);
    return [=](const ReadType &r) {
        return ReadType(GetRCRead(r.first(), rc_flags.first),
                        GetRCRead(r.second(), rc_flags.second),
                        r.orig_insert_size());
    };
}

}
