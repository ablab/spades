//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/library.hpp"

namespace io {

static inline std::pair<bool, bool> GetRCFlags(LibraryOrientation orientation) {
    switch (orientation) {
        case LibraryOrientation::RR: {
            return make_pair(true, true);
        }
        case LibraryOrientation::FR: {
            return make_pair(false, true);
        }
        case LibraryOrientation::RF: {
            return make_pair(true, false);
        }
        case LibraryOrientation::FF:
        default: {
            return make_pair(false, false);
        }
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
