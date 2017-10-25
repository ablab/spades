//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/library.hpp"

namespace io {

template<class ReadType>
using OrientationF = std::function<ReadType (const ReadType&)>;

template<typename ReadType>
OrientationF<ReadType> GetOrientationChanger(LibraryOrientation orientation) {
    switch (orientation) {
    case LibraryOrientation::RR:  {
        return [](const ReadType &r) {
            return ReadType::Create(!r.first(), !r.second(), r.orig_insert_size());
        };
    }
    case LibraryOrientation::FR:  {
        return [](const ReadType &r) {
            return ReadType::Create(r.first(), !r.second(), r.orig_insert_size());
        };
    }
    case LibraryOrientation::RF:  {
        return [](const ReadType &r) {
            return ReadType::Create(!r.first(), r.second(), r.orig_insert_size());
        };
    }
    case LibraryOrientation::FF:
    default: {
        return [](const ReadType &r) {
            return ReadType(r);
        };
    }
    }
}

}
