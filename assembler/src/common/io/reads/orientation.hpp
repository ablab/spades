//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/library.hpp"

namespace io {

template<typename ReadType>
class OrientationChanger {

public:

    virtual ReadType Perform(const ReadType& r) const = 0;

    virtual ~OrientationChanger() {
    }
};

template<typename ReadType>
class IdeticalChanger : public OrientationChanger<ReadType> {

public:

    virtual ReadType Perform(const ReadType& r) const {
        return r;
    }
};

template<typename ReadType>
class ReverseSecondChanger : public OrientationChanger<ReadType> {

public:

    virtual ReadType Perform(const ReadType& r) const {
        return ReadType(r.first(), !r.second(), r.insert_size());
    }
};

template<typename ReadType>
class ReverseFirstChanger : public OrientationChanger<ReadType> {

public:

    virtual ReadType Perform(const ReadType& r) const {
        return ReadType(!r.first(), r.second(), r.insert_size());
    }
};

template<typename ReadType>
class ReverseChanger : public OrientationChanger<ReadType> {

public:

    virtual ReadType Perform(const ReadType& r) const {
        return ReadType(!r.first(), !r.second(), r.insert_size());
    }
};

template<typename ReadType>
std::unique_ptr<OrientationChanger<ReadType>> GetOrientationChanger(LibraryOrientation orientation) {
    OrientationChanger<ReadType> * result;
    switch (orientation) {
    case LibraryOrientation::FF:  {
        result = new IdeticalChanger<ReadType>();
        break;
    }
    case LibraryOrientation::RR:  {
        result = new ReverseChanger<ReadType>();
        break;
    }
    case LibraryOrientation::FR:  {
        result = new ReverseSecondChanger<ReadType>();
        break;
    }
    case LibraryOrientation::RF:  {
        result = new ReverseFirstChanger<ReadType>();
        break;
    }
    default: {
        result = new IdeticalChanger<ReadType>();
        break;
    }
    }
    return std::unique_ptr<OrientationChanger<ReadType>>(result);
}

}
