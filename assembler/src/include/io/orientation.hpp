/*
 * orientation.hpp
 *
 *  Created on: May 16, 2013
 *      Author: andrey
 */

#ifndef ORIENTATION_HPP_
#define ORIENTATION_HPP_

#include "library.hpp"

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
OrientationChanger<ReadType> * GetOrientationChanger(LibraryOrientation orientation) {
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
    return result;
}

}


#endif /* ORIENTATION_HPP_ */
