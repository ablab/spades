//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/noncopyable.hpp>

#include "read_stream_vector.hpp"
#include "delegating_reader_wrapper.hpp"
#include "orientation.hpp"

namespace io {

/**
 * RCWrapper is the class-wrapper that gets reads and reverse
 * complimentary reads from given reader (one by one).
 */
template<typename ReadType>
class RCWrapper: public DelegatingWrapper<ReadType> {
    typedef DelegatingWrapper<ReadType> base;
public:
    explicit RCWrapper(typename base::ReadStreamT reader)
            : base(std::move(reader)), rc_read_(), was_rc_(true) {}

    bool eof() {
        return was_rc_ && base::eof();
    }

    RCWrapper& operator>>(ReadType& read) {
        if (was_rc_) {
            base::operator >>(read);
            rc_read_ = read;
        } else {
            read = !rc_read_;
        }
        was_rc_ = !was_rc_;
        return (*this);
    }

    void reset() {
        was_rc_ = true;
        base::reset();
    }

private:
    ReadType rc_read_;
    bool was_rc_;
};

template<class ReadType>
ReadStream<ReadType> RCWrap(ReadStream<ReadType> reader_ptr) {
    return RCWrapper<ReadType>(std::move(reader_ptr));
}

template<class ReadType>
ReadStreamList<ReadType> RCWrap(ReadStreamList<ReadType> readers) {
    ReadStreamList<ReadType> answer;
    for (auto &reader : readers)
        answer.push_back(RCWrap(std::move(reader)));

    return answer;
}

template<typename ReadType>
class OrientationChangingWrapper: public DelegatingWrapper<ReadType> {
    typedef DelegatingWrapper<ReadType> base;
public:

    OrientationChangingWrapper(typename base::ReadStreamT reader,
                               LibraryOrientation orientation)
            : base(std::move(reader)), changer_(orientation) {}

    OrientationChangingWrapper& operator>>(ReadType& read) {
        base::operator >>(read);
        changer_(read);
        return (*this);
    }

private:
    OrientationChanger changer_;
};

}
