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
    explicit RCWrapper(typename base::ReadStreamPtrT reader) :
            base(reader), rc_read_(), was_rc_(true) {
    }

    bool eof() override {
        return was_rc_ && base::eof();
    }

    RCWrapper& operator>>(ReadType& read) override {
        if (was_rc_) {
            base::operator >>(read);
            rc_read_ = read;
        } else {
            read = !rc_read_;
        }
        was_rc_ = !was_rc_;
        return (*this);
    }

    void reset() override {
        was_rc_ = true;
        base::reset();
    }

private:
    ReadType rc_read_;
    bool was_rc_;
};

template<class ReadType>
std::shared_ptr<ReadStream<ReadType>> RCWrap(std::shared_ptr<ReadStream<ReadType>> reader_ptr) {
    return std::make_shared<RCWrapper<ReadType>>(reader_ptr);
}

template<class ReadType>
ReadStreamList<ReadType> RCWrap(ReadStreamList<ReadType>& readers) {
    ReadStreamList<ReadType> answer;
    for (size_t i = 0; i < readers.size(); ++i) {
        answer.push_back(RCWrap<ReadType>(readers.ptr_at(i)));
    }
    return answer;
}

template<typename ReadType>
class OrientationChangingWrapper: public DelegatingWrapper<ReadType> {
    typedef DelegatingWrapper<ReadType> base;
public:

    OrientationChangingWrapper(typename base::ReadStreamPtrT reader,
                               LibraryOrientation orientation) :
            base(reader), changer_(GetOrientationChanger<ReadType>(orientation)) {
    }

    OrientationChangingWrapper& operator>>(ReadType& read) override {
        base::operator >>(read);
        read = changer_(read);
        return (*this);
    }

private:
    OrientationF<ReadType> changer_;
};

template<typename ReadType>
class RCRemovingWrapper: public DelegatingWrapper<ReadType> {
    typedef DelegatingWrapper<ReadType> base;
public:

    explicit RCRemovingWrapper(typename base::ReadStreamPtrT reader) : base(reader) {
    }

    /*virtual*/
    RCRemovingWrapper& operator>>(ReadType& read) {
        base::operator>>(read);

        VERIFY(!this->eof());
        ReadType skip;
        base::operator>>(skip);

        return *this;
    }

};

template<class ReadType>
std::shared_ptr<ReadStream<ReadType>> UnRCWrap(std::shared_ptr<ReadStream<ReadType>> reader_ptr) {
    return std::make_shared<RCRemovingWrapper<ReadType>>(reader_ptr);
}

template<class ReadType>
ReadStreamList<ReadType> UnRCWrap(ReadStreamList<ReadType>& readers) {
    ReadStreamList<ReadType> answer;
    for (size_t i = 0; i < readers.size(); ++i) {
        answer.push_back(UnRCWrap<ReadType>(readers.ptr_at(i)));
    }
    return answer;
}

}
