//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "single_read.hpp"
#include "read_stream_vector.hpp"
#include "delegating_reader_wrapper.hpp"

namespace io {

class SplittingWrapper: public DelegatingWrapper<SingleRead> {
    typedef DelegatingWrapper<SingleRead> base;
private:
    std::vector<SingleRead> buffer_;
    size_t buffer_position_;

    void FillBuffer(SingleRead& tmp_read) {
        buffer_.clear();
        for (size_t i = 0; i < tmp_read.size(); ++i) {
            size_t j = i;
            while (j < tmp_read.size() && is_nucl(tmp_read.GetSequenceString()[j])) {
                j++;
            }
            if (j > i) {
                buffer_.push_back(tmp_read.Substr(i, j));
                i = j - 1;
            }
        }
        buffer_position_ = 0;
    }

    bool Skip() {
        while (!this->reader().eof() && buffer_position_ == buffer_.size()) {
            SingleRead tmp_read;
            this->reader() >> tmp_read;
            FillBuffer(tmp_read);
        }
        return buffer_position_ != buffer_.size();
    }

public:

    explicit SplittingWrapper(base::ReadStreamPtrT reader) :
            base(reader), buffer_position_(0) {
    }

    /* virtual */
    SplittingWrapper& operator>>(SingleRead& read) {
        Skip();
        read = buffer_[buffer_position_];
        buffer_position_++;
        return *this;
    }

    //todo fix needed!!! seems that eof can't be called multiple times in a row!!!
    /* virtual */ bool eof() {
        return !Skip();
    }
};

inline std::shared_ptr<ReadStream<SingleRead>> SplittingWrap(std::shared_ptr<ReadStream<SingleRead>> reader_ptr) {
    return std::make_shared<SplittingWrapper>(reader_ptr);
}

inline ReadStreamList<SingleRead> SplittingWrap(ReadStreamList<SingleRead>& readers) {
    ReadStreamList<SingleRead> answer;
    for (size_t i = 0; i < readers.size(); ++i) {
        answer.push_back(SplittingWrap(readers.ptr_at(i)));
    }
    return answer;
}
}
