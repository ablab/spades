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

    void FillBuffer(const SingleRead& tmp_read) {
        buffer_.clear();
        buffer_position_ = 0;
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

    explicit SplittingWrapper(base::ReadStreamT reader) :
            base(std::move(reader)), buffer_position_(0) {
    }

    SplittingWrapper& operator>>(SingleRead& read) {
        Skip();
        read = buffer_[buffer_position_];
        buffer_position_++;
        return *this;
    }

    //todo fix needed!!! seems that eof can't be called multiple times in a row!!!
    bool eof() {
        return !Skip();
    }
};

inline ReadStream<SingleRead> SplittingWrap(ReadStream<SingleRead> reader_ptr) {
    return SplittingWrapper(std::move(reader_ptr));
}

inline ReadStreamList<SingleRead> SplittingWrap(ReadStreamList<SingleRead> readers) {
    ReadStreamList<SingleRead> answer;
    for (auto &reader : readers) {
        answer.push_back(SplittingWrap(std::move(reader)));
    }
    return answer;
}
}
