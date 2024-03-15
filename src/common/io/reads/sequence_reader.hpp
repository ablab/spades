//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/reads/ireader.hpp"
#include "basic/reads/single_read.hpp"

namespace io {

//todo merge with VectorReader
template<class ReadType>
class SequenceReadStream : public ReadStream<ReadType> {
public:
    explicit SequenceReadStream(const Sequence &sequence, const std::string &name = "")
            : sequence_(sequence),
              name_(name),
              opened_(true),
              eof_(false) {
    }

    virtual ~SequenceReadStream() {
    }

    virtual bool is_open() {
        return opened_;
    }

    virtual bool eof() {
        return eof_;
    }

    virtual void close() {
        opened_ = false;
    }

    void reset() {
        eof_ = false;
        opened_ = true;
    }

    SequenceReadStream &operator>>(ReadType &read);

private:
    Sequence sequence_;
    std::string name_;
    bool opened_;
    bool eof_;
};

template<>
SequenceReadStream<SingleRead> &SequenceReadStream<SingleRead>::operator>>(SingleRead &read) {
    if (!eof_) {
        read = SingleRead(name_, sequence_.str());
        eof_ = true;
    }
    return *this;
}

template<>
SequenceReadStream<SingleReadSeq> &SequenceReadStream<SingleReadSeq>::operator>>(SingleReadSeq &read) {
    if (!eof_) {
        read = SingleReadSeq(sequence_);
        eof_ = true;
    }
    return *this;
}

}
