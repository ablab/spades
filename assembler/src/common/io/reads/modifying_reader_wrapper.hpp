//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/verify.hpp"
#include "delegating_reader_wrapper.hpp"
#include "single_read.hpp"
#include "paired_readers.hpp"

#include <memory>

namespace io {

class SequenceModifier {
public:
    virtual ~SequenceModifier() {}

    SingleRead Modify(const SingleRead& read) {
        return SingleRead(read.name(), Modify(read.sequence()).str());
    }

    SingleReadSeq Modify(const SingleReadSeq& read) {
        return SingleReadSeq(Modify(read.sequence()));
    }

    virtual Sequence Modify(const Sequence& s) = 0;
};

class TrivialModifier : public SequenceModifier {
public:

    virtual Sequence Modify(const Sequence& s) {
        return s;
    }
};

/**
 * Attention!!! this class clears quality!!!
 */
template<class ReadType>
class ModifyingWrapper;

template<>
class ModifyingWrapper<SingleRead>: public DelegatingWrapper<SingleRead> {
  typedef DelegatingWrapper<SingleRead> base;
  std::shared_ptr<SequenceModifier> modifier_;

public:
    ModifyingWrapper(base::ReadStreamPtrT reader, std::shared_ptr<SequenceModifier> modifier) :
            base(reader), modifier_(modifier) {}

    ModifyingWrapper& operator>>(SingleRead& read) {
        this->reader() >> read;
        read = modifier_->Modify(read);
        return *this;
    }
};

template<>
class ModifyingWrapper<PairedRead>: public DelegatingWrapper<PairedRead> {
  typedef DelegatingWrapper<PairedRead> base;
  std::shared_ptr<SequenceModifier> modifier_;

public:
    ModifyingWrapper(base::ReadStreamPtrT reader, std::shared_ptr<SequenceModifier> modifier) :
            base(reader), modifier_(modifier) {}

    ModifyingWrapper& operator>>(PairedRead& read) {
        this->reader() >> read;
        read = PairedRead(modifier_->Modify(read.first()),
                          modifier_->Modify(read.second()),
                          read.orig_insert_size());
        return *this;
    }
};

template<>
class ModifyingWrapper<SingleReadSeq>: public DelegatingWrapper<SingleReadSeq> {
  typedef DelegatingWrapper<SingleReadSeq> base;
  std::shared_ptr<SequenceModifier> modifier_;

public:
  ModifyingWrapper(base::ReadStreamPtrT reader, std::shared_ptr<SequenceModifier> modifier) :
      base(reader), modifier_(modifier) {}

    ModifyingWrapper& operator>>(SingleReadSeq& read) {
        this->reader() >> read;
        read = SingleReadSeq(modifier_->Modify(read.sequence()), read.GetLeftOffset(), read.GetRightOffset());
        return *this;
    }
};

template<>
class ModifyingWrapper<PairedReadSeq>: public DelegatingWrapper<PairedReadSeq> {
  typedef DelegatingWrapper<PairedReadSeq> base;
  std::shared_ptr<SequenceModifier> modifier_;

public:
  ModifyingWrapper(base::ReadStreamPtrT reader, std::shared_ptr<SequenceModifier> modifier) :
            base(reader), modifier_(modifier) {}

    ModifyingWrapper& operator>>(PairedReadSeq& read) {
        this->reader() >> read;
        read = PairedReadSeq(modifier_->Modify(read.first()), modifier_->Modify(read.second()), read.orig_insert_size());
        return *this;
    }
};

}
