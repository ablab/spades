//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "verify.hpp"
#include "io/delegating_reader_wrapper.hpp"

#include <memory>

namespace io {

class SequenceModifier {
public:
	virtual ~SequenceModifier() {}

	virtual std::string Modify(const std::string& s) const {
	    return Modify(Sequence(s)).str();
    }

	virtual Sequence Modify(const Sequence& s) const = 0;
};

class TrivialModifier : public SequenceModifier {
public:

	virtual std::string Modify(const std::string& s) const {
	    return s;
    }

	virtual Sequence Modify(const Sequence& /*s*/) const {
	    VERIFY(false);
	    return Sequence();
	}
};

/**
 * Attention!!! this class clears quality!!!
 */
template<class ReadType>
class ModifyingWrapper;

template<>
class ModifyingWrapper<io::SingleRead>: public io::DelegatingReaderWrapper<io::SingleRead> {
  std::shared_ptr<const SequenceModifier> modifier_;

public:
	ModifyingWrapper(io::IReader<io::SingleRead>& reader, std::shared_ptr<const SequenceModifier> modifier) :
			io::DelegatingReaderWrapper<io::SingleRead>(reader), modifier_(modifier) {}

	ModifyingWrapper& operator>>(io::SingleRead& read) {
		this->reader() >> read;
		read = io::SingleRead(read.name(), modifier_->Modify(read.GetSequenceString()));
		return *this;
	}
};

template<>
class ModifyingWrapper<io::PairedRead>: public io::DelegatingReaderWrapper<io::PairedRead> {
  std::shared_ptr<const SequenceModifier> modifier_;

public:
	ModifyingWrapper(io::IReader<io::PairedRead>& reader, std::shared_ptr<const SequenceModifier> modifier) :
			io::DelegatingReaderWrapper<io::PairedRead>(reader), modifier_(modifier) {}

	ModifyingWrapper& operator>>(io::PairedRead& read) {
		this->reader() >> read;
		read = io::PairedRead(io::SingleRead(read.first().name(), modifier_->Modify(read.first().sequence()).str())
				, io::SingleRead(read.first().name(), modifier_->Modify(read.second().sequence()).str()), read.insert_size());
		return *this;
	}
};

template<>
class ModifyingWrapper<io::SingleReadSeq>: public io::DelegatingReaderWrapper<io::SingleReadSeq> {
  std::shared_ptr<const SequenceModifier> modifier_;

public:
  ModifyingWrapper(io::IReader<io::SingleReadSeq>& reader, std::shared_ptr<const SequenceModifier> modifier) :
      io::DelegatingReaderWrapper<io::SingleReadSeq>(reader), modifier_(modifier) {}

    ModifyingWrapper& operator>>(io::SingleReadSeq& read) {
        this->reader() >> read;
        read = io::SingleReadSeq(modifier_->Modify(read.sequence()));
        return *this;
    }
};

template<>
class ModifyingWrapper<io::PairedReadSeq>: public io::DelegatingReaderWrapper<io::PairedReadSeq> {
  std::shared_ptr<const SequenceModifier> modifier_;

public:
  ModifyingWrapper(io::IReader<io::PairedReadSeq>& reader, std::shared_ptr<const SequenceModifier> modifier) :
            io::DelegatingReaderWrapper<io::PairedReadSeq>(reader), modifier_(modifier) {}

    ModifyingWrapper& operator>>(io::PairedReadSeq& read) {
        this->reader() >> read;
        read = io::PairedReadSeq(io::SingleReadSeq(modifier_->Modify(read.first().sequence()))
            , io::SingleReadSeq(modifier_->Modify(read.second().sequence())), read.insert_size());
        return *this;
    }
};

//
//class ModifyingWrapper: public io::DelegatingReaderWrapper<io::SingleRead> {
//	SequenceModifier modifier_;
//protected:
//
//	ModifyingWrapper(io::IReader<io::SingleRead>& reader, const SequenceModifier& modifier) :
//			io::DelegatingReaderWrapper<io::SingleRead>(reader), modifier_(modifier) {}
//
//public:
//	ModifyingWrapper& operator>>(io::SingleRead& read) {
//		this->reader() >> read;
//		read = io::SingleRead(read.name(), modifier_.Modify(read.sequence()).str());
//		return *this;
//	}
//};
}
