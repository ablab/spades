#pragma once

#include "io/delegating_reader_wrapper.hpp"

class SequenceModifier {
public:
	virtual ~SequenceModifier() {}
	virtual Sequence Modify(const Sequence& s) const = 0;
};

/**
 * Attention!!! this class clears quality!!!
 */
template<class ReadType>
class ModifyingWrapper;

template<>
class ModifyingWrapper<io::SingleRead>: public io::DelegatingReaderWrapper<io::SingleRead> {
	shared_ptr<const SequenceModifier> modifier_;

public:
	ModifyingWrapper(io::IReader<io::SingleRead>& reader, shared_ptr<const SequenceModifier> modifier) :
			io::DelegatingReaderWrapper<io::SingleRead>(reader), modifier_(modifier) {}

	ModifyingWrapper& operator>>(io::SingleRead& read) {
		this->reader() >> read;
		read = io::SingleRead(read.name(), modifier_->Modify(read.sequence()).str());
		return *this;
	}
};

template<>
class ModifyingWrapper<io::PairedRead>: public io::DelegatingReaderWrapper<io::PairedRead> {
	shared_ptr<const SequenceModifier> modifier_;

public:
	ModifyingWrapper(io::IReader<io::PairedRead>& reader, shared_ptr<const SequenceModifier> modifier) :
			io::DelegatingReaderWrapper<io::PairedRead>(reader), modifier_(modifier) {}

	ModifyingWrapper& operator>>(io::PairedRead& read) {
		this->reader() >> read;
		read = io::PairedRead(io::SingleRead(read.first().name(), modifier_->Modify(read.first().sequence()).str())
				, io::SingleRead(read.first().name(), modifier_->Modify(read.second().sequence()).str()), read.insert_size());
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
