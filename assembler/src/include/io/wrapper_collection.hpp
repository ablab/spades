//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "delegating_reader_wrapper.hpp"

namespace io {

//todo refactor!!!
class IdSettingReaderWrapper: public io::DelegatingReaderWrapper<io::SingleRead> {
	typedef io::DelegatingReaderWrapper<io::SingleRead> base;
	size_t next_id_;
public:
	IdSettingReaderWrapper(io::IReader<io::SingleRead>& reader, size_t start_id = 0) :
			base(reader), next_id_(start_id) {

	}

	/* virtual */
	IdSettingReaderWrapper& operator>>(io::SingleRead& read) {
		this->reader() >> read;
		read.ChangeName(ToString(next_id_++));
		return *this;
	}
};

class PrefixAddingReaderWrapper: public io::DelegatingReaderWrapper<io::SingleRead> {
	typedef io::DelegatingReaderWrapper<io::SingleRead> base;
	string prefix_;
public:
	PrefixAddingReaderWrapper(io::IReader<io::SingleRead>& reader,
			const string& prefix) :
			base(reader), prefix_(prefix) {

	}

	/* virtual */
	PrefixAddingReaderWrapper& operator>>(io::SingleRead& read) {
		this->reader() >> read;
		read.ChangeName(prefix_ + read.name());
		return *this;
	}
};

}
