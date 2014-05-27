//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "delegating_reader_wrapper.hpp"

namespace io {

//todo refactor!!!
class IdSettingReaderWrapper: public DelegatingWrapper<SingleRead> {
	typedef DelegatingWrapper<SingleRead> base;
	size_t next_id_;
public:
	IdSettingReaderWrapper(base::ReadStreamPtrT reader, size_t start_id = 0) :
			base(reader), next_id_(start_id) {
	}

	/* virtual */
	IdSettingReaderWrapper& operator>>(SingleRead& read) {
		this->reader() >> read;
		read.ChangeName(ToString(next_id_++));
		return *this;
	}
};

class PrefixAddingReaderWrapper: public DelegatingWrapper<SingleRead> {
	typedef DelegatingWrapper<SingleRead> base;
	std::string prefix_;
public:
	PrefixAddingReaderWrapper(base::ReadStreamPtrT reader,
			const std::string& prefix) :
			base(reader), prefix_(prefix) {
	}

	/* virtual */
	PrefixAddingReaderWrapper& operator>>(SingleRead& read) {
		this->reader() >> read;
		read.ChangeName(prefix_ + read.name());
		return *this;
	}
};

}
