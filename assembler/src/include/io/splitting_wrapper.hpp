//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef N_REMOVING_WRAPPER_HPP_
#define N_REMOVING_WRAPPER_HPP_

#include "io/ireader.hpp"

namespace io {

class SplittingWrapper: public DelegatingReaderWrapper<SingleRead> {
public:
	typedef SingleRead ReadType;
private:
	std::vector<ReadType> buffer_;
	size_t buffer_position_;

	void FillBuffer(ReadType& tmp_read) {
		buffer_.clear();
		for(size_t i = 0; i < tmp_read.size(); i++) {
			size_t j = i;
			while(j < tmp_read.size() && is_nucl(tmp_read.GetSequenceString()[j])) {
				j++;
			}
			if(j > i) {
				buffer_.push_back(tmp_read.Substr(i, j));
				i = j - 1;
			}
		}
		buffer_position_ = 0;
	}

	bool Skip() {
		while(!this->reader().eof() && buffer_position_ == buffer_.size()) {
			ReadType tmp_read;
			this->reader() >> tmp_read;
			FillBuffer(tmp_read);
		}
		return buffer_position_ != buffer_.size();
	}

public:

	explicit SplittingWrapper(IReader<ReadType>& reader) :
			DelegatingReaderWrapper<ReadType>(reader), buffer_position_(0) {
	}

	/* virtual */
	SplittingWrapper& operator>>(ReadType& read) {
		Skip();
		read = buffer_[buffer_position_];
		buffer_position_++;
		return *this;
	}

	/* virtual */ bool eof() {
		return !Skip();
	}
};

}
#endif /* N_REMOVING_WRAPPER_HPP_ */
