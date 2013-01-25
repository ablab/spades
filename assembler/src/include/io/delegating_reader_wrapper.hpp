//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "io/ireader.hpp"

namespace io {

template<typename ReadType>
class DelegatingReaderWrapper: public IReader<ReadType> {
public:

	explicit DelegatingReaderWrapper(IReader<ReadType>& reader) {
		Init(reader);
	}


	/* virtual */ ~DelegatingReaderWrapper() {
	}

	/* virtual */ bool is_open() {
		return reader_->is_open();
	}

	/* virtual */ bool eof() {
		return reader_->eof();
	}

	/* virtual */ DelegatingReaderWrapper& operator>>(ReadType& read) {
		(*reader_) >> read;
		return *this;
	}

	/* virtual */
	void close() {
		reader_->close();
	}

	/*
	 * Close the stream and open it again.
	 */
	/* virtual */
	void reset() {
		reader_->reset();
	}

    ReadStat get_stat() const {
            return reader_->get_stat();
    }

protected:
	DelegatingReaderWrapper() {

	}

	IReader<ReadType>& reader() {
		return *reader_;
	}

	void Init(IReader<ReadType>& reader) {
		reader_ = &reader;
	}

private:
	IReader<ReadType>* reader_;

	explicit DelegatingReaderWrapper(
			const DelegatingReaderWrapper& reader);
	/*
	 * Hidden assign operator.
	 */
	void operator=(const DelegatingReaderWrapper& reader);
};

}
