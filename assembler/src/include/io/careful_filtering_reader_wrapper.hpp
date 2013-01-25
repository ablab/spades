//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef COMMON_IO_CAREFULFILTERINGREADERWRAPPER_HPP_
#define COMMON_IO_CAREFULFILTERINGREADERWRAPPER_HPP_

#include "io/ireader.hpp"

namespace io {

const size_t none = -1;

pair<size_t, size_t> longestValidCoords(const SingleRead& r) {
	size_t best_len = 0;
	size_t best_pos = none;
	size_t pos = none;
	std::string seq = r.GetSequenceString();
	for (size_t i = 0; i <= seq.size(); ++i) {
		if (i < seq.size() && is_nucl(seq[i])) {
			if (pos == none) {
				pos = i;
			}
		} else {
			if (pos != none) {
				size_t len = i - pos;
				if (len > best_len) {
					best_len = len;
					best_pos = pos;
				}
			}
			pos = none;
		}
	}
	if (best_len == 0) {
		return make_pair(0, 0);
	}
	return make_pair(best_pos, best_pos + best_len);
}

SingleRead longestValid(const SingleRead& r) {
	pair<size_t, size_t> p = longestValidCoords(r);
	return r.Substr(p.first, p.second);
}

PairedRead longestValid(const PairedRead& r) {
	pair<size_t, size_t> c1 = longestValidCoords(r.first());
	pair<size_t, size_t> c2 = longestValidCoords(r.second());
	size_t len1 = c1.second - c1.first;
	size_t len2 = c2.second - c2.first;
	if (len1 == 0 || len2 == 0) {
		return PairedRead();
	}
	if (len1 == r.first().size() && len2 == r.second().size()) {
		return r;
	}
	size_t is = r.insert_size() - c1.first - r.second().size() + c2.second;

	return PairedRead(r.first().Substr(c1.first, c1.second), r.second().Substr(c2.first, c2.second), is);
}

template<typename ReadType>
class CarefulFilteringReaderWrapper : public IReader<ReadType> {
};

template<>
class CarefulFilteringReaderWrapper<SingleRead> : public IReader<SingleRead> {
public:
  /*
   * Default constructor.
   *
   * @param reader Reference to any other reader (child of IReader).
   */
	explicit CarefulFilteringReaderWrapper(IReader<SingleRead>& reader) :
			reader_(reader), eof_(false) {
		StepForward();
	}

  /* 
   * Default destructor.
   */
	/* virtual */ ~CarefulFilteringReaderWrapper() {
		close();
	}

  /* 
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
	/* virtual */ bool is_open() {
		return reader_.is_open();
	}

  /* 
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
	/* virtual */ bool eof() {
		return eof_;
	}

  /*
   * Read SingleRead from stream.
   *
   * @param read The SingleRead that will store read * data.
   *
   * @return Reference to this stream.
   */
	/* virtual */ CarefulFilteringReaderWrapper& operator>>(SingleRead& read) {
		read = next_read_;
		StepForward();
		return *this;
	}

	/*
	 * Close the stream.
	 */
	/* virtual */
	void close() {
		reader_.close();
	}

	/*
	 * Close the stream and open it again.
	 */
	/* virtual */
	void reset() {
		reader_.reset();
		eof_ = false;
		StepForward();
	}

	ReadStat get_stat() const {
        return reader_.get_stat();
    }

private:
  /*
   * @variable Internal stream readers.
   */
	IReader<SingleRead>& reader_;
  /*
   * @variable Flag that shows whether the end of stream reached.
   */
	bool eof_;
  /*
   * @variable Next read to be outputted by stream.
   */
	SingleRead next_read_;

  /*
   * Read next valid read in the stream.
   */
	void StepForward() {
		while (!reader_.eof()) {
			reader_ >> next_read_;
			next_read_ = longestValid(next_read_);
			if (next_read_.IsValid()) {
				return;
			}
			
		}
		eof_ = true;
	}

	/*
	 * Hidden copy constructor.
	 */
	explicit CarefulFilteringReaderWrapper(
			const CarefulFilteringReaderWrapper<SingleRead>& reader);
	/*
	 * Hidden assign operator.
	 */
	void operator=(const CarefulFilteringReaderWrapper<SingleRead>& reader);
};

template<>
class CarefulFilteringReaderWrapper<PairedRead>: public IReader<PairedRead> {
public:
  /*
   * Default constructor.
   *
   * @param reader Reference to any other reader (child of IReader).
   */
	explicit CarefulFilteringReaderWrapper(IReader<PairedRead>& reader) :
			reader_(reader), eof_(false) {
		StepForward();
	}

  /* 
   * Default destructor.
   */
	/* virtual */ ~CarefulFilteringReaderWrapper() {
		close();
	}

  /* 
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
	/* virtual */ bool is_open() {
		return reader_.is_open();
	}

  /* 
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
	/* virtual */ bool eof() {
		return eof_;
	}

  /*
   * Read PairedRead from stream.
   *
   * @param read The PairedRead that will store read * data.
   *
   * @return Reference to this stream.
   */
	/* virtual */ CarefulFilteringReaderWrapper& operator>>(PairedRead& read) {
		read = next_read_;
		StepForward();
		return *this;
	}

	/*
	 * Close the stream.
	 */
	/* virtual */
	void close() {
		reader_.close();
	}

	/*
	 * Close the stream and open it again.
	 */
	/* virtual */
	void reset() {
		reader_.reset();
		eof_ = false;
		StepForward();
	}

	ReadStat get_stat() const {
        return reader_.get_stat();
    }

private:
  /*
   * @variable Internal stream readers.
   */
	IReader<PairedRead>& reader_;
  /*
   * @variable Flag that shows whether the end of stream reached.
   */
	bool eof_;
  /*
   * @variable Next read to be outputted by stream.
   */
	PairedRead next_read_;

  /*
   * Read next valid read in the stream.
   */
	void StepForward() {
		while (!reader_.eof()) {
			reader_ >> next_read_;
			next_read_ = longestValid(next_read_);
			if (next_read_.IsValid()) {
				return;
			}
			
		}
		eof_ = true;
	}

	/*
	 * Hidden copy constructor.
	 */
	explicit CarefulFilteringReaderWrapper(
			const CarefulFilteringReaderWrapper<PairedRead>& reader);
	/*
	 * Hidden assign operator.
	 */
	void operator=(const CarefulFilteringReaderWrapper<PairedRead>& reader);
};

} // namespace io

#endif /* COMMON_IO_CAREFULFILTERINGREADERWRAPPER_HPP_ */
