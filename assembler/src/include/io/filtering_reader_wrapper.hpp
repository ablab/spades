#ifndef COMMON_IO_FILTERINGREADERWRAPPER_HPP_
#define COMMON_IO_FILTERINGREADERWRAPPER_HPP_

#include "io/ireader.hpp"

namespace io {

template<typename ReadType>
class FilteringReaderWrapper: public IReader<ReadType> {
public:

	explicit FilteringReaderWrapper(IReader<ReadType>& reader) :
			reader_(reader), eof_(false) {
		StepForward();
	}

	/* virtual */
	~FilteringReaderWrapper() {
		close();
	}

	/* virtual */
	bool is_open() {
		return reader_.is_open();
	}

	/* virtual */
	bool eof() {
		return eof_;
	}

	/* virtual */
	FilteringReaderWrapper& operator>>(ReadType& read) {
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

private:
	IReader<ReadType>& reader_;

	bool eof_;

	ReadType next_read_;

	void StepForward() {
		while (!reader_.eof()) {
			reader_ >> next_read_;
			if (next_read_.IsValid()) {
				return;
			}
		}
		eof_ = true;
	}

	/*
	 * Hidden copy constructor.
	 */
	explicit FilteringReaderWrapper(
			const FilteringReaderWrapper<ReadType>& reader);
	/*
	 * Hidden assign operator.
	 */
	void operator=(const FilteringReaderWrapper<ReadType>& reader);
};

}

#endif /* COMMON_IO_FILTERINGREADERWRAPPER_HPP_ */
