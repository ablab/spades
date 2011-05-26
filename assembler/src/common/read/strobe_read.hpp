/*
 * strobe_read.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef STROBE_READ_HPP_
#define STROBE_READ_HPP_

#include "seq.hpp"
#include "read.hpp"
#include "paired_read.hpp"

using namespace std;

/*
 * Use single_read<size,T>::type for single read
 * Use mate_read<size,T>::type for mate read
 * Use strobe_read<size,cnt,T>::type for cnt-read
 *
 * where size -- number of nucleotides in each read
 * (we don't store info about gaps' size here)
 */
template <size_t size, size_t cnt = 1, typename T = int>
class strobe_read {
public:
	strobe_read() __attribute__ ((deprecated))  {
		// random sequences constructor
	}

	strobe_read(const string *ss) __attribute__ ((deprecated)) {
		for (size_t i = 0; i < cnt; ++i) {
			put(i, ss[i]);
		}
	}

	void put(int i, const string &s) { // by value (probably faster)
		data_[i] = Seq<size,T>(s);
	}

	Seq<size,T> operator[](size_t i) const {
		return data_[i];
	}

private:
	array<Seq<size,T>,cnt> data_;
};

// use mate_read<size, T>::type for mate reads
template <int size, typename T = char>
struct mate_read {
	typedef strobe_read<size,2,T> type; // because "There is no direct way to have templated typedefs in C++" :(
};

// use single_read<size, T>::type for single reads
template <int size, typename T = char>
struct single_read {
	typedef strobe_read<size,1,T> type; // because "There is no direct way to have templated typedefs in C++" :(
};


/////////////////////////////////



template<size_t cnt, typename T, typename TR>
class StrobeReader {
	vector<TR*> readers_;
public:
	typedef T ReadType;
	//	StrobeReader(const TR **readers) {
	//		for (size_t i = 0; i < cnt; ++i) {
	//			readers_[i] = readers[i];
	//		}
	//	}

	StrobeReader(const string filenames[]) {
		stringstream s;
		for (size_t i = 0; i < cnt; ++i) {
			readers_.push_back(new TR(filenames[i]));
		}
	}

	virtual ~StrobeReader() {
		close();
		for (size_t i = 0; i < cnt; ++i) {
			delete readers_[i];
		}
	}

	bool eof() const {
		for (size_t i = 0; i < cnt; ++i) {
			if (readers_[i]->eof()) {
				return true;
			}
		}
		return false;
	}

	StrobeReader& operator>>(vector<T>& v) {
		v.clear();
		T t;
		for (size_t i = 0; i < cnt; ++i) {
			(*readers_[i]) >> t;
			v.push_back(t);
		}
		return *this;
	}

	void reset() {
		for (size_t i = 0; i < cnt; ++i) {
			readers_[i]->reset();
		}
	}

	void close() {
		for (size_t i = 0; i < cnt; ++i) {
			readers_[i]->close();
		}
	}
};

template<typename T, typename TR>
struct MateReader {
	typedef StrobeReader<2, T, TR> type;
};

template<typename T, typename TR>
struct SingleReader {
	typedef StrobeReader<1, T, TR> type;
};

template<typename TR>
class CuttingReader {
	TR reader_;
	size_t cut_;
	size_t read_;
public:
	CuttingReader(TR reader, size_t cut = -1) : reader_(reader), cut_(cut), read_(0) {}
	virtual ~CuttingReader() {}

	bool eof() const {
		return read_ == cut_ || reader_.eof();
	}

	template<typename T>
	CuttingReader& operator>>(T& v) {
		reader_ >> v;
		++read_;
		return *this;
	}

	void reset() {
		read_ = 0;
		reader_.reset();
	}

	void close() {
		reader_.close();
	}
};

//////////////////

template <class TR>
class PairedReader {
public:
	typedef StrobeReader<2, Read, TR> InnerReader;
	typedef Read ReadType;
private:
	InnerReader& reader_;
	size_t distance_;
public:

	PairedReader(InnerReader& reader, size_t distance) : reader_(reader), distance_(distance) {

	}

	virtual ~PairedReader() {
		close();
	}

	bool eof() const {
		return reader_.eof();
	}

	PairedReader& operator>>(PairedRead& p_r) {
		vector<Read> v;
		reader_ >> v;
		p_r = PairedRead(v[0], !v[1], distance_);
		return *this;
	}

	void reset() {
		reader_.reset();
	}

	void close() {
		reader_.close();
	}
};

/////////////////////////////////


template<class Stream, typename ReadType>
class RCReaderWrapper {
public:
//	typedef typename Stream::ReadType ReadType;
private:
	Stream &inner_reader_;
	ReadType rc_result_;
	bool was_rc_;
public:
	RCReaderWrapper(Stream &reader):
		inner_reader_(reader), was_rc_(false) {
	}

	bool eof() const {
		return inner_reader_.eof() && !was_rc_;
	}

	RCReaderWrapper& operator>>(ReadType& p_r) {
		if (!was_rc_) {
			inner_reader_ >> p_r;
			rc_result_ = !p_r;
		} else {
			p_r = rc_result_;
		}
		was_rc_ = !was_rc_;
		return *this;
	}

	void reset() {
		was_rc_ = false;
		inner_reader_.reset();
	}

	void close() {
		inner_reader_.close();
	}
};

/////////////////////////////////


template<class Stream>
class SimpleReaderWrapper {
public:
//	typedef typename Stream::ReadType ReadType;
private:
	Stream &inner_reader_;
	PairedRead result_;
	size_t current_;
public:
	SimpleReaderWrapper(Stream &reader) :
		inner_reader_(reader), current_(0) {
	}

	bool eof() const {
		return inner_reader_.eof() && current_ == 0;
	}

	SimpleReaderWrapper& operator>>(Read& r) {
		if (current_ == 0) {
			inner_reader_ >> result_;
		}
		r = result_[current_];
		current_++;
		if (current_ == 2)
			current_ = 0;
		return *this;
	}

	void reset() {
		current_ = 0;
		inner_reader_.reset();
	}

	void close() {
		inner_reader_.close();
	}
};


#endif /* STROBE_READ_HPP_ */
