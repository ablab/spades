/*
 * read.hpp
 *
 *  Created on: 29.03.2011
 *      Author: vyahhi
 */

#ifndef READ_HPP_
#define READ_HPP_

#include "quality.hpp"
#include "sequence.hpp"
#include "nucl.hpp"
#include <string>
#include <iostream>
#include "simple_tools.hpp"
#include "strobe_reader.hpp"
using namespace std;

class Read {
public:
	static const int PHRED_OFFSET = 33;

	bool isValid() const {
		return valid;
	}

	Sequence getSequence() const {
		assert(valid);
		return Sequence(seq_);
	}

	Quality* createQuality() const {
		assert(valid);
		return new Quality(qual_);
	}

	const string& getSequenceString() const {
		return seq_;
	}
	const string& getQualityString() const{
		return qual_;
	}
	const string& getName() const {
		return name_;
	}
	size_t size() const {
		return seq_.size();
	}
	char operator[](size_t i) const {
		assert(is_nucl(seq_[i]));
		return dignucl(seq_[i]);
	}
	void trimNs() {
		size_t index = seq_.find('N');
		if (index != string::npos) {
			seq_.erase(seq_.begin() + index, seq_.end());
			qual_.erase(qual_.begin() + index, qual_.end());
		}
	}
	Read() : valid(false) {
		;
	}
	Read(const string &name, const string &seq, const string &qual) : name_(name), seq_(seq), qual_(qual) { // for test only!
		valid = updateValid();
	}
private:
	string name_;
	string seq_;
	string qual_;
	bool valid;
	friend class ireadstream;
	void setName(const char* s) {
		name_ = s;
	}
	void setQuality(const char* s) {
		qual_ = s;
		for (size_t i = 0; i < qual_.size(); ++i) {
			qual_[i] -= PHRED_OFFSET;
		}
	}
	void setSequence(const char* s) {
		seq_ = s;
		valid = updateValid();
	}
	const bool updateValid() const {
		if (seq_.size() == 0) {
			return false;
		}
		for (size_t i = 0; i < seq_.size(); ++i) {
			if (!is_nucl(seq_[i])) {
				return false;
			}
		}
		return true;
	}


public:
	Read operator!() const {
		string newName;
		if(name_ == "" || name_[0] != '!')
			newName = '!' + name_;
		else
			newName = name_.substr(1, name_.length());
		return Read(newName, ReverseComplement(seq_), Reverse(qual_));
	}
};

class PairedRead {
	Read first_;
	Read second_;
	size_t distance_;
public:

	PairedRead() {

	}

	PairedRead(const Read& first, const Read& second, size_t distance) : first_(first), second_(second), distance_(distance) {

	}

	const Read& first() const {
		return first_;
	}

	const Read& second() const {
		return second_;
	}

	size_t distance() const {
		return distance_;
	}

	bool IsValid() const {
		return first_.isValid() && second_.isValid();
	}

	const Read& operator[] (size_t index) const {
		if (index == 0) {
			return first_;
		} else if (index == 1) {
			return second_;
		}
		assert(false);
	}

	const PairedRead operator!() const{
		return PairedRead(!second_, !first_, distance_);
	}
};

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

template<class Stream>
class RCReaderWrapper {
public:
//	typedef typename Stream::ReadType ReadType;
private:
	Stream &inner_reader_;
	PairedRead rc_result_;
	bool was_rc_;
public:
	RCReaderWrapper(Stream &reader) :
		inner_reader_(reader), was_rc_(false) {
	}

	bool eof() const {
		return inner_reader_.eof() && !was_rc_;
	}

	RCReaderWrapper& operator>>(PairedRead& p_r) {
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
		return inner_reader_.eof();
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

#endif /* READ_HPP_ */
