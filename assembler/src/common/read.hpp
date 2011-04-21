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

#endif /* READ_HPP_ */
