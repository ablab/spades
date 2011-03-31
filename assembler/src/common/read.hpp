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
using namespace std;

class Read {
public:
	const bool isValid() const {
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

	Sequence* createSequence() const {
		return new Sequence(seq_);
	}

	Quality* createQuality() const {
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
	Read() {
		;
	}
	Read(const string &name, const string &seq, const string &qual) : name_(name), seq_(seq), qual_(qual) { // for test only!
		;
	}
private:
	static const int PHRED_OFFSET = 33;
	string name_;
	string seq_;
	string qual_;
	friend class ireadstream;
	void setName(const char* s) {
		name_ = s;
	}
	void setQuality(const char* s, bool rtl = false) {
		qual_ = s;
		// TODO: HAHA! NEXT TWO BRANCHES DO THE SAME :)
		if (rtl) {
			for (size_t i = 0; i < qual_.size(); ++i) {
				qual_[qual_.size() - i - 1] -= PHRED_OFFSET;
			}
		} else {
			for (size_t i = 0; i < qual_.size(); ++i) {
				qual_[i] -= PHRED_OFFSET;
			}
		}
	}
	void setSequence(const char* s, bool rtl = false) {
		string tmp = s;
		seq_ = tmp;
		if (rtl) {
			int len = tmp.length();
			for(int i = 0; i < len; i++) {
				switch (tmp[i]) {
					case 'A':
						seq_[len - i - 1] = 'T';
						break;
					case 'C':
						seq_[len - i - 1] = 'G';
						break;
					case 'T':
						seq_[len - i - 1] = 'A';
						break;
					case 'G':
						seq_[len - i - 1] = 'C';
						break;
					case 'N':
						seq_[len - i - 1] = 'N';
						break;
					default:
						cerr << " strange letter in read. Exiting" << tmp; // TODO: wtf?
						assert(0);
				}
			}
		}
	}
};

#endif /* READ_HPP_ */
