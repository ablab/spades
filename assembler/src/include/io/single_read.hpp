/**
 * @file    single_read.hpp
 * @author  Mariya Fomkina
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * SingleRead is a structure, where information from input files is stored.
 * It includes 3 strings: with id, sequence and quality of the input read.
 */

#ifndef COMMON_IO_SINGLEREAD_HPP_
#define COMMON_IO_SINGLEREAD_HPP_

#include <string>
#include "verify.hpp"
#include "sequence/quality.hpp"
#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"
#include "sequence/sequence_tools.hpp"
#include "simple_tools.hpp"

namespace io {

/* 
 * This enumerate contains offset type.
 * UnknownOffset is equal to "offset = 0".
 * PhredOffset is equal to "offset = 33".
 * SolexaOffset is equal to "offset = 64".
 */
enum OffsetType {
	UnknownOffset, PhredOffset, SolexaOffset
};

//todo extract code about offset from here
class SingleRead {

public:
	static string EmptyQuality(const string& seq) {
		return std::string(seq.size(), (char) 33);
	}

	static const int UNKNOWN_OFFSET = 0;
	static const int PHRED_OFFSET = 33;
	static const int SOLEXA_OFFSET = 64;
	static const int BAD_QUALITY_THRESHOLD = 2;

	/*
	 * Type of variables which will store file names for reading from
	 * Reader stream.
	 */
	typedef std::string FilenameType;

	/*
	 * Default constructor.
	 */
	SingleRead() :
			name_(""), seq_(""), qual_(""), valid_(false) {
	}

	/*
	 * Test constructor.
	 *
	 * @param name The name of the SingleRead (id in input file).
	 * @param seq The sequence of ATGC letters.
	 * @param qual The quality of the SingleRead sequence.
	 */
	SingleRead(const std::string& name, const std::string& seq,
			const std::string& qual) :
			name_(name), seq_(seq), qual_(qual) {
		Init();
	}

	SingleRead(const std::string& name, const std::string& seq) :
			name_(name), seq_(seq), qual_(EmptyQuality(seq_)) {
		Init();
	}

	/*
	 * Check whether SingleRead is valid.
	 *
	 * @return true if SingleRead is valid (there is no N in sequence
	 * and sequence size is equal to quality size), and false otherwise
	 */
	bool IsValid() const {
		return valid_;
	}

	/*
	 * Return Sequence object, got from sequence string.
	 *
	 * @return SingleRead sequence.
	 */
	Sequence sequence() const {
		VERIFY(valid_);
		return Sequence(seq_);
	}

	/*
	 * Return Quality object, got from quality string.
	 *
	 * @return SingleRead quality.
	 */
	Quality quality() const {
		VERIFY(valid_);
		return Quality(qual_);
	}

	/*
	 * Return name of single read.
	 *
	 * @return SingleRead name.
	 */
	const std::string& name() const {
		return name_;
	}

	/*
	 * Return size of SingleRead.
	 *
	 * @return The size of SingleRead sequence.
	 */
	size_t size() const {
		return seq_.size();
	}

	/*
	 * Return SingleRead sequence string (in readable form with ATGC).
	 *
	 * @return SingleRead sequence string.
	 */
	const std::string& GetSequenceString() const {
		return seq_;
	}

	/*
	 * Return SingleRead quality string (in readable form).
	 *
	 * @return SingleRead quality string.
	 */
	const std::string& GetQualityString() const {
		return qual_;
	}

	/*
	 * Return SingleRead quality string, where every quality value is
	 * increased by PhredOffset (need for normalization of quality values).
	 * Do not modify original quality values.
	 *
	 * @return Modified SingleRead quality string.
	 */
	std::string GetPhredQualityString() const {
		int offset = PHRED_OFFSET;
		std::string res = qual_;
		for (size_t i = 0; i < res.size(); ++i) {
			res[i] += offset;
		}
		return res;
	}

	/*
	 * Return ith nucleotide of SingleRead sequence in unreadable form
	 * (0, 1, 2 or 3).
	 *
	 * @param i Nucleotide index.
	 * @return Nucleotide on ith position of SingleRead sequence.
	 */
	char operator[](size_t i) const {
		VERIFY(is_nucl(seq_[i]));
		return dignucl(seq_[i]);
	}

	/*
	 * Return reversed complimentary SingleRead (SingleRead with new
	 * name, reversed complimentary sequence, and reversed quality).
	 *
	 * @return Reversed complimentary SingleRead.
	 */
	SingleRead operator!() const {
		std::string new_name;
		if (name_.length() >= 3 && name_.substr(name_.length() - 3) == "_RC") {
			new_name = name_.substr(0, name_.length() - 3);
		} else {
			new_name = name_ + "_RC";
		}
//		TODO make naming nicer
//		if (name_ == "" || name_[0] != '!') {
//			new_name = '!' + name_;
//		} else {
//			new_name = name_.substr(1, name_.length());
//		}
		return SingleRead(new_name, ReverseComplement(seq_), Reverse(qual_));
	}

	SingleRead SubstrStrict(size_t from, size_t to) const {
		size_t len = to - from;
		return SingleRead(name_, seq_.substr(from, len), qual_.substr(from, len));
//		TODO make naming nicer
//		std::string new_name = name_ + ".substr(" + ToString(from) + "," + ToString(to) + ")";
//		return SingleRead(new_name, seq_.substr(from, len), qual_.substr(from, len));
	}

	SingleRead Substr(size_t from, size_t to) const {
		size_t len = to - from;
		if (len == size()) {
			return *this;
		}
		if (len == 0) {
			return SingleRead();
		}
		return SubstrStrict(from, to);
	}

	/*
	 * Check whether two SingleReads are equal.
	 *
	 * @param singleread The SingleRead we want to compare ours with.
	 *
	 * @return true if these two single reads have similar sequences,
	 * and false otherwise.
	 */
	bool operator==(const SingleRead& singleread) const {
		return seq_ == singleread.seq_;
	}

	void ChangeName(const string& new_name) {
		name_ = new_name;
	}

//	void ClearQuality() {
//		qual_ = std::string(seq_.size(), (char) 0);
//		UpdateValid();
//	}

	static bool IsValid(const string& seq) {
		for (size_t i = 0; i < seq.size(); ++i) {
			if (!is_nucl(seq[i])) {
				return false;
			}
		}
		return true;
	}

private:
	/*
	 * @variable The name of SingleRead in input file.
	 */
	std::string name_;
	/*
	 * @variable The sequence of nucleotides.
	 */
	std::string seq_;
	/*
	 * @variable The quality of SingleRead.
	 */
	std::string qual_;
	/*
	 * @variable The flag of SingleRead correctness.
	 */
	bool valid_;

	/*
	 * todo move from here!!!
	 * Return quality offset value from offset type.
	 *
	 * @param offset_type One of possible enum values.
	 * @return Quality offset value.
	 */
	int GetOffset(OffsetType offset_type) const {
		switch (offset_type) {
		case UnknownOffset:
			return 0;
		case PhredOffset:
			return 33;
		case SolexaOffset:
			return 64;
		}
		return -1;
	}

//	/*
//	 * Set quality of SingleRead.
//	 *
//	 * @param new_quality New quality of SingleRead.
//	 * @param offset The offset of SingleRead quality
//	 * (PHRED_OFFSET by default).
//	 */
//	void SetQuality(const string& new_quality, OffsetType offset_type =
//			PhredOffset) {
//		int offset = GetOffset(offset_type);
//		qual_ = new_quality;
//		for (size_t i = 0; i < qual_.size(); ++i) { // oh, really does it work with char* copying?
//			qual_[i] -= offset;
//		}
//		cout << "Set quality" << endl;
//	}

	void Init() {
		VERIFY(seq_.size() == qual_.size());
		valid_ = SingleRead::IsValid(seq_);
	}
};

}

#endif /* COMMON_IO_SINGLEREAD_HPP_ */
