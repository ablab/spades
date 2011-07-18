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

#ifndef SINGLEREAD_HPP_
#define SINGLEREAD_HPP_

#include <string>
#include <cassert>
#include "common/sequence/quality.hpp"
#include "common/sequence/sequence.hpp"
#include "common/sequence/nucl.hpp"
#include "common/sequence/sequence_tools.hpp"
#include "common/simple_tools.hpp"

class SingleRead {
 private:
  /*
   * @variable The name of single read in input file.
   */
  std::string name_;
  /*
   * @variable The sequence of nucleotides.
   */
  std::string seq_;
  /*
   * @variable The quality of single read.
   */
  std::string qual_;
  /*
   * @variable The flag of single read correctness.
   */
  bool valid_;
  friend class ireadstream;

  /*
   * Set name of single read.
   *
   * @param new_name New name.
   */
  void set_name(const char* new_name) {
    name_ = new_name;
  }

  /*
   * Set sequence of single read.
   *
   * @param new_sequence New sequence.
   */
  void set_sequence(const char* new_sequence) {
    seq_ = new_sequence;
    valid_ = UpdateValid();
  }

  /*
   * Set quality of single read.
   *
   * @param new_quality New quality of single read.
   * @param offset The offset of single read quality 
   * (PHRED_OFFSET by default).
   */
  void set_quality(const char* new_quality, int offset = PHRED_OFFSET) {
    qual_ = new_quality;
    for (size_t i = 0; i < qual_.size(); ++i) {
      qual_[i] -= offset;
    }
  }

  /*
   * Update valid_ flag.
   *
   * @see IsValid()
   */
  const bool UpdateValid() const {
    if (seq_.size() == 0) {
      return false;
    }
    if (seq_.size() != qual_.size()) {
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
  static const int PHRED_OFFSET = 33;
  static const int BAD_QUALITY_THRESHOLD = 2;

  /*
   * Default constructor.
   */
  SingleRead() : valid_(false) {}

  /*
   * Test constructor.
   *
   * @param name The name of the single read (id in input file).
   * @param seq The sequence of ATGC letters.
   * @param qual The quality of the single read sequence.
   */
  SingleRead(const std::string& name,
             const std::string& seq,
             const std::string& qual)
    : name_(name), seq_(seq), qual_(qual) {
    valid_ = UpdateValid();
  }

  /*
   * Check whether single read is valid.
   *
   * @return true if single read is valid (there is no N in sequence
   * and sequence size is equal to quality size), and false otherwise
   */
  bool IsValid() const {
    return valid_;
  }

  /*
   * Return Sequence object, got from sequence string.
   *
   * @return Single read sequence.
   */
  Sequence sequence() const {
    assert(valid_);
    return Sequence(seq_);
  }

  /*
   * Return Quality object, got from quality string.
   *
   * @return Single read quality.
   */
  Quality quality() const {
    assert(valid_);
    return Quality(qual_);
  }

  /*
   * Return name of single read.
   *
   * @return Single read name.
   */
  const std::string& name() const {
    return name_;
  }

  /*
   * Return size of single read.
   *
   * @return The size of single read sequence.
   */
  size_t size() const {
    return seq_.size();
  }

  /*
   * Return single read sequence string (in readable form with ATGC).
   *
   * @return Single read sequence string.
   */
  const std::string& GetSequenceString() const {
    return seq_;
  }

  /*
   * Return single read quality string (in readable form).
   *
   * @return Single read quality string.
   */
  const std::string& GetQualityString() const {
    return qual_;
  }

  /*
   * Return single read quality string, where every quality value is
   * increased by offset (need for normalization of quality values).
   * Do not modify original quality values.
   *
   * @param offset The offset of single read quality (PHRED_OFFSET by default).
   * @return Modified single read quality string.
   */
  std::string GetPhredQualityString(int offset = PHRED_OFFSET) const {
    std::string res = qual_;
    for (size_t i = 0; i < res.size(); ++i) {
      res[i] += offset;
    }
    return res;
  }

  /*
   * Return ith nucleotide of single read sequence in unreadable form
   * (0, 1, 2 or 3).
   *
   * @param i Nucleotide index.
   * @return Nucleotide on ith position of single read sequence.
   */
  char operator[](size_t i) const {
    assert(is_nucl(seq_[i]));
    return dignucl(seq_[i]);
  }

  /*
   * Return reversed complimentary single read (single read with new name,
   * reversed complimentary sequence, and reversed quality).
   *
   * @return Reversed complimentary single read.
   */
  SingleRead operator!() const {
    std::string new_name;
    if (name_ == "" || name_[0] != '!') {
      new_name = '!' + name_;
    } else {
      new_name = name_.substr(1, name_.length());
    }
    return SingleRead(new_name, ReverseComplement(seq_), Reverse(qual_));
  }

  /*
   * Check whether two single reads are equal.
   *
   * @param singleread The single read we want to compare ours with.
   * @return true if these two single reads have similar sequences,
   * and false otherwise.
   */
  bool operator==(const SingleRead& singleread) const {
    return seq_ == singleread.seq_;
  }
};

#endif /* SINGLEREAD_HPP_ */
