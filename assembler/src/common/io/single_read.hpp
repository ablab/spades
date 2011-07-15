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
#include <iostream>
#include <cassert>
#include "common/sequence/quality.hpp"
#include "common/sequence/sequence.hpp"
#include "common/sequence/nucl.hpp"
#include "common/sequence/sequence_tools.hpp"
#include "simple_tools.hpp"

class SingleRead {
public:
  static const int PHRED_OFFSET = 33;
  static const int BAD_QUALITY_THRESHOLD = 2;
  
  SingleRead() : valid_(false) {}
  
  SingleRead(const std::string &name, const std::string &seq, const std::string &qual) 
    : name_(name), seq_(seq), qual_(qual) { // for test only!
    valid_ = UpdateValid();
  }

  bool IsValid() const {
    return valid_;
  }

  Sequence sequence() const {
    assert(valid_);
    return Sequence(seq_);
  }

  Quality quality() const {
    assert(valid_);
    return Quality(qual_);
  }

  const std::string& name() const {
    return name_;
  }

  size_t size() const { 
    return seq_.size();
  }

  const std::string& GetSequenceString() const { 
    return seq_;
  }

  const std::string& GetQualityString() const {  
    return qual_;
  }

  std::string GetPhredQualityString(int offset = PHRED_OFFSET) const {
    std::string res = qual_;
    for (size_t i = 0; i < res.size(); ++i) {
      res[i] += offset;
    }
    return res;
  }

  char operator[](size_t i) const {
    assert(is_nucl(seq_[i]));
    return dignucl(seq_[i]);
  }

  SingleRead operator!() const {
    std::string new_name;
    if (name_ == "" || name_[0] != '!') {
      new_name = '!' + name_;
    } else {
      new_name = name_.substr(1, name_.length());
    }
    return SingleRead(new_name, ReverseComplement(seq_), Reverse(qual_));
  }

private:
  std::string name_;
  std::string seq_;
  std::string qual_;
  bool valid_;
  friend class ireadstream;
  
  void set_name(const char* s) {
    name_ = s;
  }
  
  void set_quality(const char* s, int offset = PHRED_OFFSET) {
    qual_ = s;
    for (size_t i = 0; i < qual_.size(); ++i) {
      qual_[i] -= offset;
    }
  }
  
  void set_sequence(const char* s) {
    seq_ = s;
    valid_ = UpdateValid();
  }
  
  const bool UpdateValid() const {
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
};

#endif /* SINGLEREAD_HPP_ */
