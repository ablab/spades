/*
 * read.hpp
 *
 *  Created on: 29.03.2011
 *      Author: vyahhi
 */

#ifndef READ_HPP_
#define READ_HPP_

#include <string>
#include <iostream>
#include <cassert>
#include "common/sequence/quality.hpp"
#include "common/sequence/sequence.hpp"
#include "common/sequence/nucl.hpp"
#include "common/sequence/sequence_tools.hpp"
#include "simple_tools.hpp"

class Read {
public:
  static const int PHRED_OFFSET = 33;
  static const int BAD_QUALITY_THRESHOLD = 2;
  
  bool isValid() const {
    return valid_;
  }
  
  Sequence getSequence() const {
    assert(valid_);
    return Sequence(seq_);
  }
  Sequence getSubSequence(size_t start, size_t length) const {
    assert(length > 0 && start >= 0 && start + length <= seq_.size());
    return Sequence(seq_.substr(start, length));
  }
  
  Quality getQuality() const {
    assert(valid_);
    return Quality(qual_);
  }
  
  const std::string& getSequenceString() const {
    return seq_;
  }
  
  const std::string& getQualityString() const {
    return qual_;
  }
  
  std::string getPhredQualityString(int offset = PHRED_OFFSET) const {
    std::string res = qual_;
    for (size_t i = 0; i < res.size(); ++i) {
      res[i] += offset;
    }
    return res;
  }
  
  const std::string& getName() const {
    return name_;
  }
  
  size_t size() const {
    return seq_.size();
  }
  
  char operator[](size_t i) const {
    assert(is_nucl(seq_[i]));
    return dignucl(seq_[i]);
  }
  
  /**
   * It's actually not trim Ns, but trim everything before first 'N'
   * P.S. wtf? (Kolya)
   */
  void trimNs() __attribute__ ((deprecated)) {
    size_t index = seq_.find('N');
    if (index != std::string::npos) {
      seq_.erase(seq_.begin() + index, seq_.end());
      qual_.erase(qual_.begin() + index, qual_.end());
    }
    valid_ = updateValid();
  }
  
  /**
   * trim bad quality nucleotides from start and end of the read
   * @return size of the read left
   */
  size_t trimBadQuality() __attribute__ ((deprecated)) {
    size_t start = 0;
    for (; start < seq_.size(); ++start) {
      if (qual_[start] > BAD_QUALITY_THRESHOLD)
	break;
    }
    if (start != seq_.size()) {
      seq_.erase(seq_.begin(), seq_.begin() + start);
      qual_.erase(qual_.begin(), qual_.begin() + start);
      size_t end = seq_.size();
      for (; end > 0; --end) {
	if (qual_[end] > BAD_QUALITY_THRESHOLD)
	  break;
      }
      seq_.erase(seq_.begin() + end + 1, seq_.end());			
      qual_.erase(qual_.begin() + end + 1, qual_.end());
      valid_ = updateValid();
      return seq_.size();
    } else {
      seq_ = "";
      qual_ = "";
      valid_ = updateValid();
      return 0;
    }
  }
  /**
   * @param k k as in k-mer
   * @param start start point
   * @return the first starting point of a valid k-mer >=start; return -1 if no such place exists
   */
  int firstValidKmer(size_t start, size_t k) const __attribute__ ((deprecated)) {
    size_t curHypothesis = start;
    size_t i = start;
    for (; i < seq_.size(); ++i) {
      if (i >= k + curHypothesis)
	return curHypothesis;
      if (!is_nucl(seq_[i])) {
	curHypothesis = i + 1;
      }
    }
    if (i >= k + curHypothesis) {
      return curHypothesis;
    }
    return -1;
  }
  
  Read() :
    valid_(false) {
    ;
  }
  
  Read(const std::string &name, const std::string &seq, const std::string &qual) :
    name_(name), seq_(seq), qual_(qual) { // for test only!
    valid_ = updateValid();
  }
private:
  std::string name_;
  std::string seq_;
  std::string qual_;
  bool valid_;
  friend class ireadstream;
  void setName(const char* s) {
    name_ = s;
  }
  void setQuality(const char* s, int offset = PHRED_OFFSET) {
    qual_ = s;
    for (size_t i = 0; i < qual_.size(); ++i) {
      qual_[i] -= offset;
    }
  }
  void setSequence(const char* s) {
    seq_ = s;
    valid_ = updateValid();
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
    std::string newName;
    if (name_ == "" || name_[0] != '!') {
      newName = '!' + name_;
    } else {
      newName = name_.substr(1, name_.length());
    }
    return Read(newName, ReverseComplement(seq_), Reverse(qual_));
  }
};

// todo: put this to *.cpp
//ostream& operator<<(ostream& os, const Read& read) {
//	return os << read.getSequenceString();
//}

#endif /* READ_HPP_ */
