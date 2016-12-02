//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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
#include <fstream>
#include "utils/verify.hpp"
#include "sequence/quality.hpp"
#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"
#include "sequence/sequence_tools.hpp"
#include "utils/simple_tools.hpp"

//fixme deprecated!!! used in hammer!
class Read {
public:
    static const int PHRED_OFFSET = 33;

    bool isValid() const {
        return valid_;
    }

    Sequence getSequence() const {
        VERIFY(valid_);
        return Sequence(seq_);
    }

    Sequence getSubSequence(size_t start, size_t length) const __attribute__ ((deprecated)) {
        VERIFY(length > 0 && start + length <= seq_.size());
        return Sequence(seq_.substr(start, length));
    }

    Quality getQuality() const {
        VERIFY(valid_);
        return Quality(qual_);
    }

    const std::string &getSequenceString() const {
        return seq_;
    }

    const std::string &getQualityString() const {
        return qual_;
    }

    std::string getPhredQualityString(int offset = PHRED_OFFSET) const {
        std::string res = qual_;
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] = (char) (res[i] + offset);
        }
        return res;
    }

    const std::string &getName() const {
        return name_;
    }

    size_t size() const {
        return seq_.size();
    }

    char operator[](size_t i) const {
        VERIFY(is_nucl(seq_[i]));
        return dignucl(seq_[i]);
    }

    /**
      * trim read
      * @param ltrim first good base
      * @param rtrim last good base
      * @return whether there is anything left
      */
    bool trimLeftRight(int ltrim, int rtrim) {
        if (ltrim >= (int) seq_.size() || rtrim < 0 || rtrim < ltrim) {
            seq_ = "";
            qual_ = "";
            valid_ = false;
            return 0;
        }
        bool donesomething = false;
        if (ltrim > 0) {
            ltrim_ += ltrim;
            seq_.erase(0, ltrim);
            qual_.erase(0, ltrim);
            donesomething = true;
        }
        if (rtrim - ltrim + 1 < (int) seq_.size() && rtrim < (int) seq_.size() - ltrim - 1) {
            rtrim_ -= ((int) seq_.size() - (rtrim - ltrim + 1));
            seq_.erase(rtrim - ltrim + 1, std::string::npos);
            qual_.erase(rtrim - ltrim + 1, std::string::npos);
            donesomething = true;
        }
        if (donesomething) valid_ = updateValid();
        return true;
    }

    size_t trimNsAndBadQuality(int threshold) {
        int start = 0;
        for (; start < (int) seq_.size(); ++start) {
            if (seq_[start] != 'N' && (int) qual_[start] > threshold) break;
        }
        int end = 0;
        for (end = (int) seq_.size() - 1; end > -1; --end) {
            if (seq_[end] != 'N' && (int) qual_[end] > threshold) break;
        }
        if (!trimLeftRight(start, end)) return 0;
        else return seq_.size();
    }

    /**
     * @param k k as in k-mer
     * @param start start point
     * @return the first starting point of a valid k-mer >=start; return -1 if no such place exists
     */
    size_t firstValidKmer(size_t start, size_t k) const __attribute__ ((deprecated)) {
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
        return -1ULL;
    }

    void setSequence(const char *s, bool preserve_trimming = false) {
        seq_ = s;
        if (!preserve_trimming) {
            ltrim_ = 0;
            rtrim_ = initial_size_ = (int) seq_.size();
        }
        valid_ = updateValid();
    }

    void setQuality(const char *s, int offset = PHRED_OFFSET) {
        qual_ = s;
        for (size_t i = 0; i < qual_.size(); ++i) {
            qual_[i] = (char) (qual_[i] - offset);
        }
    }

    void setName(const char *s) {
        name_ = s;
    }

    Read()
            : valid_(false), ltrim_(0), rtrim_(0), initial_size_(0) {
        ;
    }

    Read(const std::string &name, const std::string &seq, const std::string &qual) :
            name_(name), seq_(seq), qual_(qual) {  // for test only!
        ltrim_ = 0;
        initial_size_ = rtrim_ = (int) seq_.size();
        valid_ = updateValid();
    }

    int ltrim() const { return ltrim_; }

    void set_ltrim(unsigned val) { ltrim_ = val; };

    int rtrim() const { return rtrim_; }

    int initial_size() const { return initial_size_; }

private:
    std::string name_;
    std::string seq_;
    std::string qual_;
    bool valid_;
    int ltrim_;
    int rtrim_;
    int initial_size_;

    friend class ireadstream;

    friend uint32_t TrimBadQuality(Read *, int);

    bool updateValid() const {
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

    void print(std::ostream &outf, int offset) const {
        outf << "@" << name_.c_str() << "\n";
        for (int i = 0; i < ltrim_; ++i) outf << "N";
        outf << seq_.c_str();
        for (int i = 0; i < initial_size_ - rtrim_; ++i) outf << "N";
        outf << "\n" << "+" << name_.c_str();
        if (ltrim_ > 0) outf << " ltrim=" << ltrim_;
        if (rtrim_ < initial_size_)
            outf << " rtrim=" << (initial_size_ - rtrim_);
        outf << "\n";
        char badq = (char) (offset + 2);
        for (int i = 0; i < ltrim_; ++i) outf << badq;
        outf << getPhredQualityString(offset).c_str();
        for (int i = 0; i < initial_size_ - rtrim_; ++i) outf << badq;
        outf << "\n";
    }
};

// todo: put this to *.cpp
//ostream& operator<<(ostream& os, const Read& read) {
//    return os << read.getSequenceString();
//}

#endif /* READ_HPP_ */
