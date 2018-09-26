//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/verify.hpp"
#include "sequence/quality.hpp"
#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"
#include "sequence/sequence_tools.hpp"
#include "utils/stl_utils.hpp"

#include <string>

namespace io {

/*
* This enumerate contains offset type.
* UnknownOffset is equal to "offset = 0".
* PhredOffset is equal to "offset = 33".
* SolexaOffset is equal to "offset = 64".
*/
//todo change to enum class
enum OffsetType {
    UnknownOffset = 0,
    PhredOffset = 33,
    SolexaOffset = 64
};

//todo extract code about offset from here
typedef uint16_t SequenceOffsetT;

class SingleRead {
public:

    static std::string EmptyQuality(const std::string &seq) {
        return std::string(seq.size(), (char) PhredOffset);
    }

    SingleRead() :
            name_(""), seq_(""), qual_(""), left_offset_(0), right_offset_(0), valid_(false) {
    }

    SingleRead(const std::string &name, const std::string &seq,
               const std::string &qual, OffsetType offset,
               SequenceOffsetT left_offset = 0, SequenceOffsetT right_offset = 0) :
            name_(name), seq_(seq), qual_(qual), left_offset_(left_offset), right_offset_(right_offset) {
        Init();
        for (size_t i = 0; i < qual_.size(); ++i) {
            qual_[i] = (char) (qual_[i] - offset);
        }
    }

    SingleRead(const std::string &name, const std::string &seq,
               const std::string &qual,
               SequenceOffsetT left_offset = 0, SequenceOffsetT right_offset = 0) :
            name_(name), seq_(seq), qual_(qual), left_offset_(left_offset), right_offset_(right_offset) {
        Init();
    }

    SingleRead(const std::string &name, const std::string &seq,
               SequenceOffsetT left_offset = 0, SequenceOffsetT right_offset = 0) :
            name_(name), seq_(seq), qual_(EmptyQuality(seq_)), left_offset_(left_offset),
            right_offset_(right_offset) {
        Init();
    }

    bool IsValid() const {
        return valid_;
    }

    Sequence sequence(bool rc = false) const {
        VERIFY(valid_);
        return Sequence(seq_, rc);
    }

    Quality quality() const {
        VERIFY(valid_);
        return Quality(qual_);
    }

    const std::string &name() const {
        return name_;
    }

    size_t size() const {
        return seq_.size();
    }

    size_t nucl_count() const {
        return size();
    }

    const std::string &GetSequenceString() const {
        return seq_;
    }

    const std::string &GetQualityString() const {
        return qual_;
    }

    std::string GetPhredQualityString() const {
        int offset = PhredOffset;
        std::string res = qual_;
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] = (char) (res[i] + offset);
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

    SingleRead operator!() const {
        std::string new_name;
        if (name_.length() >= 3 && name_.substr(name_.length() - 3) == "_RC") {
            new_name = name_.substr(0, name_.length() - 3);
        } else {
            new_name = name_ + "_RC";
        }
        //        TODO make naming nicer
        //        if (name_ == "" || name_[0] != '!') {
        //            new_name = '!' + name_;
        //        } else {
        //            new_name = name_.substr(1, name_.length());
        //        }
        return SingleRead(new_name, ReverseComplement(seq_), Reverse(qual_), right_offset_, left_offset_);
    }

    SingleRead Substr(size_t from, size_t to) const {
        VERIFY(from <= to && to <= size());
        size_t len = to - from;
        if (len == size()) {
            return *this;
        }
        if (len == 0) {
            return SingleRead();
        }
        return SubstrStrict(from, to);
    }

    bool operator==(const SingleRead &singleread) const {
        return seq_ == singleread.seq_;
    }

    void ChangeName(const std::string &new_name) {
        name_ = new_name;
    }

    static bool IsValid(const std::string &seq) {
        for (size_t i = 0; i < seq.size(); ++i) {
            if (!is_nucl(seq[i])) {
                return false;
            }
        }
        return true;
    }

    SequenceOffsetT GetLeftOffset() const {
        return left_offset_;
    }

    SequenceOffsetT GetRightOffset() const {
        return right_offset_;
    }

    bool BinWrite(std::ostream &file, bool rc = false) const {
        sequence(rc).BinWrite(file);
        if (rc) {
            file.write((const char *) &right_offset_, sizeof(right_offset_));
            file.write((const char *) &left_offset_, sizeof(left_offset_));
        } else {
            file.write((const char *) &left_offset_, sizeof(left_offset_));
            file.write((const char *) &right_offset_, sizeof(right_offset_));
        }
        return !file.fail();
    }

    bool BinRead(std::istream &/*file*/) {
        VERIFY(false);
        return false;
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

    //Left and right offsets with respect to original sequence
    SequenceOffsetT left_offset_;

    SequenceOffsetT right_offset_;

    bool valid_;

    void Init() {
        VERIFY(seq_.size() == qual_.size());
        valid_ = SingleRead::IsValid(seq_);
    }

    SingleRead SubstrStrict(size_t from, size_t to) const {
        size_t len = to - from;
        //        return SingleRead(name_, seq_.substr(from, len), qual_.substr(from, len));
        //        TODO remove naming?
        std::string new_name;
        if (name_.length() >= 3 && name_.substr(name_.length() - 3) == "_RC") {
            new_name = name_.substr(0, name_.length() - 3) + "_SUBSTR(" + std::to_string(size() - to) + "," +
                       std::to_string(size() - from) + ")" + "_RC";
        } else {
            new_name = name_ + "_SUBSTR(" + std::to_string(from) + "," + std::to_string(to) + ")";
        }
        return SingleRead(new_name, seq_.substr(from, len), qual_.substr(from, len),
                          SequenceOffsetT(from + (size_t) left_offset_),
                          SequenceOffsetT(size() - to + (size_t) right_offset_));
    }


};

inline std::ostream &operator<<(std::ostream &os, const SingleRead &read) {
    os << "Single read name=" << read.name() << " sequence=" << read.GetSequenceString() << std::endl;
    return os;
}

class SingleReadSeq {

public:
    explicit SingleReadSeq(const Sequence &s,
                  SequenceOffsetT left_offset = 0, SequenceOffsetT right_offset = 0) :
            seq_(s), left_offset_(left_offset), right_offset_(right_offset) {
    }

    SingleReadSeq() : seq_(), left_offset_(0), right_offset_(0) {
    }

    bool BinRead(std::istream &file) {
        seq_.BinRead(file);
        file.read((char *) &left_offset_, sizeof(left_offset_));
        file.read((char *) &right_offset_, sizeof(right_offset_));
        return !file.fail();
    }

    bool BinWrite(std::ostream &file, bool rc = false) const {
        if (rc)
            (!seq_).BinWrite(file);
        else
            seq_.BinWrite(file);
        if (rc) {
            file.write((const char *) &right_offset_, sizeof(right_offset_));
            file.write((const char *) &left_offset_, sizeof(left_offset_));
        } else {
            file.write((const char *) &left_offset_, sizeof(left_offset_));
            file.write((const char *) &right_offset_, sizeof(right_offset_));
        }
        return !file.fail();
    }

    //    SingleReadSeq(std::istream& file): seq_(file, true) {
    //    }

    bool operator==(const SingleReadSeq &singleread) const {
        return seq_ == singleread.seq_;
    }

    const Sequence sequence() const {
        return seq_;
    }

    size_t size() const {
        return seq_.size();
    }

    size_t nucl_count() const {
        return size();
    }

    SingleReadSeq operator!() const {
        return SingleReadSeq(!seq_, right_offset_, left_offset_);
    }

    SequenceOffsetT GetLeftOffset() const {
        return left_offset_;
    }

    SequenceOffsetT GetRightOffset() const {
        return right_offset_;
    }

private:
    Sequence seq_;

    //Left and right offsets with respect to original sequence
    SequenceOffsetT left_offset_;

    SequenceOffsetT right_offset_;
};

inline std::ostream &operator<<(std::ostream &os, const SingleReadSeq &read) {
    os << "Single read sequence=" << read.sequence() << std::endl;
    return os;
}

}
