//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "file_read_flags.hpp"
#include "sequence/quality.hpp"
#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"
#include "sequence/sequence_tools.hpp"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"

#include <string>

namespace io {

//todo extract code about offset from here
typedef uint16_t SequenceOffsetT;

class SingleRead {
public:

// Silence bogus GCC warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"

    SingleRead()
            : name_{}, comment_{},
              seq_{}, qual_{},
              left_offset_(0), right_offset_(0), valid_(false) { }

    SingleRead(std::string name, std::string comment,
               std::string seq,
               std::string qual, OffsetType offset,
               SequenceOffsetT left_offset = 0, SequenceOffsetT right_offset = 0,
               bool validate = true)
            : name_{std::move(name)}, comment_{std::move(comment)},
              seq_{std::move(seq)}, qual_{std::move(qual)},
              left_offset_(left_offset), right_offset_(right_offset),
              valid_(false) {
        Init(validate);
        for (size_t i = 0; i < qual_.size(); ++i) {
            qual_[i] = (char) (qual_[i] - offset);
        }
    }

    SingleRead(std::string name, std::string comment,
               std::string seq,
               std::string qual,
               SequenceOffsetT left_offset = 0, SequenceOffsetT right_offset = 0,
               bool validate = true)
            : name_{std::move(name)}, comment_{std::move(comment)},
              seq_{std::move(seq)}, qual_{std::move(qual)},
              left_offset_(left_offset), right_offset_(right_offset), valid_(false) {
        Init(validate);
    }

    SingleRead(std::string name, std::string comment,
               std::string seq,
               SequenceOffsetT left_offset = 0, SequenceOffsetT right_offset = 0,
               bool validate = true)
            : name_{std::move(name)}, comment_{std::move(comment)},
              seq_{std::move(seq)}, qual_{},
              left_offset_(left_offset), right_offset_(right_offset), valid_(false) {
        Init(validate);
    }

    SingleRead(std::string name,
               std::string seq,
               SequenceOffsetT left_offset = 0, SequenceOffsetT right_offset = 0,
               bool validate = true)
            : name_{std::move(name)}, comment_{},
              seq_{std::move(seq)}, qual_{},
              left_offset_(left_offset), right_offset_(right_offset), valid_(false) {
        Init(validate);
    }

    SingleRead(std::string seq,
               SequenceOffsetT left_offset = 0, SequenceOffsetT right_offset = 0,
               bool validate = true)
            : name_{}, comment_{},
              seq_{std::move(seq)}, qual_{},
              left_offset_(left_offset), right_offset_(right_offset), valid_(false) {
        Init(validate);
    }

#pragma GCC diagnostic pop

    bool IsValid() const {
        return valid_;
    }

    void validate() {
        CHECK_FATAL_ERROR(!qual_.size() || seq_.size() == qual_.size(),
                     "Invalid read: length of sequence should equal to length of quality line");
        valid_ = SingleRead::IsValid(seq_);
    }

    Sequence sequence(bool rc = false) const {
        return Sequence(seq_, rc);
    }

    Quality quality() const {
        VERIFY(valid_);
        return Quality(qual_);
    }

    const std::string &name() const {
        return name_;
    }

    const std::string &comment() const {
        return comment_;
    }

    size_t size() const {
        return seq_.size();
    }

    size_t qual_size() const {
        return qual_.size();
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
        } else if (name_.size()) {
            new_name = name_ + "_RC";
        }

        // Decide on ctor to call
        if (qual_.size())
            return SingleRead(new_name, comment_, ReverseComplement(seq_), Reverse(qual_), right_offset_, left_offset_);
        else if (new_name.size())
            return SingleRead(new_name, comment_, ReverseComplement(seq_), right_offset_, left_offset_);
        else
            return SingleRead(ReverseComplement(seq_), right_offset_, left_offset_);
    }

    SingleRead Substr(size_t from, size_t to, bool validate = true) const {
        VERIFY(from <= to && to <= size());
        size_t len = to - from;
        if (len == size())
            return *this;

        if (len == 0)
            return SingleRead();

        return SubstrStrict(from, to, validate);
    }

    bool operator==(const SingleRead &singleread) const {
        return seq_ == singleread.seq_;
    }

    void ChangeName(const std::string &new_name) {
        name_ = new_name;
    }

    static bool IsValid(const std::string &seq) {
        size_t sz = seq.size();
        for (size_t i = 0; i < sz; ++i) {
            if (!is_nucl(seq[i]))
                return false;
        }
        return true;
    }

    SequenceOffsetT GetLeftOffset() const {
        return left_offset_;
    }

    SequenceOffsetT GetRightOffset() const {
        return right_offset_;
    }

    bool BinWrite(std::ostream &file, bool rc = false, uint64_t tag = 0) const {
        sequence(rc).BinWrite(file);
        SequenceOffsetT left_offset = left_offset_, right_offset = right_offset_;
        if (rc) {
            file.write((const char *) &right_offset, sizeof(right_offset));
            file.write((const char *) &left_offset, sizeof(left_offset));
        } else {
            file.write((const char *) &left_offset, sizeof(left_offset));
            file.write((const char *) &right_offset, sizeof(right_offset));
        }

        file.write((const char *) &tag, sizeof(tag));
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
     * @variable The comment of SingleRead in input file.
     */
    std::string comment_;
    /*
     * @variable The sequence of nucleotides.
     */
    std::string seq_;
    /*
     * @variable The quality of SingleRead.
     */
    std::string qual_;

    // Left and right offsets with respect to original sequence
    SequenceOffsetT left_offset_  : 15;
    SequenceOffsetT right_offset_ : 15;
    bool valid_                   : 1;

    void Init(bool v) {
        if (v)
            validate();
    }

    SingleRead SubstrStrict(size_t from, size_t to, bool validate = true) const {
        size_t len = to - from;
        //        return SingleRead(name_, seq_.substr(from, len), qual_.substr(from, len));
        //        TODO remove naming?
        std::string new_name;
        if (name_.length() >= 3 && name_.substr(name_.length() - 3) == "_RC") {
            new_name = name_.substr(0, name_.length() - 3) + "_SUBSTR(" + std::to_string(size() - to) + "," +
                       std::to_string(size() - from) + ")" + "_RC";
        } else if (name_.length()) {
            new_name = name_ + "_SUBSTR(" + std::to_string(from) + "," + std::to_string(to) + ")";
        }

        if (qual_.length())
            return SingleRead(new_name, comment_, seq_.substr(from, len), qual_.substr(from, len),
                              SequenceOffsetT(from + (size_t) left_offset_),
                              SequenceOffsetT(size() - to + (size_t) right_offset_),
                              validate);
        else if (new_name.length())
            return SingleRead(new_name, comment_, seq_.substr(from, len),
                              SequenceOffsetT(from + (size_t) left_offset_),
                              SequenceOffsetT(size() - to + (size_t) right_offset_),
                              validate);
        else if (comment_.length())
            return SingleRead("", comment_, seq_.substr(from, len),
                              SequenceOffsetT(from + (size_t) left_offset_),
                              SequenceOffsetT(size() - to + (size_t) right_offset_),
                              validate);
        else
            return SingleRead(seq_.substr(from, len),
                              SequenceOffsetT(from + (size_t) left_offset_),
                              SequenceOffsetT(size() - to + (size_t) right_offset_),
                              validate);
    }
};

inline std::ostream &operator<<(std::ostream &os, const SingleRead &read) {
    os << "Single read name=" << (read.name().length() ? read.name() : "(empty)")
       << " comment= " << (read.comment().length() ? read.comment() : "(empty)")
       << " sequence=" << read.GetSequenceString() << std::endl;
    return os;
}

class SingleReadSeq {
public:
    explicit SingleReadSeq(const Sequence &s,
                           SequenceOffsetT left_offset = 0, SequenceOffsetT right_offset = 0,
                           uint64_t tag = -1ULL)
            : seq_(s), left_offset_(left_offset), right_offset_(right_offset), tag_(tag) {  }

    SingleReadSeq()
            : seq_(), left_offset_(0), right_offset_(0), tag_(-1ULL) {}

    bool BinRead(std::istream &file) {
        seq_.BinRead(file);
        file.read((char *) &left_offset_, sizeof(left_offset_));
        file.read((char *) &right_offset_, sizeof(right_offset_));
        file.read((char *) &tag_, sizeof(tag_));
        return !file.fail();
    }

    bool BinWrite(std::ostream &file, bool rc = false) const {
        if (rc) {
            (!seq_).BinWrite(file);
            file.write((const char *) &right_offset_, sizeof(right_offset_));
            file.write((const char *) &left_offset_, sizeof(left_offset_));
            file.write((const char *) &tag_, sizeof(tag_));
        } else {
            seq_.BinWrite(file);
            file.write((const char *) &left_offset_, sizeof(left_offset_));
            file.write((const char *) &right_offset_, sizeof(right_offset_));
            file.write((const char *) &tag_, sizeof(tag_));
        }
        return !file.fail();
    }

    bool BinWrite(std::ostream &file, bool rc, uint64_t tag) const {
        if (rc) {
            (!seq_).BinWrite(file);
            file.write((const char *) &right_offset_, sizeof(right_offset_));
            file.write((const char *) &left_offset_, sizeof(left_offset_));
            file.write((const char *) &tag, sizeof(tag));
        } else {
            seq_.BinWrite(file);
            file.write((const char *) &left_offset_, sizeof(left_offset_));
            file.write((const char *) &right_offset_, sizeof(right_offset_));
            file.write((const char *) &tag, sizeof(tag));
        }
        return !file.fail();
    }

    bool operator==(const SingleReadSeq &singleread) const {
        return seq_ == singleread.seq_;
    }

    Sequence sequence() const {
        return seq_;
    }

    size_t size() const {
        return seq_.size();
    }

    size_t nucl_count() const {
        return size();
    }

    uint64_t tag() const {
        return tag_;
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

    // Left and right offsets with respect to original sequence
    SequenceOffsetT left_offset_;
    SequenceOffsetT right_offset_;

    uint64_t tag_;
};

inline std::ostream &operator<<(std::ostream &os, const SingleReadSeq &read) {
    os << "Single read sequence=" << read.sequence() << std::endl;
    return os;
}

}
