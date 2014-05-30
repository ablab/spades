//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * oreadstream.hpp
 *
 *  Created on: 23.06.2011
 *      Author: vyahhi
 */

#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "single_read.hpp"
#include "paired_read.hpp"

namespace io {

inline std::string MakeContigId(int number, size_t length) {
    return  "NODE_" + ToString(number) + "_length_" + ToString(length);
}

inline std::string MakeContigId(int number, size_t length, double coverage) {
    return "NODE_" + ToString(number)  + "_length_" + ToString(length) + "_cov_" + ToString(coverage);
}

inline std::string MakeContigId(int number, size_t length, double coverage, size_t id) {
    return "NODE_" + ToString(number)  + "_length_" + ToString(length) + "_cov_" + ToString(coverage)  + "_ID_" +  ToString(id);
}

class osequencestream {
protected:
    std::ofstream ofstream_;

    int id_;

    void write_str(const std::string& s) {
        size_t cur = 0;
        while (cur < s.size()) {
            ofstream_ << s.substr(cur, 60) << std::endl;
            cur += 60;
        }
    }

    virtual void write_header(const std::string& s) {
        // Velvet format: NODE_1_length_24705_cov_358.255249
        ofstream_ << ">" << MakeContigId(id_++, s.size()) << std::endl;
    }

public:
    osequencestream(const std::string& filename): id_(0) {
        ofstream_.open(filename.c_str());
    }

    virtual ~osequencestream() {
        ofstream_.close();
    }

    osequencestream& operator<<(const std::string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

    osequencestream& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

    /**
     * Has different way of making headers
     * Doesn't increase counters, don't mix with other methods!
     */
    osequencestream& operator<<(const SingleRead& read) {
        ofstream_ << ">" << read.name() << std::endl;
        size_t cur = 0;
        std::string s = read.GetSequenceString();
        while (cur < s.size()) {
            ofstream_ << s.substr(cur, 60) << std::endl;
            cur += 60;
        }
        return *this;
    }
};

class PairedOutputSequenceStream {
protected:
    std::ofstream ofstreaml_;
    std::ofstream ofstreamr_;

  static void write(const SingleRead& read, std::ofstream& stream) {
    stream << ">" << read.name() << std::endl;
    size_t cur = 0;
    std::string s = read.GetSequenceString();
    while (cur < s.size()) {
      stream << s.substr(cur, 60) << std::endl;
      cur += 60;
    }
  }

public:
    PairedOutputSequenceStream(const std::string& filename1, const std::string &filename2) {
      ofstreaml_.open(filename1);
      ofstreamr_.open(filename2);
    }

    virtual ~PairedOutputSequenceStream() {
        ofstreaml_.close();
        ofstreamr_.close();
    }

    PairedOutputSequenceStream& operator<<(const PairedRead& read) {
        write(read.first(), ofstreaml_);
        write(read.second(), ofstreamr_);
        return *this;
    }
};


class osequencestream_cov: public osequencestream {
protected:
    double coverage_;

    virtual void write_header(const std::string& s) {
        // Velvet format: NODE_1_length_24705_cov_358.255249
        ofstream_ << ">" << MakeContigId(id_++, s.size(), coverage_) << std::endl;
    }


public:
    osequencestream_cov(const std::string& filename)
        : osequencestream(filename), coverage_(0.) { }

    virtual ~osequencestream_cov() {
        ofstream_.close();
    }

    osequencestream_cov& operator<<(double coverage) {
        coverage_ = coverage;
        return *this;
    }

    osequencestream_cov& operator<<(const std::string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

    osequencestream_cov& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

};


class osequencestream_with_id: public osequencestream {
protected:
    int uid_;

    double cov_;

  virtual void write_header(const std::string& s) {
        ofstream_ << ">" << MakeContigId(++id_, s.size(), cov_, uid_) << std::endl;
    }

public:
    osequencestream_with_id(const std::string& filename)
        : osequencestream(filename), uid_(0), cov_(0.0) { }

    virtual ~osequencestream_with_id() {
        ofstream_.close();
    }

    void setCoverage(double c) {
        cov_ = c;
    }

    void setID(int uid) {
        uid_ = uid;
    }

    osequencestream_with_id& operator<<(const std::string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }
    osequencestream_with_id& operator<<(double coverage) {
        cov_ = coverage;
        return *this;
    }

    osequencestream_with_id& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

};


class osequencestream_with_data_for_scaffold: public osequencestream_with_id  {
protected:
    std::ofstream scstream_;

    virtual void write_header(const std::string& s) {
        scstream_ << id_ << "\tNODE_" << id_ << "\t" << s.size() << "\t" << (int) round(cov_) << std::endl;
        ofstream_ << ">" << MakeContigId(id_++, s.size(), cov_, uid_) << std::endl;
    }

public:
    osequencestream_with_data_for_scaffold(const std::string& filename): osequencestream_with_id(filename) {
        id_ = 1;
        std::string sc_filename = filename + ".info";
        scstream_.open(sc_filename.c_str());
    }

    virtual ~osequencestream_with_data_for_scaffold() {
        ofstream_.close();
        scstream_.close();
    }

    osequencestream_with_data_for_scaffold& operator<<(const std::string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

    osequencestream_with_data_for_scaffold& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

};

class osequencestream_for_fastg: public osequencestream_with_id  {
protected:
    std::string header_;

    virtual void write_header(const std::string& s) {
        ofstream_ << ">" << s;
    }

public:
    osequencestream_for_fastg(const std::string& filename): osequencestream_with_id(filename) {
        id_ = 1;
    }

    virtual ~osequencestream_for_fastg() {
        ofstream_.close();
    }

    void set_header(const std::string& h) {
        header_=  h;
    }

    osequencestream_for_fastg& operator<<(const std::vector<std::string>& v) {
        write_header(header_);
        for (size_t i = 0; i < v.size(); ++i) {
            ofstream_ << ":" << v[i];
        }
        ofstream_ << ";" << std::endl;
        return *this;
    }

    osequencestream_for_fastg& operator<<(const std::string& s) {
        write_str(s);
        return *this;
    }

    osequencestream_for_fastg& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

};

}
