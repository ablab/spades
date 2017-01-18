//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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

inline std::string MakeContigId(size_t number, size_t length, const std::string& prefix = "NODE") {
    return prefix + "_" + ToString(number) + "_length_" + ToString(length);
}

inline std::string MakeContigId(size_t number, size_t length, double coverage, const std::string& prefix = "NODE") {
    return MakeContigId(number, length, prefix) + "_cov_" + ToString(coverage);
}

inline std::string MakeContigId(size_t number, size_t length, double coverage, size_t id, const std::string& prefix = "NODE") {
    return MakeContigId(number, length, coverage, prefix) + "_ID_" +  ToString(id);
}

inline std::string MakeRNAContigId(size_t number, size_t length, double coverage, size_t gene_id, size_t isoform_id, const std::string& prefix = "NODE") {
    return MakeContigId(number, length, coverage, prefix) + "_g" + ToString(gene_id)  + "_i" + ToString(isoform_id);
}

inline std::string MakeContigComponentId(size_t number, size_t length, double coverage, size_t component_id, const std::string& prefix = "NODE") {
    return MakeContigId(number, length, coverage, prefix)  + "_component_" + ToString(component_id);
}

class osequencestream {
protected:
    std::ofstream ofstream_;

    size_t id_;

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
    osequencestream(const std::string& filename): id_(1) {
            ofstream_.open(filename.c_str());
    }


    virtual ~osequencestream() {
        ofstream_.close();
    }

    virtual osequencestream& operator<<(const std::string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

    virtual osequencestream& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

    /**
     * Has different way of making headers
     * Doesn't increase counters, don't mix with other methods!
     */
    virtual osequencestream& operator<<(const SingleRead& read) {
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


class osequencestream_simple: public osequencestream {
protected:
    std::string header_;

    double cov_;

    virtual void write_header(const std::string& /*s*/) {
        ofstream_ << ">" << header_ << std::endl;
    }

public:
    osequencestream_simple(const std::string& filename)
            : osequencestream(filename), header_("") { }

    virtual ~osequencestream_simple() {
        ofstream_.close();
    }

    void set_header(const std::string &header) {
        header_ = header;
    }

    osequencestream_simple& operator<<(const std::string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

    osequencestream_simple& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

};

class osequencestream_with_id: public osequencestream {
protected:
    size_t uid_;

    double cov_;

    virtual void write_header(const std::string& s) {
        ofstream_ << ">" << GetId(s) << std::endl;
        id_++;
    }

public:
    osequencestream_with_id(const std::string& filename)
        : osequencestream(filename), uid_(0), cov_(0.0) { }

    virtual ~osequencestream_with_id() {
        ofstream_.close();
    }

    std::string GetId(const std::string& s) const {
        return MakeContigId(id_, s.size(), cov_, uid_);
    }

    void setCoverage(double c) {
        cov_ = c;
    }

    void setID(size_t uid) {
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

class osequencestream_with_manual_node_id: public osequencestream_with_id {
    bool is_id_set_;
    virtual void write_header(const std::string& s) {
        //for manual NODE ID setting osequencestream need to chech that node ID is really manually set
        if (!is_id_set_) {
            WARN ("NODE ID is not set manually, setting to 0");
            id_ = 0;
        }
        ofstream_ << ">" << MakeContigId(id_, s.size(), cov_, uid_) << std::endl;
        is_id_set_ = false;
    }

public:
//unfortunately constructor inheritance is supported only since g++4.8
    osequencestream_with_manual_node_id(const std::string& filename): osequencestream_with_id(filename) {
        is_id_set_ = false;
    }

    void setNodeID(int id) {
        id_ = id;
        is_id_set_ = true;
    }

    osequencestream_with_manual_node_id& operator<<(const std::string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

    osequencestream_with_manual_node_id& operator<<(const Sequence& seq) {
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
    osequencestream_for_fastg(const std::string& filename):
            osequencestream_with_id(filename) {
        id_ = 1;
    }

    virtual ~osequencestream_for_fastg() {
        ofstream_.close();
    }

    void set_header(const std::string& h) {
        header_=  h;
    }

    osequencestream_for_fastg& operator<<(const std::set<std::string>& s) {
        write_header(header_);
        if (s.size() > 0) {
            auto iter = s.begin();
            ofstream_ << ":" << *iter;
            ++iter;
            while (iter != s.end()) {
                ofstream_ << "," << *iter;
                ++iter;
            }
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
