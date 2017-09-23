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

inline std::string MakeContigId(size_t number, const std::string& prefix = "NODE") {
    return prefix.empty() ? std::to_string(number) : (prefix + "_" + std::to_string(number));
}

inline std::string MakeContigId(size_t number, size_t length, const std::string& prefix = "NODE") {
    return MakeContigId(number, prefix) + "_length_" + std::to_string(length);
}

inline std::string MakeContigId(size_t number, size_t length, double coverage, const std::string& prefix = "NODE") {
    return MakeContigId(number, length, prefix) + "_cov_" + std::to_string(coverage);
}

inline std::string MakeContigId(size_t number, size_t length, double coverage, size_t id, const std::string& prefix = "NODE") {
    return MakeContigId(number, length, coverage, prefix) + "_ID_" +  std::to_string(id);
}

inline std::string MakeRNAContigId(size_t number, size_t length, double coverage, size_t gene_id, size_t isoform_id, const std::string& prefix = "NODE") {
    return MakeContigId(number, length, coverage, prefix) + "_g" + std::to_string(gene_id)  + "_i" + std::to_string(isoform_id);
}

inline std::string AddComponentId(const string& s, size_t component_id) {
    return s + "_component_" + std::to_string(component_id);
}

inline void WriteWrapped(const std::string &s, ostream &os, size_t max_width = 60) {
    size_t cur = 0;
    while (cur < s.size()) {
        os << s.substr(cur, max_width) << "\n";
        cur += max_width;
    }
}

class osequencestream {
protected:
    std::ofstream ofstream_;
    size_t id_;

    void write_str(const std::string& s) {
        WriteWrapped(s, ofstream_);
    }

    virtual void write_header(const std::string& s) {
        // Velvet format: NODE_1_length_24705_cov_358.255249
        ofstream_ << ">" << MakeContigId(id_++, s.size()) << std::endl;
    }

public:
    osequencestream(const std::string& filename):
            ofstream_(filename), id_(1) {
    }

    virtual ~osequencestream() {}

    osequencestream& operator<<(const std::string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

    osequencestream& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
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

    osequencestream_cov& operator<<(double coverage) {
        coverage_ = coverage;
        return *this;
    }

    using osequencestream::operator<<;

};

class OutputSequenceStream {
    std::ofstream ofstream_;
public:
    typedef SingleRead ReadT;

    OutputSequenceStream(const std::string& filename):
            ofstream_(filename) {
    }

    OutputSequenceStream& operator<<(const SingleRead& read) {
        ofstream_ << ">" << read.name() << "\n";
        WriteWrapped(read.GetSequenceString(), ofstream_);
        return *this;
    }
};

class PairedOutputSequenceStream {
    OutputSequenceStream os_l_;
    OutputSequenceStream os_r_;

public:
    typedef PairedRead ReadT;

    PairedOutputSequenceStream(const std::string& filename1,
                               const std::string &filename2) :
            os_l_(filename1),
            os_r_(filename2) {
    }

    PairedOutputSequenceStream& operator<<(const PairedRead& read) {
        os_l_ << read.first();
        os_r_ << read.second();
        return *this;
    }
};

//FIXME reduce code duplication
class OSingleReadStream {
    std::ofstream os_;

public:
    typedef SingleRead ReadT;

    OSingleReadStream(const std::string& fn) :
            os_(fn) {
    }

    OSingleReadStream& operator<<(const SingleRead& read) {
        os_ << "@" << read.name() << std::endl;
        os_ << read.GetSequenceString() << std::endl;
        os_ << "+" << std::endl;
        os_ << read.GetPhredQualityString() << std::endl;
        return *this;
    }

    void close() {
        os_.close();
    }
};

class OPairedReadStream {
    OSingleReadStream l_os_;
    OSingleReadStream r_os_;

public:
    typedef PairedRead ReadT;

    OPairedReadStream(const std::string& l_fn, const std::string& r_fn) :
            l_os_(l_fn), r_os_(r_fn) {
    }

    OPairedReadStream& operator<<(const PairedRead& read) {
        l_os_ << read.first();
        r_os_ << read.second();
        return *this;
    }

    void close() {
        l_os_.close();
        r_os_.close();
    }
};

}
