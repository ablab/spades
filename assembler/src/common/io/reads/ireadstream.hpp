//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
* ifastqstream.hpp
*
*  Created on: 03.03.2011
*      Author: vyahhi
*/

#ifndef IREADSTREAM_HPP_
#define IREADSTREAM_HPP_

#include "kseq/kseq.h"
#include <zlib.h>
#include "utils/verify.hpp"
#include "read.hpp"
#include "sequence/nucl.hpp"

// Silence bogus gcc warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)
#pragma GCC diagnostic pop

/*
* Read name, seq and qual strings from FASTQ data (one by one)
*/
//fixme deprecated!!! remove usages!
class ireadstream {

public:
typedef Read ReadT;

ireadstream(const std::string &filename) : filename_(filename), offset_(Read::PHRED_OFFSET)  {
    is_open_ = open(filename);
}

ireadstream(const std::string &filename, int offset) : filename_(filename), offset_(offset) {
    is_open_ = open(filename);
}

virtual ~ireadstream() {
    close();
}

bool is_open() const {
    return is_open_;
}

bool eof() const {
    return eof_;
}

static std::vector <Read> *readAll(std::string filename, int cnt = -1) {
    ireadstream irs(filename);
    VERIFY(irs.is_open());
    std::vector <Read> *res = new std::vector<Read>();
    Read r;
    while (cnt-- && irs.is_open() && !irs.eof()) {
        irs >> r;
        if (!r.isValid()) {
            cnt++;
            continue;
        }
        res->push_back(r);
    }
    irs.close();
    return res;
}

static void readAllNoValidation(std::vector <Read> *res, std::string filename, uint64_t *totalsize,
                                int qvoffset = Read::PHRED_OFFSET, int trim_quality = -1, int cnt = -1) {
    ireadstream irs(filename, qvoffset);
    VERIFY(irs.is_open());
    *totalsize = 0;
    Read r;
    while (cnt-- && irs.is_open() && !irs.eof()) {
        irs >> r;
        size_t read_size = r.trimNsAndBadQuality(trim_quality);
        res->push_back(r);
        *totalsize += read_size;
    }
    irs.close();
}

ireadstream &operator>>(Read &r) {
    VERIFY(is_open());
    VERIFY(!eof());
    if (!is_open() || eof()) {
        return *this;
    }
    r.setName(seq_->name.s);
    if (seq_->qual.s) {
        r.setQuality(seq_->qual.s, offset_);
    }
    r.setSequence(seq_->seq.s);
    read_ahead(); // make actual read for the next result
    return *this;
}

void close() {
    if (is_open()) {
        kseq_destroy(seq_); // STEP 5: destroy seq
        gzclose(fp_); // STEP 6: close the file handler
        is_open_ = false;
    }
}

void reset() {
    close();
    open(filename_);
}

private:
std::string filename_;
gzFile fp_;
kseq_t *seq_;
bool is_open_;
bool eof_;
int offset_;

/*
 * open i's file with FASTQ reads,
 * return true if it opened file, false otherwise
 */
bool open(std::string filename) {
    fp_ = gzopen(filename.c_str(), "r"); // STEP 2: open the file handler
    if (!fp_) {
        return false;
    }
    is_open_ = true;
    seq_ = kseq_init(fp_); // STEP 3: initialize seq
    eof_ = false;
    read_ahead();
    return true;
}

void read_ahead() {
    VERIFY(is_open());
    VERIFY(!eof());
    if (kseq_read(seq_) < 0) {
        eof_ = true;
    }
}
};

//return -1 if failed to determine offset
inline int determine_offset(const std::string &filename) {
ireadstream stream(filename, 0);
size_t count = 0;
Read r;
while (!stream.eof() && count++ < 10000) {
    stream >> r;
    std::string q_str = r.getQualityString();
    for (size_t i = 0; i < q_str.size(); ++i) {
        int q_val = q_str[i];
        if (q_val < ';')
            return 33;
        if (q_val > 'K')
            return 64;
    }
}
return -1;
}

#endif /* IREADSTREAM_HPP_ */
