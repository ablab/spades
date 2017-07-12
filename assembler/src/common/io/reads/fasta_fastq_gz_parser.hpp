//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/**
 * @file    fastqgz_parser.hpp
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
 * FastaFastqGzParser is the parser stream that reads data from .fastq.gz
 * files.
 */

#ifndef COMMON_IO_FASTAFASTQGZPARSER_HPP
#define COMMON_IO_FASTAFASTQGZPARSER_HPP

#include <zlib.h>
#include <string>
#include "kseq/kseq.h"
#include "utils/verify.hpp"
#include "single_read.hpp"
#include "io/reads/parser.hpp"
#include "sequence/quality.hpp"
#include "sequence/nucl.hpp"

namespace io {

namespace fastafastqgz {
// STEP 1: declare the type of file handler and the read() function
// Silence bogus gcc warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)
#pragma GCC diagnostic pop
}

class FastaFastqGzParser: public Parser {
public:
    /*
     * Default constructor.
     *
     * @param filename The name of the file to be opened.
     * @param offset The offset of the read quality.
     */
    FastaFastqGzParser(const std::string& filename, OffsetType offset_type =
            PhredOffset) :
            Parser(filename, offset_type), fp_(), seq_(NULL) {
        open();
    }

    /*
     * Default destructor.
     */
    /* virtual */
    ~FastaFastqGzParser() {
        close();
    }

    /*
     * Read SingleRead from stream.
     *
     * @param read The SingleRead that will store read data.
     *
     * @return Reference to this stream.
     */
    /* virtual */
    FastaFastqGzParser& operator>>(SingleRead& read) {
        if (!is_open_ || eof_) {
            return *this;
        }
        //todo offset_type_ should be used in future
        if (seq_->qual.s) {
            read = SingleRead(seq_->name.s, seq_->seq.s, seq_->qual.s, offset_type_);
        } else {
            read = SingleRead(seq_->name.s, seq_->seq.s);
//            size_t len = strlen(seq_->seq.s);
//            char* qual = (char*) malloc(len + 1);
//            char q = '\2' + 64;
//            for (size_t i = 0; i < len; ++i) {
//                qual[i] = q;
//            }
//            qual[len] = '\0';
//            read.SetAll(seq_->name.s, seq_->seq.s, qual, SolexaOffset);
//            free(qual);
        }
        ReadAhead();
        return *this;
    }

    /*
     * Close the stream.
     */
    /* virtual */
    void close() {
        if (is_open_) {
            // STEP 5: destroy seq
            fastafastqgz::kseq_destroy(seq_);
            // STEP 6: close the file handler
            gzclose(fp_);
            is_open_ = false;
            eof_ = true;
        }
    }

private:
    /*
     * @variable File that is associated with gzipped data file.
     */
    gzFile fp_;
    /*
     * @variable Data element that stores last SingleRead got from
     * stream.
     */
    fastafastqgz::kseq_t* seq_;

    /*
     * Open a stream.
     */
    /* virtual */
    void open() {
        // STEP 2: open the file handler
        fp_ = gzopen(filename_.c_str(), "r");
        if (!fp_) {
            is_open_ = false;
            return;
        }
        // STEP 3: initialize seq
        seq_ = fastafastqgz::kseq_init(fp_);
        eof_ = false;
        is_open_ = true;
        ReadAhead();
    }

    /*
     * Read next SingleRead from file.
     */
    void ReadAhead() {
        VERIFY(is_open_);
        VERIFY(!eof_);
        if (fastafastqgz::kseq_read(seq_) < 0) {
            eof_ = true;
        }
    }

    /*
     * Hidden copy constructor.
     */
    FastaFastqGzParser(const FastaFastqGzParser& parser);
    /*
     * Hidden assign operator.
     */
    void operator=(const FastaFastqGzParser& parser);
};

}

#endif /* COMMON_IO_FASTAFASTQGZPARSER_HPP */
