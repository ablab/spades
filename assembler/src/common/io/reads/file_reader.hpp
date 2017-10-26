//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/**

* Reader<SingleRead> is the very base class that reads from one file
* through Parser object.
* Reader<PairedRead> is the class that reads data from two input
* files and gets PairedReads using this data and distance information.
*/

#pragma once

#include "ireader.hpp"
#include "single_read.hpp"
#include "parser.hpp"
#include "utils/filesystem/path_helper.hpp"

namespace io {

class FileReadStream : public ReadStream<SingleRead> {
public:
    /*
     * Default constructor.
     *
     * @param filename The name of the file to be opened.
     * @param distance Doesn't have any sense here, but necessary for
     * wrappers.
     * @param offset The offset of the read quality.
     */
    explicit FileReadStream(const std::string &filename,
                            OffsetType offset_type = PhredOffset)
            : filename_(filename), offset_type_(offset_type), parser_(NULL) {
        fs::CheckFileExistenceFATAL(filename_);
        parser_ = SelectParser(filename_, offset_type_);
    }

    /*
     * Default destructor.
     */
    /* virtual */ ~FileReadStream() {
        close();
        delete parser_;
    }

    /*
     * Check whether the stream is opened.
     *
     * @return true of the stream is opened and false otherwise.
     */
    /* virtual */ bool is_open() {
        if (parser_ != NULL) {
            return parser_->is_open();
        } else {
            return false;
        }
    }

    /*
     * Check whether we've reached the end of stream.
     *
     * @return true if the end of stream is reached and false
     * otherwise.
     */
    /* virtual */ bool eof() {
        if (parser_ != NULL) {
            return parser_->eof();
        } else {
            return true;
        }
    }

    /*
     * Read SingleRead from stream.
     *
     * @param singleread The SingleRead that will store read data.
     *
     * @return Reference to this stream.
     */
    /* virtual */ FileReadStream &operator>>(SingleRead &singleread) {
        if (parser_ != NULL) {
            (*parser_) >> singleread;
        }
        return *this;
    }

    /*
     * Close the stream.
     */
    /* virtual */ void close() {
        if (parser_ != NULL) {
            parser_->close();
        }
    }

    /*
     * Close the stream and open it again.
     */
    /* virtual */ void reset() {
        if (parser_ != NULL) {
            parser_->reset();
        }
    }

private:
    /*
     * @variable The name of the file which stream reads from.
     */
    std::string filename_;
    /*
     * @variable Quality offset type.
     */
    OffsetType offset_type_;
    /*
     * @variable Internal stream that reads from file.
     */
    Parser *parser_;

};

}
