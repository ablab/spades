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
            : filename_(filename), offset_type_(offset_type), parser_(nullptr) {
        fs::CheckFileExistenceFATAL(filename_);
        parser_.reset(SelectParser(filename_, offset_type_));
    }

    /*
     * Default destructor.
     */
    ~FileReadStream() {
        close();
    }

    /*
     * Check whether the stream is opened.
     *
     * @return true of the stream is opened and false otherwise.
     */
    bool is_open() override {
        if (parser_)
            return parser_->is_open();

        return false;
    }

    /*
     * Check whether we've reached the end of stream.
     *
     * @return true if the end of stream is reached and false
     * otherwise.
     */
    bool eof() override {
        if (parser_)
            return parser_->eof();

        return true;
    }

    /*
     * Read SingleRead from stream.
     *
     * @param singleread The SingleRead that will store read data.
     *
     * @return Reference to this stream.
     */
    FileReadStream &operator>>(SingleRead &singleread) override {
        if (parser_)
            (*parser_) >> singleread;

        return *this;
    }

    /*
     * Close the stream.
     */
    void close() override {
        if (parser_)
            parser_->close();
    }

    /*
     * Close the stream and open it again.
     */
    void reset() override {
        if (parser_)
            parser_->reset();
    }

private:
    /* @variable The name of the file which stream reads from. */
    std::string filename_;
    /* @variable Quality offset type. */
    OffsetType offset_type_;
    /* @variable Internal stream that reads from file. */
    std::unique_ptr<Parser> parser_;

};

}
