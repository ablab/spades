//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
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

#include "parser.hpp"
#include "read_stream.hpp"
#include "single_read.hpp"

#include "utils/logger/logger.hpp"

#include <filesystem>

namespace io {

class FileReadStream {
public:
    typedef SingleRead ReadT;
    
    /*
     * Default constructor.
     *
     * @param filename The name of the file to be opened.
     * @param distance Doesn't have any sense here, but necessary for
     * wrappers.
     * @param offset The offset of the read quality.
     */
    explicit FileReadStream(const std::filesystem::path &filename,
                            FileReadFlags flags = FileReadFlags())
            : filename_(filename), flags_(flags), parser_(nullptr) {
        CHECK_FATAL_ERROR(exists(filename), "File " << filename << " doesn't exist or can't be read!");
        parser_.reset(SelectParser(filename_, flags_));
    }

    /*
     * Check whether the stream is opened.
     *
     * @return true of the stream is opened and false otherwise.
     */
    bool is_open() {
        if (!parser_)
            return false;

        return parser_->is_open();
    }

    /*
     * Check whether we've reached the end of stream.
     *
     * @return true if the end of stream is reached and false
     * otherwise.
     */
    bool eof() {
        if (!parser_)
            return true;

        return parser_->eof();
    }

    /*
     * Read SingleRead from stream.
     *
     * @param singleread The SingleRead that will store read data.
     *
     * @return Reference to this stream.
     */
    FileReadStream &operator>>(SingleRead &singleread) {
        if (parser_)
            (*parser_) >> singleread;

        return *this;
    }

    /*
     * Close the stream.
     */
    void close() {
        if (!parser_)
            return;

        parser_->close();
    }

    /*
     * Close the stream and open it again.
     */
    void reset() {
        if (!parser_)
            return;

        parser_->reset();
    }

private:
    /* @variable The name of the file which stream reads from. */
    std::filesystem::path filename_;
    /* @variable Flags */
    FileReadFlags flags_;
    /* @variable Internal stream that reads from file. */
    std::unique_ptr<Parser> parser_;

};

}
