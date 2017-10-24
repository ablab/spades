//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/**
 * @file    filtering_reader_wrapper.hpp
 * @author  Sergey Nurk
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
 * FilteringReaderWrapper is the class-wrapper that gets only valid
 * reads. 
 */

#ifndef COMMON_IO_FILTERINGREADERWRAPPER_HPP_
#define COMMON_IO_FILTERINGREADERWRAPPER_HPP_

#include "io/ireader.hpp"

namespace io {

template<typename ReadType>
class FilteringReaderWrapper: public IReader<ReadType> {
public:
  /*
   * Default constructor.
   *
   * @param reader Reference to any other reader (child of IReader).
   */
    explicit FilteringReaderWrapper(IReader<ReadType>& reader) :
            reader_(reader), eof_(false) {
        StepForward();
    }

  /* 
   * Default destructor.
   */
    /* virtual */ ~FilteringReaderWrapper() {
        close();
    }

  /* 
   * Check whether the stream is opened.
   *
   * @return true of the stream is opened and false otherwise.
   */
    bool is_open() override {
        return reader_.is_open();
    }

  /* 
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of stream is reached and false
   * otherwise.
   */
    bool eof() override {
        return eof_;
    }

  /*
   * Read SingleRead or PairedRead from stream (according to ReadType).
   *
   * @param read The SingleRead or PairedRead that will store read
   * data.
   *
   * @return Reference to this stream.
   */
    FilteringReaderWrapper& operator>>(ReadType& read) override {
        read = next_read_;
        StepForward();
        return *this;
    }

    /*
     * Close the stream.
     */
    void close() override {
        reader_.close();
    }

    /*
     * Close the stream and open it again.
     */
    void reset() override {
        reader_.reset();
        eof_ = false;
        StepForward();
    }

private:
  /*
   * @variable Internal stream readers.
   */
    IReader<ReadType>& reader_;
  /*
   * @variable Flag that shows whether the end of stream reached.
   */
    bool eof_;
  /*
   * @variable Next read to be outputted by stream.
   */
    ReadType next_read_;

  /*
   * Read next valid read in the stream.
   */
    void StepForward() {
        while (!reader_.eof()) {
            reader_ >> next_read_;
            if (next_read_.IsValid()) {
                return;
            }
        }
        eof_ = true;
    }

    /*
     * Hidden copy constructor.
     */
    explicit FilteringReaderWrapper(
            const FilteringReaderWrapper<ReadType>& reader);
    /*
     * Hidden assign operator.
     */
    void operator=(const FilteringReaderWrapper<ReadType>& reader);
};

}

#endif /* COMMON_IO_FILTERINGREADERWRAPPER_HPP_ */
