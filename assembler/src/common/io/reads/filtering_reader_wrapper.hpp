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

#pragma once

#include <boost/optional.hpp>
#include "delegating_reader_wrapper.hpp"
#include "read_stream_vector.hpp"

namespace io {

template<typename ReadType>
class FilteringReaderWrapper : public DelegatingWrapper<ReadType> {
    typedef DelegatingWrapper<ReadType> base;
    //allows read to be modified
public:
    typedef std::function<bool (ReadType&)> FilterF;
//    static constexpr FilterF VALIDITY_FILTER = [] (ReadType& r) -> bool { return r.IsValid(); };

    static bool VALIDITY_FILTER(ReadType& r) {
        return r.IsValid();
    };

  /*
   * Default constructor.
   *
   * @param reader Reference to any other reader (child of IReader).
   */
    explicit FilteringReaderWrapper(typename base::ReadStreamPtrT reader_ptr,
                                    FilterF filter = VALIDITY_FILTER) :
            base(reader_ptr), filter_f_(filter), eof_(false) {
        StepForward();
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
     * Close the stream and open it again.
     */
    void reset() override {
        base::reset();
        eof_ = false;
        StepForward();
    }

private:

    FilterF filter_f_;

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
        while (!base::eof()) {
            base::operator>>(next_read_);

            if (filter_f_(next_read_)) {
                return;
            }
        }
        eof_ = true;
    }

};

template<class ReadType>
std::shared_ptr<ReadStream<ReadType>> FilteringWrap(std::shared_ptr<ReadStream<ReadType>> reader_ptr,
                                                    typename FilteringReaderWrapper<ReadType>::FilterF filter =
                                                            FilteringReaderWrapper<ReadType>::VALIDITY_FILTER) {
    return std::make_shared<FilteringReaderWrapper<ReadType>>(reader_ptr, filter);
}

template<class ReadType>
ReadStreamList<ReadType> FilteringWrap(const ReadStreamList<ReadType>& readers,
                                       typename FilteringReaderWrapper<ReadType>::FilterF filter =
                                            FilteringReaderWrapper<ReadType>::VALIDITY_FILTER) {
    ReadStreamList<ReadType> answer;
    for (size_t i = 0; i < readers.size(); ++i) {
        answer.push_back(FilteringWrap<ReadType>(readers.ptr_at(i), filter));
    }
    return answer;
}

}
