//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/**
 * @file    multifile_reader.hpp
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
 * MultifileReader is the stream that gets data from number of files,
 * given in a constructor.
 */

#ifndef COMMON_IO_MULTIFILEREADER_HPP_
#define COMMON_IO_MULTIFILEREADER_HPP_

#include <vector>
#include "io/ireader.hpp"
#include "io/reader.hpp"

namespace io {

template<typename ReadType>
class MultifileReader: public IReader<ReadType> {
public:
	/*
	 * Default constructor.
	 *
	 * @param filenames The names of the files to be opened. Wrapper
	 * tries to open all the files of this list and proceeds only those
	 * of them which present. The names of non-existing files are
	 * ignored.
	 * @param distance Distance between parts of PairedReads (or useless
	 * parameter when we work with SingleReads).
	 * @param offset The offset of the read quality.
	 */
	MultifileReader(const vector<IReader<ReadType>*>& readers,
			bool destroy_readers = false) : /*filenames_(filenames), */
			current_reader_index_(0), destroy_readers_(destroy_readers) {
		for (size_t i = 0; i < readers.size(); ++i) {
			VERIFY(readers[i]->is_open());
			readers_.push_back(readers[i]);
		}
	}

	MultifileReader(IReader<ReadType>& reader_1, IReader<ReadType>& reader_2,
			bool destroy_readers = false) : /*filenames_(filenames), */
			current_reader_index_(0), destroy_readers_(destroy_readers) {
		VERIFY(reader_1.is_open() && reader_2.is_open());
		readers_.push_back(&reader_1);
		readers_.push_back(&reader_2);
	}

	/*
	 * Default destructor.
	 */
	/*virtual*/
	~MultifileReader() {
		if (destroy_readers_) {
			close();
			for (size_t i = 0; i < readers_.size(); ++i) {
				delete readers_[i];
			}
		}
	}

	/*
	 * Check whether the stream is opened.
	 *
	 * @return true if the stream is opened and false otherwise.
	 */
	/* virtual */
	bool is_open() {
		return (readers_.size() > 0) && readers_[0]->is_open();
	}

	/*
	 * Check whether we've reached the end of stream.
	 *
	 * @return true if the end of the stream is reached and false
	 * otherwise.
	 */
	/* virtual */
	bool eof() {
		while ((current_reader_index_ < readers_.size()) && readers_[current_reader_index_]->eof()) {
			++current_reader_index_;
		}
		return current_reader_index_ == readers_.size();
	}

	/*
	 * Read SingleRead or PairedRead from stream (according to ReadType).
	 *
	 * @param read The SingleRead or PairedRead that will store read
	 * data.
	 *
	 * @return Reference to this stream.
	 */
	/* virtual */
	MultifileReader& operator>>(ReadType& read) {
		if (!eof()) {
			(*readers_[current_reader_index_]) >> read;
		}
		return (*this);
	}

	/*
	 * Close the stream.
	 */
	/* virtual */
	void close() {
		for (size_t i = 0; i < readers_.size(); ++i) {
			readers_[i]->close();
		}
	}

	/*
	 * Close the stream and open it again.
	 */
	/* virtual */
	void reset() {
		for (size_t i = 0; i < readers_.size(); ++i) {
			readers_[i]->reset();
		}
		current_reader_index_ = 0;
	}

	ReadStat get_stat() const {
	    ReadStat stat;

        for (size_t i = 0; i < readers_.size(); ++i) {
            stat.merge(readers_[i]->get_stat());
        }
	    return stat;
	}

private:
	/*
	 * @variable Internal stream readers.
	 */
	vector<IReader<ReadType>*> readers_;
	/*
	 * @variable The index of the file that is currently read from.
	 */
	size_t current_reader_index_;

	bool destroy_readers_;
}
;

}

#endif /* COMMON_IO_MULTIFILEREADER_HPP_ */
