//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * binary_io.hpp
 *
 *  Created on: Apr 12, 2012
 *      Author: andrey
 */
#pragma once

#include <fstream>

#include "utils/verify.hpp"
#include "ireader.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"
#include "orientation.hpp"
#include "pipeline/library.hpp"

namespace io {

template<class Read>
class ReadBinaryWriter {
    bool rc_;

public:

    ReadBinaryWriter(bool rc = false) : rc_(rc) {
    }

    bool Write(std::ostream& file, const Read& r) const {
        return r.BinWrite(file, rc_);
    }
};

template<class Read>
class PairedReadBinaryWriter {
    bool rc1_;
    bool rc2_;

public:
    PairedReadBinaryWriter(LibraryOrientation orientation = LibraryOrientation::Undefined) {
        std::tie(rc1_, rc2_) = GetRCFlags(orientation);
    }

    bool Write(std::ostream& file, const Read& r) const {
        return r.BinWrite(file, rc1_, rc2_);
    }
};

class BinaryWriter {
    const std::string file_name_prefix_;
    size_t file_num_;
    std::vector<std::ofstream*> file_ds_;
    size_t buf_size_;

    template<class Writer, class Read>
    void FlushBuffer(const std::vector<Read>& buffer, const Writer& read_writer, std::ostream& file) {
        for (const Read &r : buffer) {
            read_writer.Write(file, r);
        }
    }

    template<class Writer, class Read>
    ReadStreamStat ToBinary(const Writer &writer, io::ReadStream<Read> &stream, size_t buf_size) {
        size_t buffer_reads = buf_size / (sizeof (Read) * 4);
        size_t reads_to_flush = buffer_reads * file_num_;

        std::vector<std::vector<Read>> buf(file_num_, std::vector<Read>(buffer_reads) );
        std::vector<ReadStreamStat> read_stats(file_num_);
        std::vector<size_t> current_buf_sizes(file_num_, 0);
        size_t read_count = 0;

        for (size_t i = 0; i < file_num_; ++i) {
            file_ds_[i]->seekp(0);
            read_stats[i].write(*file_ds_[i]);
        }

        size_t buf_index;
        while (!stream.eof()) {
            buf_index = read_count % file_num_;

            Read& r = buf[buf_index][current_buf_sizes[buf_index]];
            stream >> r;
            read_stats[buf_index].increase(r);

            ++current_buf_sizes[buf_index];
            VERBOSE_POWER(++read_count, " reads processed");

            if (read_count % reads_to_flush == 0) {
                for (size_t i = 0; i < file_num_; ++i) {
                    for (const Read &read : buf[i]) {
                        writer.Write(*file_ds_[i], read);
                    }
                    current_buf_sizes[i] = 0;
                }
            }
        }

        ReadStreamStat result;
        for (size_t i = 0; i < file_num_; ++i) {
            buf[i].resize(current_buf_sizes[i]);
            for (const Read &r : buf[i]) {
                writer.Write(*file_ds_[i], r);
            }

            file_ds_[i]->seekp(0);
            read_stats[i].write(*file_ds_[i]);
            result.merge(read_stats[i]);
        }

        INFO(read_count << " reads written");
        return result;
    }

public:

    BinaryWriter(const std::string& file_name_prefix, size_t file_num,
            size_t buf_size):
                file_name_prefix_(file_name_prefix),
                file_num_(file_num),
                buf_size_(buf_size) {

        for (size_t i = 0; i < file_num_; ++i) {
            std::string fname = file_name_prefix_ + "_" + std::to_string(i) + ".seq";
            file_ds_.push_back(new std::ofstream(fname, std::ios_base::binary));
        }
    }

    ~BinaryWriter() {
        for (size_t i = 0; i < file_num_; ++i) {
            if (file_ds_[i]->is_open()) {
                file_ds_[i]->close();
            }
            delete file_ds_[i];
        }
    }


    ReadStreamStat ToBinary(io::ReadStream<io::SingleReadSeq>& stream) {
        ReadBinaryWriter<io::SingleReadSeq> read_writer;
        return ToBinary(read_writer, stream, buf_size_ / file_num_);
    }

    ReadStreamStat ToBinary(io::ReadStream<io::SingleRead>& stream) {
        ReadBinaryWriter<io::SingleRead> read_writer;
        return ToBinary(read_writer, stream, buf_size_ / file_num_);
    }

    ReadStreamStat ToBinary(io::ReadStream<io::PairedReadSeq>& stream,
                            LibraryOrientation orientation = LibraryOrientation::Undefined) {
        PairedReadBinaryWriter<io::PairedReadSeq> read_writer(orientation);
        return ToBinary(read_writer, stream, buf_size_ / (2 * file_num_));
    }

    ReadStreamStat ToBinary(io::ReadStream<io::PairedRead>& stream,
                            LibraryOrientation orientation = LibraryOrientation::Undefined) {
        PairedReadBinaryWriter<io::PairedRead> read_writer(orientation);
        return ToBinary(read_writer, stream, buf_size_ / (2 * file_num_));
    }

};

}
