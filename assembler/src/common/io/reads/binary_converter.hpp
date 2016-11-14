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

#ifndef BINARY_IO_HPP_
#define BINARY_IO_HPP_

#include <fstream>

#include "utils/verify.hpp"
#include "ireader.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"
#include "pipeline/library.hpp"

namespace io {

template<class Read>
class ReadBinaryWriter {

public:

    ReadBinaryWriter(LibraryOrientation /*orientation*/ = LibraryOrientation::Undefined) {
    }

    bool Write(std::ostream& file, const Read& r) const {
        return r.BinWrite(file);
    }
};

template<>
class ReadBinaryWriter<PairedRead> {

private:

    bool rc1_;

    bool rc2_;

public:

    ReadBinaryWriter(LibraryOrientation orientation) {
        switch (orientation) {
        case LibraryOrientation::FF:  {
            rc1_ = false;
            rc2_ = false;
            break;
        }
        case LibraryOrientation::RR:  {
            rc1_ = true;
            rc2_ = true;
            break;
        }
        case LibraryOrientation::FR:  {
            rc1_ = false;
            rc2_ = true;
            break;
        }
        case LibraryOrientation::RF:  {
            rc1_ = true;
            rc2_ = false;
            break;
        }
        default: {
            rc1_ = false;
            rc2_ = false;
            break;
        }
        }

    }

    bool Write(std::ostream& file, const PairedRead& r) const {
        return r.BinWrite(file, rc1_, rc2_);
    }
};


class BinaryWriter {

private:
    const std::string file_name_prefix_;

    size_t file_num_;

    std::vector<std::ofstream*> file_ds_;

    size_t buf_size_;

    template<class Read>
    void FlushBuffer(const std::vector<Read>& buffer, const ReadBinaryWriter<Read>& read_writer, std::ostream& file, size_t from, size_t to) {
        for (size_t i = from; i < to; ++i) {
            read_writer.Write(file, buffer[i]);
        }
    }

    template<class Read>
    void FlushBuffer(const std::vector<Read>& buffer, const ReadBinaryWriter<Read>& read_writer, std::ostream& file) {
        FlushBuffer(buffer, read_writer, file, 0, buffer.size());
    }

    template<class Read>
    ReadStreamStat ToBinary(io::ReadStream<Read>& stream, size_t buf_size,
            LibraryOrientation orientation) {

        ReadBinaryWriter<Read> read_writer(orientation);
        size_t buffer_reads = buf_size / (sizeof (Read) * 4);
        size_t reads_to_flush = buffer_reads * file_num_;

        std::vector< std::vector<Read> > buf(file_num_, std::vector<Read>(buffer_reads) );
        std::vector< ReadStreamStat > read_stats(file_num_);
        std::vector< size_t > current_buf_sizes(file_num_, 0);
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
                    FlushBuffer(buf[i], read_writer, *file_ds_[i]);
                    current_buf_sizes[i] = 0;
                }
            }
        }

        ReadStreamStat result;
        for (size_t i = 0; i < file_num_; ++i) {
            buf[i].resize(current_buf_sizes[i]);
            FlushBuffer(buf[i], read_writer, *file_ds_[i]);

            file_ds_[i]->seekp(0);
            read_stats[i].write(*file_ds_[i]);
            result.merge(read_stats[i]);
        }

        INFO(read_count << " reads written");
        return result;
    }


    template<class Read>
    ReadStreamStat ToBinaryForThread(io::ReadStream<Read>& stream, size_t buf_size,
            size_t thread_num, LibraryOrientation orientation) {

        ReadBinaryWriter<Read> read_writer(orientation);
        size_t buffer_reads = buf_size / (sizeof (Read) * 4);
        std::vector<Read> buf(buffer_reads);

        ReadStreamStat stat;
        file_ds_[thread_num]->seekp(0);
        stat.write(*file_ds_[thread_num]);

        size_t current = 0;

        while (!stream.eof()) {
            Read& r = buf[current];
            stream >> r;
            stat.increase(r);
            ++current;

            if (stat.read_count_ % buffer_reads == 0) {
                FlushBuffer(buf, read_writer, *file_ds_[thread_num]);
                current = 0;
            }
        }

        buf.resize(current);
        FlushBuffer(buf, read_writer, *file_ds_[thread_num]);

        file_ds_[thread_num]->seekp(0);
        stat.write(*file_ds_[thread_num]);

        return stat;
    }


public:

    BinaryWriter(const std::string& file_name_prefix, size_t file_num,
            size_t buf_size):
                file_name_prefix_(file_name_prefix), file_num_(file_num),
                file_ds_(), buf_size_(buf_size) {

        std::string fname;
        for (size_t i = 0; i < file_num_; ++i) {
            fname = file_name_prefix_ + "_" + ToString(i) + ".seq";
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
        return ToBinary(stream, buf_size_ / file_num_, LibraryOrientation::Undefined);
    }

    ReadStreamStat ToBinary(io::ReadStream<io::SingleRead>& stream) {
        return ToBinary(stream, buf_size_ / file_num_, LibraryOrientation::Undefined);
    }

    ReadStreamStat ToBinary(io::ReadStream<io::PairedReadSeq>& stream) {
        return ToBinary(stream, buf_size_ / (2 * file_num_), LibraryOrientation::Undefined);
    }

    ReadStreamStat ToBinary(io::ReadStream<io::PairedRead>& stream, LibraryOrientation orientation) {
        return ToBinary(stream, buf_size_ / (2 * file_num_), orientation);
    }

    ReadStreamStat ToBinaryForThread(io::ReadStream<io::SingleReadSeq>& stream, size_t thread_num) {
        return ToBinaryForThread(stream, buf_size_ / file_num_, thread_num, LibraryOrientation::Undefined);
    }

    ReadStreamStat ToBinaryForThread(io::ReadStream<io::SingleRead>& stream, size_t thread_num) {
        return ToBinaryForThread(stream, buf_size_ / file_num_, thread_num, LibraryOrientation::Undefined);
    }

    ReadStreamStat ToBinaryForThread(io::ReadStream<io::PairedReadSeq>& stream, size_t thread_num) {
        return ToBinaryForThread(stream, buf_size_ / (2 * file_num_), thread_num, LibraryOrientation::Undefined);
    }

    ReadStreamStat ToBinaryForThread(io::ReadStream<io::PairedRead>& stream, size_t thread_num, LibraryOrientation orientation) {
        return ToBinaryForThread(stream, buf_size_ / (2 * file_num_), thread_num, orientation);
    }

};


}


#endif /* BINARY_IO_HPP_ */
