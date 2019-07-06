//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binary_converter.hpp"

#include "ireader.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"
#include "orientation.hpp"
#include "pipeline/library.hpp"

#include "utils/logger/logger.hpp"
#include "utils/verify.hpp"

namespace io {

template<class Read>
class ReadBinaryWriter {
    bool rc_;

public:
    ReadBinaryWriter(bool rc = false)
            : rc_(rc) {}

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

template<class Writer, class Read>
ReadStreamStat BinaryWriter::ToBinary(const Writer &writer, io::ReadStream<Read> &stream) {
    std::vector<Read> buf;
    const size_t buf_reads = buf_size_ / (sizeof(Read) * 4);
    DEBUG("Reserving a buffer for " << buf_reads << " reads");
    buf.reserve(buf_reads);

    // Reserve space for stats
    ReadStreamStat read_stats;
    read_stats.write(*file_ds_);

    size_t rest = 1;
    auto flush_buffer = [&](){
        for (const Read &read : buf) {
            if (!--rest) {
                auto offset = (size_t)file_ds_->tellp();
                offset_ds_->write(reinterpret_cast<const char*>(&offset), sizeof(offset));
                rest = CHUNK;
            }
            writer.Write(*file_ds_, read);
        }
        buf.clear();
    };

    size_t read_count = 0;
    Read read;
    while (!stream.eof()) {
        stream >> read;
        read_stats.increase(read);
        buf.push_back(read);

        VERBOSE_POWER(++read_count, " reads processed");

        if (buf.size() == buf.capacity())
            flush_buffer();
    }
    flush_buffer(); //Write leftovers

    // Rewrite the reserved space with actual stats
    file_ds_->seekp(0);
    read_stats.write(*file_ds_);

    INFO(read_count << " reads written");
    return read_stats;
}

BinaryWriter::BinaryWriter(const std::string &file_name_prefix, size_t buf_size)
            : file_name_prefix_(file_name_prefix), buf_size_(buf_size),
              file_ds_(new std::ofstream(file_name_prefix_ + ".seq", std::ios_base::binary)),
              offset_ds_(new std::ofstream(file_name_prefix_ + ".off", std::ios_base::binary))
{}

ReadStreamStat BinaryWriter::ToBinary(io::ReadStream<io::SingleReadSeq>& stream) {
    ReadBinaryWriter<io::SingleReadSeq> read_writer;
    return ToBinary(read_writer, stream);
}

ReadStreamStat BinaryWriter::ToBinary(io::ReadStream<io::SingleRead>& stream) {
    ReadBinaryWriter<io::SingleRead> read_writer;
    return ToBinary(read_writer, stream);
}

ReadStreamStat BinaryWriter::ToBinary(io::ReadStream<io::PairedReadSeq>& stream,
                                      LibraryOrientation orientation) {
    PairedReadBinaryWriter<io::PairedReadSeq> read_writer(orientation);
    return ToBinary(read_writer, stream);
}

ReadStreamStat BinaryWriter::ToBinary(io::ReadStream<io::PairedRead>& stream,
                                      LibraryOrientation orientation) {
    PairedReadBinaryWriter<io::PairedRead> read_writer(orientation);
    return ToBinary(read_writer, stream);
}

}
