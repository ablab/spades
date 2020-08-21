//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_stream.hpp"
#include "single_read.hpp"
#include "paired_read.hpp"
#include "binary_converter.hpp"

#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/filesystem/file_opener.hpp"

#include <fstream>

namespace io {

template<typename SeqT>
class BinaryFileStream {
protected:
    std::ifstream stream_;

    virtual bool ReadImpl(SeqT &read) = 0;

private:
    size_t offset_, count_, current_;

    void Init() {
        stream_.clear();
        stream_.seekg(offset_);
        VERIFY_MSG(stream_.good(), "Stream is not good(), offset_ " << offset_ << " count_ " << count_);
        current_ = 0;
    }

public:
    /**
     * @brief Constructs a reader of a portion of reads.
     * @param file_name_prefix
     * @param portion_count Total number of (roughly equal) portions.
     * @param portion_num Index of the portion (0..portion_count - 1).
     */
    BinaryFileStream(const std::string &file_name_prefix, size_t portion_count, size_t portion_num) {
        DEBUG("Preparing binary stream #" << portion_num << "/" << portion_count);
        VERIFY(portion_num < portion_count);
        const std::string fname = file_name_prefix + ".seq";
        stream_.open(fname, std::ios_base::binary | std::ios_base::in);
        ReadStreamStat stat;
        stat.read(stream_);

        const std::string offset_name = file_name_prefix + ".off";
        const size_t chunk_count = fs::filesize(offset_name) / sizeof(size_t);

        // We split all read chunks into portion_count portions
        // Portion could have size (chunk_count / portion_count) or (chunk_count / portion_count + 1)
        const size_t small_portion_size = chunk_count / portion_count;
        const size_t big_portion_size = small_portion_size + 1;
        const size_t big_portion_count = chunk_count % portion_count;
        VERIFY(big_portion_count * big_portion_size + (portion_count - big_portion_count) * small_portion_size == chunk_count);
        // We suppose that all small portions are placed at the end
        // Note that small portion could have size 0

        // Compute the number of big portions preceding the current portion
        const size_t big_portion_before = std::min(portion_num, big_portion_count);

        // At last, compute the number of the first chunk in the current portion
        const size_t chunk_num = big_portion_before * (big_portion_size - small_portion_size) + portion_num * small_portion_size;
        VERIFY_MSG(chunk_num <= chunk_count, "chunk_num " << chunk_num << " chunk_count " << chunk_count << " big_portion_before " << big_portion_before << " big_portion_size " << big_portion_size << " small_portion_size " << small_portion_size << " portion_num " << portion_num);

        if (chunk_num < chunk_count) {  // if we start from existing chunk
            // Calculating the absolute offset in the reads file
            auto offset_stream = fs::open_file(offset_name, std::ios_base::binary | std::ios_base::in);
            offset_stream.seekg(chunk_num * sizeof(size_t));
            offset_stream.read(reinterpret_cast<char *>(&offset_), sizeof(offset_));
            DEBUG("Offset read: " << offset_ << " chunk_count " << chunk_count << " chunk_num " << chunk_num << " portion_count " << portion_count << " portion_num " << portion_num << " prefix " << file_name_prefix << " name " << offset_name);
            VERIFY(offset_stream);
            const bool is_big_portion = portion_num < big_portion_count;
            const size_t start_num = chunk_num * BinaryWriter::CHUNK;
            // Last chunk could be incomplete => we should truncate count_ for last portions
            count_ = std::min(stat.read_count - start_num,
                              (is_big_portion ? big_portion_size : small_portion_size) * BinaryWriter::CHUNK);

            DEBUG("Reads " << start_num << "-" << start_num + count_ << "/" << stat.read_count << " from " << offset_);
        } else {  // current portion has size 0 (the case of chunk_count == 0 is also included here)
            // Setup safe offset value
            offset_ = sizeof(ReadStreamStat);
            count_ = 0;
            DEBUG("Empty BinaryFileStream constructed");
        }

        Init();
    }

    /**
     * @brief Constructs a reader of all available reads.
     */
    BinaryFileStream(const std::string &file_name_prefix)
            : BinaryFileStream(file_name_prefix, 1, 0) {}

    BinaryFileStream<SeqT>& operator>>(SeqT &read) {
        ReadImpl(read);
        VERIFY(current_ < count_);
        ++current_;
        return *this;
    }

    bool is_open() {
        return stream_.is_open();
    }

    bool eof() {
        return current_ >= count_;
    }

    void close() {
        current_ = 0;
        stream_.close();
    }

    void reset() {
        Init();
    }

};

class BinaryFileSingleStream : public BinaryFileStream<SingleReadSeq>  {
protected:
    bool ReadImpl(SingleReadSeq &read) override;
public:
    BinaryFileSingleStream(const std::string &file_name_prefix, size_t portion_count, size_t portion_num);
};

class BinaryFilePairedStream: public BinaryFileStream<PairedReadSeq> {
    size_t insert_size_;
protected:
    bool ReadImpl(PairedReadSeq& read) override;
public:
    BinaryFilePairedStream(const std::string &file_name_prefix, size_t insert_size,
                           size_t portion_count, size_t portion_num);
};

// returns FF oriented paired reads
class BinaryUnmergingPairedStream {
    BinaryFileSingleStream stream_;
    size_t insert_size_;
    size_t read_length_;

    PairedReadSeq Convert(const SingleReadSeq &read) const;

public:
    BinaryUnmergingPairedStream(const std::string& file_name_prefix, size_t insert_size, size_t read_length,
                                size_t portion_count, size_t portion_num);

    bool is_open() { return stream_.is_open(); }
    bool eof() { return stream_.eof(); }
    void close() { stream_.close(); }
    void reset() { stream_.reset(); }

    BinaryUnmergingPairedStream& operator>>(PairedReadSeq& read);

};

}
