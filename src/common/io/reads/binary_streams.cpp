//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binary_streams.hpp"

#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"

#include <fstream>

namespace io {

bool BinaryFileSingleStream::ReadImpl(SingleReadSeq &read) {
    return read.BinRead(stream_);
}

BinaryFileSingleStream::BinaryFileSingleStream(const std::filesystem::path &file_name_prefix, size_t portion_count, size_t portion_num)
        : BinaryFileStream(file_name_prefix, portion_count, portion_num) {}

bool BinaryFilePairedStream::ReadImpl(PairedReadSeq& read) {
    return read.BinRead(stream_, insert_size_);
}

BinaryFilePairedStream::BinaryFilePairedStream(const std::filesystem::path &file_name_prefix, size_t insert_size,
                                               size_t portion_count, size_t portion_num)
        : BinaryFileStream(file_name_prefix, portion_count, portion_num), insert_size_ (insert_size) {}

PairedReadSeq BinaryUnmergingPairedStream::Convert(const SingleReadSeq &read) const {
    if (read.GetLeftOffset() >= read_length_ ||
        read.GetRightOffset() >= read_length_) {
        return PairedReadSeq();
    }
//        VERIFY(read_length_ >= read.GetLeftOffset() &&
//                       read_length_ >= read.GetRightOffset());

    size_t left_length = std::min(read.size(), read_length_ - read.GetLeftOffset());
    size_t right_length = std::min(read.size(), read_length_ - read.GetRightOffset());
    return PairedReadSeq(SingleReadSeq(read.sequence().Subseq(0, left_length), read.GetLeftOffset(), 0),
                         SingleReadSeq(read.sequence().Subseq(read.size() - right_length), 0, read.GetRightOffset()),
                         insert_size_);
}

BinaryUnmergingPairedStream::BinaryUnmergingPairedStream(const std::filesystem::path& file_name_prefix, size_t insert_size, size_t read_length,
                                                         size_t portion_count, size_t portion_num)
        : stream_(file_name_prefix, portion_count, portion_num),
          insert_size_(insert_size),
          read_length_(read_length) {}

BinaryUnmergingPairedStream& BinaryUnmergingPairedStream::operator>>(PairedReadSeq& read) {
    SingleReadSeq single_read;
    stream_ >> single_read;
    read = Convert(single_read);
    return *this;
}

}
