//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * read_converter.hpp
 *
 *  Created on: Apr 13, 2012
 *      Author: andrey
 */

#ifndef READ_CONVERTER_HPP_
#define READ_CONVERTER_HPP_

#include <fstream>

#include "io/binary_io.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/read_stream_vector.hpp"
#include "dataset_readers.hpp"
#include "simple_tools.hpp"

namespace debruijn_graph {

typedef io::IReader<io::SingleReadSeq> SequenceSingleReadStream;
typedef io::IReader<io::PairedReadSeq> SequencePairedReadStream;

class ReadConverter {

private:
    const static size_t current_bianry_format_verstion = 5;

    void convert_reads_to_binary() {

        if (FileExists(cfg::get().temp_bin_reads_info)) {
            std::ifstream info;
            info.open(cfg::get().temp_bin_reads_info.c_str(), std::ios_base::in);

            size_t thread_num = 0;
            size_t format = 0;

            info >> format;
            if (!info.eof()) {
                info >> thread_num;
            }

            info.close();

            if (thread_num == cfg::get().max_threads && format == current_bianry_format_verstion) {
                INFO("Binary reads detected");
                return;
            }
        }

        std::ofstream info;
        info.open(cfg::get().temp_bin_reads_info.c_str(), std::ios_base::out);
        info << "0 0";
        info.close();

        INFO("Converting paired reads to binary format (takes a while)");
        auto_ptr<PairedReadStream> paired_reader = paired_easy_reader(false, 0);
        io::BinaryWriter paired_converter(cfg::get().paired_read_prefix, cfg::get().max_threads, cfg::get().buffer_size);
        io::ReadStat paired_stat = paired_converter.ToBinary(*paired_reader);

        INFO("Converting single reads to binary format (takes a while)");
        auto_ptr<SingleReadStream> single_reader = single_easy_reader(false, false);
        io::BinaryWriter single_converter(cfg::get().single_read_prefix, cfg::get().max_threads, cfg::get().buffer_size);
        io::ReadStat single_stat = single_converter.ToBinary(*single_reader);

        paired_stat.read_count_ *= 2;
        paired_stat.merge(single_stat);

        info.open(cfg::get().temp_bin_reads_info.c_str(), std::ios_base::out);

        size_t paired_libs_count = 1;
        info << current_bianry_format_verstion << " " << cfg::get().max_threads << " " << paired_libs_count << " " <<
                paired_stat.read_count_ << " " << paired_stat.max_len_ << " " << paired_stat.total_len_;
        info.close();
    }

public:
    ReadConverter() {
        convert_reads_to_binary();
    }
};


void convert_if_needed() {
    static ReadConverter converter;
}



std::vector< SequenceSingleReadStream* > apply_single_wrappers(bool followed_by_rc,
        std::vector< SequenceSingleReadStream* >& single_readers,
        std::vector< SequencePairedReadStream* > * paired_readers = 0) {

    VERIFY(single_readers.size() != 0);
    size_t size = single_readers.size();
    std::vector<SequenceSingleReadStream*> raw_readers(size);

    if (paired_readers != 0) {
        VERIFY(single_readers.size() == paired_readers->size());

        for (size_t i = 0; i < size; ++i) {
            SequenceSingleReadStream * single_stream = single_readers.at(i);
            SequencePairedReadStream * paired_stream = paired_readers->at(i);
            io::CleanSeqSingleReadStreamWrapper * single_wrapper = new io::CleanSeqSingleReadStreamWrapper(paired_stream);

            raw_readers[i] = new io::MultifileReader<io::SingleReadSeq>(*single_wrapper, *single_stream, true);
        }
    }
    else {
       for (size_t i = 0; i < size; ++i) {
           raw_readers[i] = single_readers.at(i);
       }
    }

    if (followed_by_rc) {
        std::vector<SequenceSingleReadStream*> rc_readers(size);
        for (size_t i = 0; i < size; ++i) {
            rc_readers[i] = new io::CleanRCReaderWrapper<io::SingleReadSeq>(raw_readers[i]);
        }
        return rc_readers;
    } else {
        return raw_readers;
    }
}


std::vector< SequencePairedReadStream* > apply_paired_wrappers(bool followed_by_rc,
        std::vector< SequencePairedReadStream* >& paired_readers) {

    VERIFY(paired_readers.size() != 0);
    size_t size = paired_readers.size();

    if (followed_by_rc) {
        std::vector<SequencePairedReadStream*> rc_readers(size);
        for (size_t i = 0; i < size; ++i) {
            rc_readers[i] = new io::CleanRCReaderWrapper<io::PairedReadSeq>(paired_readers[i]);
        }
        return rc_readers;
    } else {
        return paired_readers;
    }
}


std::vector< SequenceSingleReadStream* > raw_single_binary_readers(bool followed_by_rc, bool including_paired_reads) {
    convert_if_needed();

    std::vector<SequenceSingleReadStream*> single_streams(cfg::get().max_threads);
    for (size_t i = 0; i < cfg::get().max_threads; ++i) {
        single_streams[i] = new io::SeqSingleReadStream(cfg::get().single_read_prefix, i);
    }

    if (including_paired_reads) {
        std::vector<SequencePairedReadStream*> paired_streams(cfg::get().max_threads);
        for (size_t i = 0; i < cfg::get().max_threads; ++i) {
            paired_streams[i] = new io::SeqPairedReadStream(cfg::get().paired_read_prefix, i, 0);
        }
        return apply_single_wrappers(followed_by_rc, single_streams, &paired_streams);
    }
    else {
        return apply_single_wrappers(followed_by_rc, single_streams);
    }
}

std::shared_ptr<io::ReadStreamVector<SequenceSingleReadStream>> single_binary_readers(bool followed_by_rc, bool including_paired_reads) {
    return std::make_shared<io::ReadStreamVector<SequenceSingleReadStream>>(raw_single_binary_readers(followed_by_rc, including_paired_reads));
}

std::vector< SequencePairedReadStream* > raw_paired_binary_readers(bool followed_by_rc, size_t insert_size) {
    convert_if_needed();

    std::vector<SequencePairedReadStream*> paired_streams(cfg::get().max_threads);
    for (size_t i = 0; i < cfg::get().max_threads; ++i) {
        paired_streams[i] = new io::SeqPairedReadStream(cfg::get().paired_read_prefix, i, insert_size);
    }
    return apply_paired_wrappers(followed_by_rc, paired_streams);
}

std::shared_ptr<io::ReadStreamVector<SequencePairedReadStream>> paired_binary_readers(bool followed_by_rc, size_t insert_size) {
    return std::make_shared<io::ReadStreamVector<SequencePairedReadStream>>(raw_paired_binary_readers(followed_by_rc, insert_size));
}


auto_ptr<SequenceSingleReadStream> single_binary_multireader(bool followed_by_rc, bool including_paired_reads) {
    convert_if_needed();

    return auto_ptr<SequenceSingleReadStream>(new io::MultifileReader<io::SingleReadSeq>(raw_single_binary_readers(followed_by_rc, including_paired_reads), true));
}


auto_ptr<SequencePairedReadStream> paired_binary_multireader(bool followed_by_rc, size_t insert_size) {
    convert_if_needed();

    return auto_ptr<SequencePairedReadStream>(new io::MultifileReader<io::PairedReadSeq>(raw_paired_binary_readers(followed_by_rc, insert_size), true));
}

/*

class BufferedReadersStorage {

private:

    std::vector< SequenceSingleReadStream* > * single_streams_;

    std::vector< SequencePairedReadStream* > * paired_streams_;

    BufferedReadersStorage() {
        INFO("Creating buffered read storage");

        INFO("Buffering single reads... (takes a while)");
        single_streams_ = new std::vector< SequenceSingleReadStream* >(cfg::get().max_threads);
        for (size_t i = 0; i < cfg::get().max_threads; ++i) {
            io::PredictableIReader<io::SingleReadSeq> * s_stream = new io::SeqSingleReadStream(cfg::get().single_read_prefix, i);
            single_streams_->at(i) = new io::ReadBufferedStream<io::SingleReadSeq> (*s_stream);
        }

        INFO("Buffering paired reads... (takes a while)");
        paired_streams_ = new std::vector< SequencePairedReadStream* >(cfg::get().max_threads);
        for (size_t i = 0; i < cfg::get().max_threads; ++i) {
            io::PredictableIReader<io::PairedReadSeq> * p_stream = new io::SeqPairedReadStream(cfg::get().paired_read_prefix, i, 0);
            paired_streams_->at(i) = new io::ReadBufferedStream<io::PairedReadSeq> (*p_stream);
        }
    }

    BufferedReadersStorage(const BufferedReadersStorage&);

    BufferedReadersStorage& operator=(const BufferedReadersStorage&);

public:

    static BufferedReadersStorage * GetInstance() {
        static BufferedReadersStorage instance;
        return &instance;
    }


    std::vector< SequenceSingleReadStream* > * GetSingleReaders() const {
        return single_streams_;
    }

    std::vector< SequencePairedReadStream* > * GetPairedReaders() const {
        return paired_streams_;
    }

};


std::vector< SequenceSingleReadStream* > single_buffered_binary_readers(bool followed_by_rc, bool including_paired_reads) {
    convert_if_needed();

    BufferedReadersStorage * storage = BufferedReadersStorage::GetInstance();

    if (including_paired_reads) {
        return apply_single_wrappers(followed_by_rc, *(storage->GetSingleReaders()), storage->GetPairedReaders());
    }
    else {
        return apply_single_wrappers(followed_by_rc, *(storage->GetSingleReaders()));
    }
}

std::vector< SequencePairedReadStream* > paired_buffered_binary_readers(bool followed_by_rc, size_t insert_size) {
    convert_if_needed();

    BufferedReadersStorage * storage = BufferedReadersStorage::GetInstance();

    std::vector<SequencePairedReadStream*> paired_streams(cfg::get().max_threads);
    for (size_t i = 0; i < cfg::get().max_threads; ++i) {
        paired_streams[i] = new io::InsertSizeModifyingWrapper(*(storage->GetPairedReaders()->at(i)), insert_size);
    }
    return apply_paired_wrappers(followed_by_rc, paired_streams);
}

auto_ptr<SequenceSingleReadStream> single_buffered_binary_multireader(bool followed_by_rc, bool including_paired_reads) {
    convert_if_needed();

    return auto_ptr<SequenceSingleReadStream>(new io::MultifileReader<io::SingleReadSeq>(single_buffered_binary_readers(followed_by_rc, including_paired_reads)));
}

auto_ptr<SequencePairedReadStream> paired_buffered_binary_multireader(bool followed_by_rc, size_t insert_size) {
    convert_if_needed();

    return auto_ptr<SequencePairedReadStream>(new io::MultifileReader<io::PairedReadSeq>(paired_buffered_binary_readers(followed_by_rc, insert_size)));
}
*/

}

#endif /* READ_CONVERTER_HPP_ */
