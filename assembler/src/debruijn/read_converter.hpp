//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * read_converter.hpp
 *
 *  Created on: Apr 13, 2012
 *      Author: andrey
 */

#pragma once

#include "io/binary_converter.hpp"
#include "io/io_helper.hpp"
#include "dataset_readers.hpp"
#include "simple_tools.hpp"

#include <fstream>

namespace debruijn_graph {

class ReadConverter {

private:
    const static size_t current_binary_format_version = 8;

    void convert_reads_to_binary() {

        if (path::FileExists(cfg::get().temp_bin_reads_info)) {
            std::ifstream info;
            info.open(cfg::get().temp_bin_reads_info.c_str(), std::ios_base::in);

            size_t thread_num = 0;
            size_t format = 0;
            size_t lib_count = 0;

            info >> format;
            if (!info.eof()) {
                info >> thread_num;
            }
            if (!info.eof()) {
                info >> lib_count;
            }

            if (thread_num == cfg::get().max_threads && format == current_binary_format_version  && lib_count == cfg::get().ds.reads.lib_count()) {
                INFO("Binary reads detected");

                io::ReadStreamStat stat;
                info >> stat.read_count_;
                info >> stat.max_len_;
                info >> stat.total_len_;

                auto &dataset = cfg::get_writable().ds.reads;
                for (size_t i = 0; i < dataset.lib_count(); ++i) {
                    info >> dataset[i].data().read_length;
                    info >> dataset[i].data().total_nucls;

                    dataset[i].data().thread_num = cfg::get().max_threads;
                    dataset[i].data().paired_read_prefix = cfg::get().paired_read_prefix + "_" + ToString(i);
                    dataset[i].data().single_read_prefix = cfg::get().single_read_prefix + "_" + ToString(i);
                }
                info.close();
                return;
            }
            info.close();
        }

        std::ofstream info;
        info.open(cfg::get().temp_bin_reads_info.c_str(), std::ios_base::out);
        info << "0 0";
        info.close();

        io::ReadStreamStat total_stat;
        auto& dataset = cfg::get_writable().ds.reads;

        INFO("Converting reads to binary format (takes a while)");
        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            INFO("Paired reads for library #" << i);
            dataset[i].data().thread_num = cfg::get().max_threads;
            dataset[i].data().paired_read_prefix = cfg::get().paired_read_prefix + "_" + ToString(i);

            io::PairedStreamPtr paired_reader = paired_easy_reader(dataset[i], false, 0, false, false);
            io::BinaryWriter paired_converter(dataset[i].data().paired_read_prefix, cfg::get().max_threads, cfg::get().buffer_size);
            io::ReadStreamStat paired_stat = paired_converter.ToBinary(*paired_reader, dataset[i].orientation());
            paired_stat.read_count_ *= 2;
            total_stat.merge(paired_stat);

            INFO("Single reads for library #" << i);
            dataset[i].data().single_read_prefix = cfg::get().single_read_prefix + "_" + ToString(i);
            io::SingleStreamPtr single_reader = single_easy_reader(dataset[i], false, false);
            io::BinaryWriter single_converter(dataset[i].data().single_read_prefix, cfg::get().max_threads, cfg::get().buffer_size);
            io::ReadStreamStat single_stat = single_converter.ToBinary(*single_reader);
            total_stat.merge(single_stat);

            paired_stat.merge(single_stat);
            dataset[i].data().read_length = paired_stat.max_len_;
            dataset[i].data().total_nucls = paired_stat.total_len_;
        }
        info.open(cfg::get().temp_bin_reads_info.c_str(), std::ios_base::out);
        info << current_binary_format_version << " " << cfg::get().max_threads << " " << cfg::get().ds.reads.lib_count() << " " <<
                total_stat.read_count_ << " " << total_stat.max_len_ << " " << total_stat.total_len_ << "\n";

        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            info << dataset[i].data().read_length << " " << dataset[i].data().total_nucls << "\n";
        }
        info.close();
    }

public:
    ReadConverter() {
        convert_reads_to_binary();
    }
};


inline
void convert_if_needed() {
    static ReadConverter converter;
}

inline
io::BinaryPairedStreams raw_paired_binary_readers(const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                                                   bool followed_by_rc,
                                                                   size_t insert_size = 0) {
    convert_if_needed();

    io::ReadStreamList<io::PairedReadSeq> paired_streams;
    for (size_t i = 0; i < lib.data().thread_num; ++i) {
        paired_streams.push_back(make_shared<io::BinaryFilePairedStream>(lib.data().paired_read_prefix, i, insert_size));
    }
    return io::apply_paired_wrappers(followed_by_rc, paired_streams);
}

inline
io::BinarySingleStreams raw_single_binary_readers(const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                                                   bool followed_by_rc,
                                                                   bool including_paired_reads) {
    convert_if_needed();

    io::BinarySingleStreams single_streams;
    for (size_t i = 0; i < lib.data().thread_num; ++i) {
        single_streams.push_back(make_shared<io::BinaryFileSingleStream>(lib.data().single_read_prefix, i));
    }
    if (including_paired_reads) {
        io::BinaryPairedStreams paired_streams;
        for (size_t i = 0; i < lib.data().thread_num; ++i) {
            paired_streams.push_back(make_shared<io::BinaryFilePairedStream>(lib.data().paired_read_prefix, i, 0));
        }

        return io::apply_single_wrappers(followed_by_rc, single_streams, &paired_streams);
    }
    else {
        return io::apply_single_wrappers(followed_by_rc, single_streams);
    }
}


inline
io::BinaryPairedStreams paired_binary_readers(const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                                                       bool followed_by_rc,
                                                                       size_t insert_size = 0) {
    convert_if_needed();
    return raw_paired_binary_readers(lib, followed_by_rc, insert_size);
}


inline
io::BinarySingleStreams single_binary_readers(const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                                                       bool followed_by_rc,
                                                                       bool including_paired_reads) {
    convert_if_needed();
    return raw_single_binary_readers(lib, followed_by_rc, including_paired_reads);
}


inline
//todo simplify
io::BinaryPairedStreams paired_binary_readers_for_libs(const std::vector<size_t>& libs,
                                                   bool followed_by_rc,
                                                   size_t insert_size = 0) {
    convert_if_needed();

    std::vector<io::BinaryPairedStreams> streams(cfg::get().max_threads);
    for (size_t i = 0; i < libs.size(); ++i) {
      io::BinaryPairedStreams lib_streams = raw_paired_binary_readers(cfg::get().ds.reads[libs[i]], followed_by_rc, insert_size);

      for (size_t j = 0; j < cfg::get().max_threads; ++j) {
          streams[j].push_back(lib_streams.ptr_at(j));
      }
    }

    io::BinaryPairedStreams joint_streams;
    for (size_t j = 0; j < cfg::get().max_threads; ++j) {
      joint_streams.push_back(io::MultifileWrap<io::PairedReadSeq>(streams[j]));
    }
    return joint_streams;
}

inline
io::BinarySingleStreams single_binary_readers_for_libs(const std::vector<size_t>& libs,
                                                   bool followed_by_rc,
                                                   bool including_paired_reads) {
    convert_if_needed();

    std::vector<io::BinarySingleStreams> streams(cfg::get().max_threads);
    for (size_t i = 0; i < libs.size(); ++i) {
      io::BinarySingleStreams lib_streams = raw_single_binary_readers(cfg::get().ds.reads[libs[i]], followed_by_rc, including_paired_reads);

      for (size_t j = 0; j < cfg::get().max_threads; ++j) {
          streams[j].push_back(lib_streams.ptr_at(j));
      }
    }

    io::BinarySingleStreams joint_streams;
    for (size_t j = 0; j < cfg::get().max_threads; ++j) {
      joint_streams.push_back(io::MultifileWrap<io::SingleReadSeq>(streams[j]));
    }
    return joint_streams;
}

inline
io::BinaryPairedStreams  paired_binary_readers(bool followed_by_rc,
                                                   size_t insert_size = 0) {
  std::vector<size_t> all_libs(cfg::get().ds.reads.lib_count());
  for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
      all_libs[i] = i;
  }
  return paired_binary_readers_for_libs(all_libs, followed_by_rc, insert_size);
}

inline
io::BinarySingleStreams single_binary_readers(bool followed_by_rc,
                                                   bool including_paired_reads) {
  std::vector<size_t> all_libs(cfg::get().ds.reads.lib_count());
  for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
    all_libs[i] = i;
  }
  return single_binary_readers_for_libs(all_libs, followed_by_rc, including_paired_reads);
}

inline
io::BinarySingleStreamPtr single_binary_multireader(bool followed_by_rc, bool including_paired_reads) {
    return io::MultifileWrap<io::SingleReadSeq>(single_binary_readers(followed_by_rc, including_paired_reads));
}

inline
io::BinaryPairedStreamPtr paired_binary_multireader(bool followed_by_rc, size_t insert_size = 0) {
    return io::MultifileWrap<io::PairedReadSeq>(paired_binary_readers(followed_by_rc, insert_size));
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
