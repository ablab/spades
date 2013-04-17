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

#include "io/binary_converter.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/read_stream_vector.hpp"
#include "dataset_readers.hpp"
#include "simple_tools.hpp"

namespace debruijn_graph {

typedef io::IReader<io::SingleReadSeq> SequenceSingleReadStream;
typedef io::IReader<io::PairedReadSeq> SequencePairedReadStream;

class ReadConverter {

private:
    const static size_t current_bianry_format_verstion = 6;

    void convert_reads_to_binary() {

        if (FileExists(cfg::get().temp_bin_reads_info)) {
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

            info.close();

            if (thread_num == cfg::get().max_threads && format == current_bianry_format_verstion && lib_count == cfg::get().ds.reads.lib_count()) {
                INFO("Binary reads detected");

                auto &dataset = cfg::get_writable().ds.reads;
                for (size_t i = 0; i < dataset.lib_count(); ++i) {
                    dataset[i].data().thread_num = cfg::get().max_threads;
                    dataset[i].data().paired_read_prefix = cfg::get().paired_read_prefix + "_" + ToString(i);
                    dataset[i].data().single_read_prefix = cfg::get().single_read_prefix + "_" + ToString(i);
                }
                return;
            }
        }

        std::ofstream info;
        info.open(cfg::get().temp_bin_reads_info.c_str(), std::ios_base::out);
        info << "0 0";
        info.close();

        io::ReadStat total_stat;
        auto &dataset = cfg::get_writable().ds.reads;

        INFO("Converting reads to binary format (takes a while)");
        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            INFO("Paired reads for library #" << i);
            dataset[i].data().thread_num = cfg::get().max_threads;
            dataset[i].data().paired_read_prefix = cfg::get().paired_read_prefix + "_" + ToString(i);

            auto_ptr<PairedReadStream> paired_reader = paired_easy_reader(dataset[i], false, 0);
            io::BinaryWriter paired_converter(dataset[i].data().paired_read_prefix, cfg::get().max_threads, cfg::get().buffer_size);
            io::ReadStat paired_stat = paired_converter.ToBinary(*paired_reader);
            paired_stat.read_count_ *= 2;
            total_stat.merge(paired_stat);

            INFO("Single reads for library #" << i);
            dataset[i].data().single_read_prefix = cfg::get().single_read_prefix + "_" + ToString(i);

            auto_ptr<SingleReadStream> single_reader = single_easy_reader(dataset[i], false, false);
            io::BinaryWriter single_converter(dataset[i].data().single_read_prefix, cfg::get().max_threads, cfg::get().buffer_size);
            io::ReadStat single_stat = single_converter.ToBinary(*single_reader);
            total_stat.merge(single_stat);
        }

        info.open(cfg::get().temp_bin_reads_info.c_str(), std::ios_base::out);
        info << current_bianry_format_verstion << " " << cfg::get().max_threads << " " << cfg::get().ds.reads.lib_count() << " " <<
                total_stat.read_count_ << " " << total_stat.max_len_ << " " << total_stat.total_len_;
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

std::vector< SequencePairedReadStream* > raw_paired_binary_readers(const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                                                   bool followed_by_rc,
                                                                   size_t insert_size = 0) {
    convert_if_needed();

    std::vector<SequencePairedReadStream*> paired_streams(lib.data().thread_num);
    for (size_t i = 0; i < lib.data().thread_num; ++i) {
        paired_streams[i] = new io::SeqPairedReadStream(lib.data().paired_read_prefix, i, insert_size);
    }
    return io::apply_paired_wrappers(followed_by_rc, paired_streams);
}

std::vector< SequenceSingleReadStream* > raw_single_binary_readers(const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                                                   bool followed_by_rc,
                                                                   bool including_paired_reads) {
    convert_if_needed();

    std::vector<SequenceSingleReadStream*> single_streams(lib.data().thread_num);
    for (size_t i = 0; i < lib.data().thread_num; ++i) {
        single_streams[i] = new io::SeqSingleReadStream(lib.data().single_read_prefix, i);
    }
    if (including_paired_reads) {
        std::vector<SequencePairedReadStream*> paired_streams(lib.data().thread_num);
        for (size_t i = 0; i < lib.data().thread_num; ++i) {
            paired_streams[i] = new io::SeqPairedReadStream(lib.data().paired_read_prefix, i, 0);
        }

        return io::apply_single_wrappers(followed_by_rc, single_streams, &paired_streams);
    }
    else {
        return io::apply_single_wrappers(followed_by_rc, single_streams);
    }
}


io::ReadStreamVector< SequencePairedReadStream > paired_binary_readers(const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                                                       bool followed_by_rc,
                                                                       size_t insert_size = 0) {
    convert_if_needed();
    return io::ReadStreamVector< SequencePairedReadStream >(raw_paired_binary_readers(lib, followed_by_rc, insert_size));
}


io::ReadStreamVector< SequenceSingleReadStream > single_binary_readers(const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                                                       bool followed_by_rc,
                                                                       bool including_paired_reads) {
    convert_if_needed();
    return io::ReadStreamVector< SequenceSingleReadStream >(raw_single_binary_readers(lib, followed_by_rc, including_paired_reads));
}


io::ReadStreamVector< SequencePairedReadStream > paired_binary_readers_for_libs(const std::vector<size_t>& libs,
                                                   bool followed_by_rc,
                                                   size_t insert_size = 0) {
    convert_if_needed();

    std::vector< std::vector< SequencePairedReadStream* > > streams(cfg::get().max_threads);
    for (size_t i = 0; i < libs.size(); ++i) {
      std::vector< SequencePairedReadStream* > lib_streams = raw_paired_binary_readers(cfg::get().ds.reads[libs[i]], followed_by_rc, insert_size);

      for (size_t j = 0; j < cfg::get().max_threads; ++j) {
          streams[j].push_back(lib_streams[j]);
      }
    }

    std::vector< SequencePairedReadStream* > joint_streams;
    for (size_t j = 0; j < cfg::get().max_threads; ++j) {
      joint_streams.push_back(new io::MultifileReader<io::PairedReadSeq>(streams[j], true));
    }
    return io::ReadStreamVector< SequencePairedReadStream >(joint_streams);
}


io::ReadStreamVector< SequenceSingleReadStream > single_binary_readers_for_libs(const std::vector<size_t>& libs,
                                                   bool followed_by_rc,
                                                   bool including_paired_reads) {
    convert_if_needed();

    std::vector< std::vector< SequenceSingleReadStream* > > streams(cfg::get().max_threads);
    for (size_t i = 0; i < libs.size(); ++i) {
      std::vector< SequenceSingleReadStream* > lib_streams = raw_single_binary_readers(cfg::get().ds.reads[libs[i]], followed_by_rc, including_paired_reads);

      for (size_t j = 0; j < cfg::get().max_threads; ++j) {
          streams[j].push_back(lib_streams[j]);
      }
    }

    std::vector< SequenceSingleReadStream* > joint_streams;
    for (size_t j = 0; j < cfg::get().max_threads; ++j) {
      joint_streams.push_back(new io::MultifileReader<io::SingleReadSeq>(streams[j], true));
    }
    return io::ReadStreamVector< SequenceSingleReadStream >(joint_streams);
}


io::ReadStreamVector< SequencePairedReadStream > paired_binary_readers(bool followed_by_rc,
                                                   size_t insert_size = 0) {
  std::vector<size_t> all_libs(cfg::get().ds.reads.lib_count());
  for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
      all_libs[i] = i;
  }
  return paired_binary_readers_for_libs(all_libs, followed_by_rc, insert_size);
}


io::ReadStreamVector< SequenceSingleReadStream > single_binary_readers(bool followed_by_rc,
                                                   bool including_paired_reads) {
  std::vector<size_t> all_libs(cfg::get().ds.reads.lib_count());
  for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
    all_libs[i] = i;
  }
  return single_binary_readers_for_libs(all_libs, followed_by_rc, including_paired_reads);
}


auto_ptr<SequenceSingleReadStream> single_binary_multireader(bool followed_by_rc, bool including_paired_reads) {
    auto readers = single_binary_readers(followed_by_rc, including_paired_reads);
    readers.release();
    return auto_ptr<SequenceSingleReadStream>(new io::MultifileReader<io::SingleReadSeq>(readers.get(), true));
}


auto_ptr<SequencePairedReadStream> paired_binary_multireader(bool followed_by_rc, size_t insert_size = 0) {
    auto readers = paired_binary_readers(followed_by_rc, insert_size);
    readers.release();
    return auto_ptr<SequencePairedReadStream>(new io::MultifileReader<io::PairedReadSeq>(readers.get(), true));
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
