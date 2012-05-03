/*
 * read_converter.hpp
 *
 *  Created on: Apr 13, 2012
 *      Author: andrey
 */

#ifndef READ_CONVERTER_HPP_
#define READ_CONVERTER_HPP_

#include "io/binary_io.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "dataset_readers.hpp"

namespace debruijn_graph {

typedef io::IReader<io::SingleReadSeq> SequenceSingleReadStream;
typedef io::IReader<io::PairedReadSeq> SequencePairedReadStream;

void convert_reads_to_binary() {
    INFO("Converting paired reads to binary format (takes a while)");
    auto_ptr<PairedReadStream> paired_reader = paired_easy_reader(false, 0);
    io::BinaryWriter paired_converter(cfg::get().paired_read_prefix, cfg::get().thread_number, cfg::get().buffer_reads);
    paired_converter.ToBinary(*paired_reader);

    INFO("Converting single reads to binary format (takes a while)");
    auto_ptr<SingleReadStream> single_reader = single_easy_reader(false, false);
    io::BinaryWriter single_converter(cfg::get().single_read_prefix, cfg::get().thread_number, cfg::get().buffer_reads);
    single_converter.ToBinary(*paired_reader);
}


std::vector<SequenceSingleReadStream*> single_binary_readers(bool followed_by_rc, bool including_paired_reads) {
    std::vector<SequenceSingleReadStream*> raw_readers(cfg::get().thread_number);

    if (including_paired_reads) {
        for (size_t i = 0; i < cfg::get().thread_number; ++i) {
            io::SeqSingleReadStream * single_stream = new io::SeqSingleReadStream(cfg::get().single_read_prefix, i);
            io::SeqPairedReadStream * paired_stream = new io::SeqPairedReadStream(cfg::get().paired_read_prefix, i, 0);
            io::SeqSingleReadStreamWrapper * single_wrapper = new io::SeqSingleReadStreamWrapper(*paired_stream);

            raw_readers[i] = new io::MultifileReader<io::SingleReadSeq>(*single_wrapper, *single_stream);
        }
    }
    else {
       for (size_t i = 0; i < cfg::get().thread_number; ++i) {
           raw_readers[i] = new io::SeqSingleReadStream(cfg::get().single_read_prefix, i);
       }
    }

    if (followed_by_rc) {
        std::vector<SequenceSingleReadStream*> rc_readers(cfg::get().thread_number);
        for (size_t i = 0; i < cfg::get().thread_number; ++i) {
            rc_readers[i] = new io::RCReaderWrapper<io::SingleReadSeq>(*raw_readers[i]);
        }
        return rc_readers;
    } else {
        return raw_readers;
    }
}

std::vector<SequencePairedReadStream*> paired_binary_readers(bool followed_by_rc, size_t insert_size) {
    std::vector<SequencePairedReadStream*> raw_readers(cfg::get().thread_number);

    for (size_t i = 0; i < cfg::get().thread_number; ++i) {
        raw_readers[i] = new io::SeqPairedReadStream(cfg::get().paired_read_prefix, i, insert_size);
    }

    if (followed_by_rc) {
        std::vector<SequencePairedReadStream*> rc_readers(cfg::get().thread_number);
        for (size_t i = 0; i < cfg::get().thread_number; ++i) {
            rc_readers[i] = new io::RCReaderWrapper<io::PairedReadSeq>(*raw_readers[i]);
        }
        return rc_readers;
    } else {
        return raw_readers;
    }
}

auto_ptr<SequenceSingleReadStream> single_binary_multireader(bool followed_by_rc, bool including_paired_reads) {
    return auto_ptr<SequenceSingleReadStream>(new io::MultifileReader<io::SingleReadSeq>(single_binary_readers(followed_by_rc, including_paired_reads)));
}

auto_ptr<SequencePairedReadStream> paired_binary_multireader(bool followed_by_rc, size_t insert_size) {
    return auto_ptr<SequencePairedReadStream>(new io::MultifileReader<io::PairedReadSeq>(paired_binary_readers(followed_by_rc, insert_size)));
}

}



#endif /* READ_CONVERTER_HPP_ */
