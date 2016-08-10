#ifndef __IO_KMER_ITERATOR_HPP__
#define __IO_KMER_ITERATOR_HPP__

#include "io/kmers/mmapped_reader.hpp"
#include <string>

namespace io {

template<class Seq>
using raw_kmer_iterator = MMappedFileRecordArrayIterator<typename Seq::DataType>;

template<class Seq>
raw_kmer_iterator<Seq> make_kmer_iterator(const std::string &FileName,
                                          unsigned K) {
    return raw_kmer_iterator<Seq>(FileName, Seq::GetDataSize(K));
}

template<class Seq>
std::vector<raw_kmer_iterator<Seq>> make_kmer_iterator(const std::string &FileName,
                                                       size_t K, size_t amount) {
    std::vector<raw_kmer_iterator<Seq>> res;
    if (amount == 1) {
        res.emplace_back(FileName, Seq::GetDataSize(K));
        return res;
    }

    // Determine the file size
    struct stat buf;
    VERIFY_MSG(stat(FileName.c_str(), &buf) != -1,
               "stat(2) failed. Reason: " << strerror(errno) << ". Error code: " << errno);
    size_t file_size = buf.st_size;

    // Now start creating the iterators keeping in mind, that offset should be
    // multiple of page size.
    size_t chunk = round_up(file_size / amount,
                            getpagesize() * Seq::GetDataSize(K) * sizeof(typename Seq::DataType));
    size_t offset = 0;
    if (chunk > file_size)
        chunk = file_size;

    while (offset < file_size) {
        res.emplace_back(FileName, Seq::GetDataSize(K),
                         offset,
                         offset + chunk > file_size ? file_size - offset : chunk);
        offset += chunk;
    }

    return res;
}


};

#endif
