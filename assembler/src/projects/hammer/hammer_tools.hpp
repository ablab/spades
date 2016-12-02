//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef HAMMER_TOOLS_HPP
#define HAMMER_TOOLS_HPP

#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include "io/reads/read.hpp"
#include "io/reads/ireadstream.hpp"
#include "sequence/seq.hpp"
#include "globals.hpp"
#include "kmer_stat.hpp"
#include "io/kmers/mmapped_reader.hpp"

namespace hammer {

/// initialize subkmer positions and log about it
void InitializeSubKMerPositions();

/// parallel correction of batch of reads
void CorrectReadsBatch(std::vector<bool> &res, std::vector<Read> &reads, size_t buf_size,
                       size_t &changedReads, size_t &changedNucleotides, size_t &uncorrectedNucleotides, size_t &totalNucleotides,
                       const KMerData &data);

/// correct reads in a given file
void CorrectReadFile(const KMerData &data,
                     size_t &changedReads, size_t &changedNucleotides, size_t &uncorrectedNucleotides, size_t &totalNucleotides,
                     const std::string &fname,
                     std::ofstream *outf_good, std::ofstream *outf_bad);

/// correct reads in a given pair of files
void CorrectPairedReadFiles(const KMerData &data,
                            size_t &changedReads, size_t &changedNucleotides, size_t &uncorrectedNucleotides, size_t &totalNucleotides,
                            const std::string &fnamel, const std::string &fnamer,
                            std::ofstream * ofbadl, std::ofstream * ofcorl, std::ofstream * ofbadr, std::ofstream * ofcorr, std::ofstream * ofunp);
/// correct all reads
size_t CorrectAllReads();

std::string getFilename(const std::string & dirprefix, const std::string & suffix );
std::string getFilename(const std::string & dirprefix, unsigned iter_count, const std::string & suffix );
std::string getFilename(const std::string & dirprefix, int iter_count, const std::string & suffix, int suffix_num );
std::string getFilename(const std::string & dirprefix, int iter_count, const std::string & suffix, int suffix_num, const std::string & suffix2 );
std::string getFilename(const std::string & dirprefix, const std::string & suffix, int suffix_num );
std::string getReadsFilename(const std::string & dirprefix, const std::string &fname, unsigned iter_no, const std::string & suffix);
};



#endif
