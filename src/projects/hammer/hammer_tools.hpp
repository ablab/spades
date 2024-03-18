//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef HAMMER_TOOLS_HPP
#define HAMMER_TOOLS_HPP

#include "globals.hpp"
#include "kmer_stat.hpp"

#include "io/kmers/mmapped_reader.hpp"
#include "io/reads/ireadstream.hpp"
#include "io/reads/read.hpp"
#include "sequence/seq.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <unordered_map>

namespace hammer {

/// initialize subkmer positions and log about it
void InitializeSubKMerPositions(int tau);

struct CorrectionStats {
  size_t changedReads;
  size_t changedNucleotides;
  size_t uncorrectedNucleotides;
  size_t totalNucleotides;
  CorrectionStats() : changedReads(0),
                      changedNucleotides(0),
                      uncorrectedNucleotides(0),
                      totalNucleotides(0) {}

  CorrectionStats& operator +=(const CorrectionStats &rhs) {
    changedReads += rhs.changedReads;
    changedNucleotides += rhs.changedNucleotides;
    uncorrectedNucleotides += rhs.uncorrectedNucleotides;
    totalNucleotides += rhs.totalNucleotides;
    return *this;
  }
};

/// parallel correction of batch of reads
CorrectionStats CorrectReadsBatch(std::vector<bool> &res, std::vector<Read> &reads, size_t buf_size,
                       size_t &changedReads, size_t &changedNucleotides, size_t &uncorrectedNucleotides, size_t &totalNucleotides,
                       const KMerData &data);

/// correct reads in a given file
CorrectionStats CorrectReadFile(const KMerData &data,
                         size_t &changedReads, size_t &changedNucleotides, size_t &uncorrectedNucleotides, size_t &totalNucleotides,
                         const std::filesystem::path &fname,
                         std::ofstream *outf_good, std::ofstream *outf_bad);

/// correct reads in a given pair of files
CorrectionStats CorrectPairedReadFiles(const KMerData &data,
                            size_t &changedReads, size_t &changedNucleotides, size_t &uncorrectedNucleotides, size_t &totalNucleotides,
                            const std::filesystem::path &fnamel, const std::string &fnamer,
                            std::ofstream * ofbadl, std::ofstream * ofcorl, std::ofstream * ofbadr, std::ofstream * ofcorr, std::ofstream * ofunp);
/// correct all reads
size_t CorrectAllReads();

std::filesystem::path getFilename(const std::filesystem::path & dirprefix, const std::string & suffix );
std::filesystem::path getFilename(const std::filesystem::path & dirprefix, unsigned iter_count, const std::string & suffix );
std::filesystem::path getFilename(const std::filesystem::path & dirprefix, int iter_count, const std::string & suffix, int suffix_num );
std::filesystem::path getFilename(const std::filesystem::path & dirprefix, int iter_count, const std::string & suffix, int suffix_num, const std::string & suffix2 );
std::filesystem::path getFilename(const std::filesystem::path & dirprefix, const std::string & suffix, int suffix_num );
std::filesystem::path getReadsFilename(const std::filesystem::path & dirprefix, const std::string &fname, unsigned iter_no, const std::string & suffix);
};



#endif
