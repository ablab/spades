//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_GLOBALS_HPP_
#define HAMMER_GLOBALS_HPP_

#include "kmer_stat.hpp"

class KMerData;

struct Globals {
  static int iteration_no;
  static std::vector<std::string> input_filenames;
  static std::vector<std::string> input_filename_bases;
  static std::vector<hint_t> input_file_blob_positions;
  static std::vector<size_t> input_file_sizes;

  static std::vector<uint32_t> * subKMerPositions;
  static char* blob;
  static char* blobquality;
  static std::vector<PositionRead> * pr;
  static KMerData *kmer_data;
  static hint_t blob_max_size;
  static hint_t blob_size;
  static hint_t revNo;

  static char char_offset;
  static bool char_offset_user;

  static bool use_common_quality;
  static char common_quality;
  static double common_kmer_errprob;
  static double quality_probs[256];
  static double quality_lprobs[256];

  static void writeBlob( const char * fname );
  static void readBlob( const char * fname );
};

#endif //  HAMMER_GLOBALS_HPP_

