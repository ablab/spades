//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include "getopt_pp/getopt_pp_standalone.h"
#include "options.hpp"
using quake_enhanced::Options;
using GetOpt::GetOpt_pp;
using GetOpt::Option;
using GetOpt::OptionPresent;
using GetOpt::Include_Environment;
using std::string;
using std::vector;
using std::endl;
using std::ostringstream;
using std::setw;

Options::Options(int argc, char **argv) :
    read_file(""),
    corrected_read_file(""),
    help_message(""),
    kmer_count_file("kmer.count"),
    hash_file_prefix("kmer_"),
    hash_file_number(1000),
    quality_offset(33),
    quality_threshold(2),
    hist_file(""),
    trusted_hist_file(""),
    bad_hist_file(""),
    top_threshold(5),
    average_min(0.9),
    limits_file(""),
    bad_threshold(0.1),
    trusted_kmer_file(""),
    bad_kmer_file("") {
  string help_module = "";
  bool need_help;
  vector<string> global_options;
  help_builder << "Usage: " << argv[0] << 
      " --read-file <file> --corrected-read-file <file> --trusted-kmer-file <file>[options]\n";
  GetOpt_pp options(argc, argv, Include_Environment);
  // Help Options
  options >> OptionPresent('h', "help", need_help);
  options >> Option('\0', "help-module", help_module);
  // General Options
  options >> Option('\0', "read-file", read_file, read_file);
  options >> Option('\0', "corrected-read-file", 
                    corrected_read_file, corrected_read_file);
  // Count Options
  options >> Option('\0', "hash-file-number",
                    hash_file_number, hash_file_number);
  options >> Option('\0', "hash-file-prefix",
                    hash_file_prefix, hash_file_prefix);
  options >> Option('\0', "quality-offset",
                    quality_offset, quality_offset);
  options >> Option('\0', "quality-threshold", 
                    quality_threshold, quality_threshold);
  options >> Option('\0', "kmer-count-file", 
                    kmer_count_file, kmer_count_file);
  // PrepareHist Options
  options >> Option('\0', "hist-file",
                    hist_file, hist_file);
  options >> Option('\0', "trusted-hist-file",
                    trusted_hist_file, trusted_hist_file);
  options >> Option('\0', "bad-hist-file",
                    bad_hist_file, bad_hist_file);
  options >> Option('\0', "top-threshold", 
                    top_threshold, top_threshold);
  options >> Option('\0', "average-min", 
                    average_min, average_min);
  // PrepareLimits Options
  options >> Option('\0', "limits-file", 
            limits_file, limits_file);
  options >> Option('\0', "bad-threshold", 
            bad_threshold, bad_threshold);
  // FilterTrusted Options
  options >> Option('\0', "trusted-kmer-file", 
            trusted_kmer_file, trusted_kmer_file);
  options >> Option('\0', "bad-kmer-file",
            bad_kmer_file, bad_kmer_file);
  if (need_help || help_module != "") {
    valid = false;
  } else {
    Validate();
  }
  help_builder << std::left << endl;
  if (!valid) {
    help_builder << 
 "General options:                                                           \n"
 "--read-file <str>                   file with reads to correct in one of   \n"
 "                                    supported formats: fastq, fasta        \n" 
 "--corrected-read-file <str>         fasta file, where corrected reads will \n"
 "                                    be written                             \n" 
 "--help-module <str>                 produce a help for a given module,     \n"
 "                                    module can be: count, prepare_hist     \n"
 "                                    prepare_limits, filter_trusted         \n";

    if (help_module == "count") {
      help_builder << 
 "Count options:                                                             \n" 
 "--kmer-count-file <str>             file where kmer count info will be     \n"
 "                                    written, default kmer.count            \n"
 "--hash-file-prefix <str>            prefix for hash_file, default: kmer_   \n"
 "--hash-file-number <int(>0)>        number of hash_files, default: 1000.   \n"
 "                                    Generally the greater this number is,  \n"
 "                                    the faster will program work, but there\n"
 "                                    is a risk of running out of file       \n"
 "                                    descriptors                            \n"
 "--quality-offset <int([0..255])>    offset of quality values (for fastq    \n"
 "                                    files). It's usually 33 or 64,         \n"
 "                                    default: 33                            \n"
 "--quality-threshold <int([0..255])> nucleotides with quality lower than    \n"
 "                                    threshold will be cut from the ends of \n"
 "                                    the read, default: 2                   \n";

    } else if (help_module == "prepare_hist") {      
      help_builder << 
 "PrepareHist options:                                                       \n" 
 "--hist-file <str>                   file where k-mer histogram will be     \n"
 "                                    written, default \"\" - no histogram   \n"
 "--trusted-hist <str>                file where trusted k-mer histogram will\n"
 "                                    be written, default \"\" - no histogram\n"
 "--bad-hist <str>                    file where bad k-mer histogram will be \n"
 "                                    written, default \"\" - no histogram   \n"
 "--top-threshold <int(>0)>           we will look for maximum which is at   \n"
 "                                    least top_threshold times higher than  \n"
 "                                    previous, default 5                    \n"
 "--average-min <float([0..1])>       trying to find Gauss's average we will \n"
 "                                    go to the left and to the right until  \n"
 "                                    we rich coverage average_min * max     \n";
    } else if (help_module == "prepare_limits") {      
      help_builder << 
 "PrepareLimits options:                                                     \n" 
 "--limits-file <str>                 file where 1-value limits for every    \n"
 "                                    k-value will be written,               \n"
 "                                    default \"\" - not to save limits      \n"
 "--bad-threshold <float(>0)>         k-mer will be considered untrusted if  \n"
 "                                    its probability of being bad is at     \n"
 "                                    least bad-threshold times greater then \n"
 "                                    probability of being good              \n";
    }  else if (help_module == "filter_trusted") {      
      help_builder << 
 "FilterTrusted options:                                                     \n" 
 "--trusted-kmer-file <str>           file where trusted k-mer will be       \n"
 "                                    written                                \n"
 "--bad--kmer-fil <str>               file where trusted k-mer will be       \n"
 "                                    written, default \"\" - no file        \n";
    }

  }
  help_message += help_builder.str();
}

void Options::Validate() {
  // General Validation
  if (read_file == "") {
    help_builder << 
        "Error: You must provide read_file\n";
    valid = false;
  }
  if (corrected_read_file == "") {
    help_builder << 
        "Error: You must provide corrected_read_file\n";
    valid = false;
  }
  // Count Validation
  if (hash_file_number < 1) {
    help_builder << 
        "Error: hash_file_number can not be lesser than one\n";
    valid = false;
  }
  if (quality_offset < 0 || quality_offset > 255) {
    help_builder << 
        "Error: quality_offset must be in 0..255\n";
    valid = false;
  }
  if (quality_threshold < 0 || quality_threshold > 255) {
    help_builder << 
        "Error: quality_threshold must be in 0..255\n";
    valid = false;
  }
  // PrepareHist Validation
  if (average_min < 0 || average_min > 1) {
    help_builder << 
        "Error: average_min must be in 0..1\n";
    valid = false;
  }
  // PrepareLimits Validation
  if (bad_threshold < 0) {
    help_builder << 
      "Error: bad_threshold must be in 0..*\n";
    valid = false;
  }
  // FilterTrusted Validation
  if (trusted_kmer_file == "") {
    help_builder << "Error: trusted_kmer_file must be provided\n";
    valid = false;
  }
}
