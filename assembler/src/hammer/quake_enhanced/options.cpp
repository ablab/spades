#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include "libs/getopt_pp/getopt_pp_standalone.h"
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
      " --read-file file --corrected-read-file file [options]" << endl;
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
  if (need_help || help_module != "") {
    valid = false;
  } else {
    Validate();
  }
  help_builder << std::left << endl;
  const int kTab = 35;
  if (!valid) {
    help_builder << "General options:"
                 <<  endl;
    help_builder << setw(kTab) 
                 << "--read-file str" 
                 << "file with reads to correct in one of supported formats:" 
                 << endl
                 << setw(kTab + 2) 
                 << " " 
                 << "fastq, fasta" 
                 <<  endl;
    help_builder << setw(kTab) 
                 << "--corrected-read-file str" 
                 << "fasta file, where corrected reads will be written" 
                 << endl;
    help_builder << setw(kTab) 
                 << "--help-module str" 
                 << "produce a help for a given module, module can be:" 
                 << endl 
                 << setw(kTab + 2) 
                 << " " 
                 << "count" 
                 << endl;
    if (help_module == "count") {
      help_builder << "Count options:" 
                   <<  endl;
      help_builder << setw(kTab)
                   << "--hash-file-prefix str"
                   << "prefix for hash_file, default: kmer_"
                   << endl;
      help_builder << setw(kTab)
                   << "--hash-file-number int(>0)"
                   << "number of hash_files, default: 1000. Generally the greater"
                   << endl
                   << setw(kTab + 2)
                   << ""
                   << "this number is, the faster will program work, but there is"
                   << endl
                   << setw(kTab + 2)
                   << ""
                   << "a risk of running out of file descriptors"
                   << endl;
      help_builder << setw(kTab)
                   << "--quality-offset int([0..255])"
                   << "offset of quality values (for fastq files)."
                   << endl
                   << setw(kTab + 2)
                   << ""
                   << "It's usually 33 or 64, default: 33"
                   << endl;
      help_builder << setw(kTab)
                   << "--quality-threshold int([0..255])"
                   << "nucleotides with quality lower than threshold will be cut"
                   << endl 
                   << setw(kTab + 2) 
                   << ""
                   << "from the ends of the read, default: 2"
                   << endl;      
    } else if (help_module == "prepare_hist") {      
    }
  }
  help_message += help_builder.str();
}

void Options::Validate() {
  // General Validation
  if (read_file == "") {
    help_builder << 
        "Error: You must provide read_file" << 
        std::endl;
    valid = false;
  }
  if (corrected_read_file == "") {
    help_builder << 
        "Error: You must provide corrected_read_file" << 
        std::endl;
    valid = false;
  }
  // Count Validation
  if (hash_file_number < 1) {
    help_builder << 
        "Error: hash_file_number can not be lesser than one" << 
        std::endl;
    valid = false;
  }
  if (quality_offset < 0 || quality_offset > 255) {
    help_builder << 
        "Error: quality_offset must be in 0..255" << 
        std::endl;
    valid = false;
  }
  if (quality_threshold < 0 || quality_threshold > 255) {
    help_builder << 
        "Error: quality_threshold must be in 0..255" << 
        std::endl;
    valid = false;
  }

}
