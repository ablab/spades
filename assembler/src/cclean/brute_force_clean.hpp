#ifndef BRUTE_FORCE_CLEAN_HPP
#define BRUTE_FORCE_CLEAN_HPP

#include <iostream>

#include "adapter_index.hpp"
#include "io/ireadstream.hpp"
#include "sequence/sequence_tools.hpp"
#include "io/read.hpp"
#include <config_struct_cclean.hpp>
#include <ssw/ssw_cpp.h> // Striped Smith-Waterman aligner

#include "additional.cpp"

class BruteForceClean
{
  // Class that get read with oper() and clean it, if that possible
  public:
    BruteForceClean(std::ostream& aligned_output, std::ostream& output,
                    std::ostream& bed,const std::string &db,
                    const std::vector<std::string> &gen,
                    const additional::WorkModeType &brute,
                    const std::unordered_map<std::string, std::string> &options)
      : aligned_output_stream_(aligned_output), bed_stream_(bed),
        adap_seqs_(gen), brute_(brute), output_stream_(output),
        threshold_(cfg::get().mismatch_threshold), db_name_(db),
        aligned_part_fraction_(cfg::get().aligned_part_fraction),
        options_(options) { read_mlen_ = atoi(options_["mlen"].c_str()); }

    // ReadProcessor class put each read in this operator
    bool operator()(const Read &read);
    inline void BruteSimple(const Read &read);
    inline void BruteDeep(const Read &read);
    inline int aligned() { return cuted_; }

  private:
    static int cuted_;

    std::unordered_map<std::string, std::string> options_;
    const additional::WorkModeType &brute_;
    const std::vector<std::string> &adap_seqs_;
    const double aligned_part_fraction_;
    const std::string &db_name_;

    std::ostream &output_stream_;
    std::ostream &aligned_output_stream_;
    std::ostream &bed_stream_;
    unsigned threshold_;
    unsigned read_mlen_;
};

#endif // BRUTE_FORCE_CLEAN_HPP
