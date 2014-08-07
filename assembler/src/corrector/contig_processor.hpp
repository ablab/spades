/*
 * contig_processor.hpp
 *
 *  Created on: Jun 27, 2014
 *      Author: lab42
 */

#pragma once
#include "sam_reader.hpp"
#include "read.hpp"
#include "interesting_pos_processor.hpp"
#include "positional_read.hpp"

#include "io/library.hpp"
#include "openmp_wrapper.h"

#include <string>
#include <vector>
#include <unordered_map>

namespace corrector {
typedef std::vector<std::pair<string, io::LibraryType> > sam_files_type;
class ContigProcessor {
    sam_files_type sam_files_;
    std::string contig_file_;
    std::string contig_name_;
    std::string output_contig_file_;
    std::string contig_;
    bool debug_info_;
    std::vector<position_description> charts_;
    InterestingPositionProcessor ipp_;
    std::vector<int> error_counts_;

    const size_t kMaxErrorNum = 20;
public:
    // WTF: Argument name! Use brain for naming, not search & replace
    ContigProcessor(const sam_files_type &sam_files_, const std::string &contig_file)
            : sam_files_(sam_files_), contig_file_(contig_file) {
        ReadContig();
        ipp_.set_contig(contig_);
        debug_info_ = (contig_.length() > 20000);

        //Re: or just remove all info's from parallel section?
        // WTF: All the printing should be localized to outside of parallel section as much as possible
        if (debug_info_) {
#pragma omp critical
            {
                INFO("Processing contig_ " <<contig_name_ << " in thread " << omp_get_thread_num());
            }
        }
    }
    void ProcessMultipleSamFiles();
private:
    void ReadContig();
    // WTF: Is this function const?
    void UpdateOneRead(const SingleSamRead &tmp, MappedSamStream &sm);
    //returns: number of changed nucleotides;
    // WTF: Is this function const?
    size_t UpdateOneBase(size_t i, std::stringstream &ss, const std::unordered_map<size_t, position_description> &interesting_positions);

};
}
;
