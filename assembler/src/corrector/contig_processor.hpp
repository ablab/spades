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
typedef std::vector<std::pair<std::string, io::LibraryType> > sam_files_type;
class ContigProcessor {
    sam_files_type sam_files_;
    std::string contig_file_;
    std::string contig_name_;
    std::string output_contig_file_;
    std::string contig_;
    std::vector<position_description> charts_;
    InterestingPositionProcessor ipp_;
    std::vector<int> error_counts_;

    const size_t kMaxErrorNum = 20;

public:
    ContigProcessor(const sam_files_type &sam_files, const std::string &contig_file)
            : sam_files_(sam_files), contig_file_(contig_file) {
        ReadContig();
        ipp_.set_contig(contig_);
    }
    size_t ProcessMultipleSamFiles();
private:
    void ReadContig();

    void UpdateOneRead(const SingleSamRead &tmp, MappedSamStream &sm);
    //returns: number of changed nucleotides;

    size_t UpdateOneBase(size_t i, std::stringstream &ss, const std::unordered_map<size_t, position_description> &interesting_positions) const ;

};
}
;
