//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * contig_processor.hpp
 *
 *  Created on: Jun 27, 2014
 *      Author: lab42
 */

#pragma once
#include "interesting_pos_processor.hpp"
#include "positional_read.hpp"
#include "utils/parallel/openmp_wrapper.h"

#include <io/sam/sam_reader.hpp>
#include <io/sam/read.hpp>
#include "pipeline/library.hpp"

#include <string>
#include <vector>
#include <unordered_map>

namespace corrector {

using namespace sam_reader;

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
    int interesting_weight_cutoff;
protected:
    DECL_LOGGER("ContigProcessor")
public:
    ContigProcessor(const sam_files_type &sam_files, const std::string &contig_file)
            : sam_files_(sam_files), contig_file_(contig_file) {
        ReadContig();
        ipp_.set_contig(contig_);
//At least three reads to believe in inexact repeats heuristics.
        interesting_weight_cutoff = 2;
    }
    size_t ProcessMultipleSamFiles();
private:
    void ReadContig();
//Moved from read.hpp
    bool CountPositions(const SingleSamRead &read, std::unordered_map<size_t, position_description> &ps) const;
    bool CountPositions(const PairedSamRead &read, std::unordered_map<size_t, position_description> &ps) const;

    void UpdateOneRead(const SingleSamRead &tmp, MappedSamStream &sm);
    //returns: number of changed nucleotides;

    size_t UpdateOneBase(size_t i, std::stringstream &ss, const std::unordered_map<size_t, position_description> &interesting_positions) const ;

};
}
;
