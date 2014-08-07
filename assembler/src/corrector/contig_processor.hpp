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
    // WTF: Coding style! Name member vars properly!
    sam_files_type sam_files;
    std::string sam_file;
    std::string contig_file;
    std::string contig_name;
    std::string output_contig_file;
    std::string contig;
    bool debug_info;
    std::vector<position_description> charts;
    InterestingPositionProcessor ipp;
    std::vector<int> error_counts;
    const size_t max_error_num = 20;
public:
    ContigProcessor(const sam_files_type &sam_files_, const std::string &contig_file)
            : sam_files(sam_files_), contig_file(contig_file) {
        ReadContig();
        ipp.set_contig(contig);
        debug_info = (contig.length() > 20000);
        // WTF: race on INFO()
        // Re: atomic INFO's or no INFO's in parallel code?
        if (debug_info) {
            INFO("Processing contig " <<contig_name << " in thread " << omp_get_thread_num());
        }
        DEBUG("Processing contig " <<contig_name << " in thread " << omp_get_thread_num());

    }
    void ProcessMultipleSamFiles();
private:
    void ReadContig();
    void UpdateOneRead(const SingleSamRead &tmp, MappedSamStream &sm);
    //returns: number of changed nucleotides;
    size_t UpdateOneBase(size_t i, std::stringstream &ss, const std::unordered_map<size_t, position_description> &interesting_positions);

};
}
;
