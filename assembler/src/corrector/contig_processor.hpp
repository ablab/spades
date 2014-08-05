/*
 * contig_processor.hpp
 *
 *  Created on: Jun 27, 2014
 *      Author: lab42
 */

#pragma once
#include "sam_reader.hpp"
#include "read.hpp"
#include "include.hpp"
#include "interesting_pos_processor.hpp"
#include "positional_read.hpp"

#include <omp.h>

namespace corrector {
typedef std::vector<std::pair<string, io::LibraryType> > sam_files_type;
class ContigProcessor {
    sam_files_type sam_files;
    string sam_file;
    string contig_file;
    string contig_name;
    string output_contig_file;
    string contig;
    bool debug_info;
    vector<position_description> charts;
    InterestingPositionProcessor ipp;
    vector<int> error_counts;
public:
    ContigProcessor(sam_files_type &sam_files_, string contig_file)
            : contig_file(contig_file) {
        ReadContig();
        ipp.set_contig(contig);
        debug_info = (contig.length() > 20000);
        for (auto sf : sam_files_)
            sam_files.push_back(sf);
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
    int UpdateOneBase(size_t i, stringstream &ss, const unordered_map<size_t, position_description> &interesting_positions);

};
}
;
