/*
 * contig_processor.hpp
 *
 *  Created on: Jun 27, 2014
 *      Author: lab42
 */

#pragma once
#include "sam_reader.hpp"
#include "read.hpp"
#include <omp.h>
#include "include.hpp"
#include "interesting_pos_processor.hpp"

namespace corrector {
	typedef std::vector<std::pair<string, string> > sam_files_type;
class ContigProcessor {
	sam_files_type sam_files;
	string sam_file;
	string contig_file;
	string contig_name;
	string output_contig_file;
//TODO: readlength?
	string contig;
	bool debug_info;
	//MappedSamStream sm;

	//bam_header_t *bam_header;
	vector<position_description> charts;
	InterestingPositionProcessor ipp;
	vector<int> error_counts;
public:
/*	ContigProcessor(string sam_file, string contig_file):sam_file(sam_file), contig_file(contig_file), sm(sam_file){
		bam_header = sm.ReadHeader();
		read_contig();
		ipp.set_contig(contig);
	}*/
	ContigProcessor(sam_files_type &sam_files_, string contig_file): contig_file(contig_file){
		//bam_header = sm.ReadHeader();
		read_contig();
		ipp.set_contig(contig);
		debug_info = (contig.length() > 20000);
		for (auto sf : sam_files_)
			sam_files.push_back(sf);
	}
	void read_contig();
	void UpdateOneRead(const SingleSamRead &tmp, MappedSamStream &sm);
	//returns: number of changed nucleotides;
	int UpdateOneBase(size_t i, stringstream &ss, const unordered_map<size_t, position_description> &interesting_positions);
	void process_sam_file ();
	void process_multiple_sam_files();
	//string seq, cigar; pair<size_t, size_t> borders;
};
};
