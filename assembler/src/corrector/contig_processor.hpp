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

namespace corrector {
class ContigProcessor {
	string sam_file;
	string contig_file;
	string contig_name;
	string output_contig_file;
//TODO: readlength?
	string contig;

	MappedSamStream sm;

	bam_header_t *bam_header;
	vector<position_description> charts;
	InterestingPositionProcessor ipp;
	vector<int> error_counts;
public:
	ContigProcessor(string sam_file, string contig_file):sam_file(sam_file), contig_file(contig_file), sm(sam_file){
		INFO("CP creating..");
		bam_header = sm.ReadHeader();
		read_contig();
		ipp.set_contig(contig);
		INFO("CP created..");
	}
	void read_contig();
	void UpdateOneRead(SingleSamRead &tmp);
	//returns: number of changed nucleotides;
	int UpdateOneBase(size_t i, stringstream &ss, unordered_map<size_t, position_description> &interesting_positions);
	void process_sam_file ();

	//string seq, cigar; pair<size_t, size_t> borders;
};
};
