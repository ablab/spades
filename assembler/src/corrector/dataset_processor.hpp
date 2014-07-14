#pragma once

#include "include.hpp"
#include "contig_processor.hpp"
#include "sam_reader.hpp"
#include "read.hpp"
#include "interesting_pos_processor.hpp"
#include "path_helper.hpp"
namespace corrector {

struct OneContigDescription{
	string input_contig_filename;
	string output_contig_filename;
//TODO: place for multilib
	string sam_filename;
};
typedef unordered_map<string, OneContigDescription> ContigInfoMap;

class DatasetProcessor {
	string sam_file;
	string genome_file;
	string output_contig_file;
//TODO: readlength?
	ContigInfoMap all_contigs;
//	MappedSamStream sm;
	vector<position_description> charts;
	InterestingPositionProcessor ipp;
	vector<int> error_counts;
	string work_dir;
	map<string, std::ofstream*> all_writers;
public:
	DatasetProcessor(string sam_file, string genome_file):sam_file(sam_file), genome_file(genome_file){
		work_dir = "/home/lab42/work/someshit";
		//path::make_dir(work_dir);
		output_contig_file = "corrected.fasta";
	}
	void OutputRead(string &read, string &contig_name);
	void PrepareWriters(ContigInfoMap &all_contigs);
	void SplitGenome(const string &genome, const string &genome_splitted_dir, ContigInfoMap &all_contigs);
	void GetAlignedContigs(string &read, set<string> &contigs);
	void SplitSingleLibrary();
	void SplitPairedLibrary(string &all_reads, ContigInfoMap &all_contigs);
	void GlueSplittedContigs(string &out_contigs_filename, ContigInfoMap &all_contigs);
	void ProcessSplittedLibrary();
	void ProcessLibrary(string &sam_file);
	void SplitHeaders(string &all_reads_filename, ContigInfoMap &all_contigs);
	void CloseWriters(ContigInfoMap &all_contigs);
};
};
