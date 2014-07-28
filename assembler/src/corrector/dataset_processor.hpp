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
	size_t contig_length;
//TODO: place for multilib
	string sam_filename;
	vector<string> buffered_reads;

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
	//map<string, std::ofstream*> all_writers;
	int nthreads;
	int buffered_count ;
	const int buff_size = 1000000;
public:
	DatasetProcessor(string sam_file, string genome_file, string output_dir, string work_dir):sam_file(sam_file), genome_file(genome_file), work_dir(work_dir){
		//path::make_dir(work_dir);
		output_contig_file = output_dir + "/corrected_contigs.fasta";
		nthreads = corr_cfg::get().max_nthreads;
		buffered_count = 0;
	}
	void OutputRead(string &read, string &contig_name);
//	void PrepareWriters();
	void SplitGenome(const string &genome, const string &genome_splitted_dir);
	void FlushAll();
	void BufferedOutputRead(string &read, string &contig_name);
	void GetAlignedContigs(string &read, set<string> &contigs);
	void SplitSingleLibrary(string &out_contigs_filename);
	void SplitPairedLibrary(string &all_reads);
	void GlueSplittedContigs(string &out_contigs_filename);
	void ProcessSplittedLibrary();
	void ProcessLibrary(string &sam_file);
	void SplitHeaders(string &all_reads_filename);
//	void CloseWriters();
	void ProcessDataset();

};
};
