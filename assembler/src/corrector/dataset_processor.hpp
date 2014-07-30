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
	sam_files_type sam_filenames;
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
	sam_files_type unsplitted_sam_files;
	string work_dir;
	//map<string, std::ofstream*> all_writers;
	int nthreads;
	int buffered_count ;
	const int buff_size = 100000;
public:
	DatasetProcessor(string sam_file, string genome_file, string output_dir, string work_dir):sam_file(sam_file), genome_file(genome_file), work_dir(work_dir){
		//path::make_dir(work_dir);
		output_contig_file = output_dir + "/corrected_contigs.fasta";
		nthreads = corr_cfg::get().max_nthreads;
		buffered_count = 0;
	}
	//void OutputRead(string &read, string &contig_name);
//	void PrepareWriters();
	void SplitGenome(const string &genome, const string &genome_splitted_dir);
	void FlushAll(const size_t lib_count) ;
	void BufferedOutputRead(const string &read, const string &contig_name, const size_t lib_count) ;
	void GetAlignedContigs(const string &read, set<string> &contigs) const ;
	std::string GetLibDir(const size_t lib_count) const;
	void SplitSingleLibrary(const string &out_contigs_filename, const size_t lib_count);
	void SplitPairedLibrary(const string &all_reads, const size_t lib_count);
	void GlueSplittedContigs(string &out_contigs_filename);
	std::string RunPairedBwa(const string &left, const string &right, const size_t lib) const;
	std::string RunSingleBwa(const string &single, const size_t lib) const;
	void PrepareContigDirs(const size_t lib_count);
//	void ProcessLibrary(string &sam_file);
	//void SplitHeaders(string &all_reads_filename);
//	void CloseWriters();
	void ProcessDataset();

};
};
