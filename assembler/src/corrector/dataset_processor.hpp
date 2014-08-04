#pragma once

// WTF: get rid of include-all-you-can-ever-imagine-header
#include "include.hpp"
#include "contig_processor.hpp"
#include "sam_reader.hpp"
#include "read.hpp"
#include "interesting_pos_processor.hpp"
#include "path_helper.hpp"

// FIXME: EVERYWHERE: USE SPACES, NOT TABS! FIX ALL THE CODING STYLE PROBLEMS EVERYWHERE

namespace corrector {

struct OneContigDescription{
	string input_contig_filename;
	string output_contig_filename;
	size_t contig_length;
	sam_files_type sam_filenames;
	string sam_filename;
	vector<string> buffered_reads;

};
typedef unordered_map<string, OneContigDescription> ContigInfoMap;

class DatasetProcessor {
	string genome_file;
	string output_contig_file;
	ContigInfoMap all_contigs;
	vector<position_description> charts;
	InterestingPositionProcessor ipp;
	vector<int> error_counts;
	sam_files_type unsplitted_sam_files;
	string work_dir;
    // WTF: Why these vars (3 of them) int? Can you have -1 threads?
	int nthreads;
	int buffered_count ;
	const int buff_size = 100000;
public:
	DatasetProcessor(string genome_file): genome_file(genome_file), work_dir(corr_cfg::get().work_dir){
		//path::make_dir(work_dir);
        // WTF: Use stuff from path to form filenames properly
		output_contig_file = corr_cfg::get().output_dir + "/corrected_contigs.fasta";
		nthreads = corr_cfg::get().max_nthreads;
		buffered_count = 0;
	}
	void ProcessDataset();
private:
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

};
};
