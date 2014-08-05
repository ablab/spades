#pragma once

#include "include.hpp"
//#include "io/file_reader.hpp"
#include "contig_processor.hpp"

#include "path_helper.hpp"


namespace corrector {

struct OneContigDescription {
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
    size_t nthreads;
    size_t buffered_count;
    const size_t buff_size = 100000;
public:
    DatasetProcessor(string genome_file)
            : genome_file(genome_file),
              work_dir(corr_cfg::get().work_dir) {
        output_contig_file = path::append_path(corr_cfg::get().output_dir , "corrected_contigs.fasta");
        nthreads = corr_cfg::get().max_nthreads;
        buffered_count = 0;
    }
    void ProcessDataset();
private:
    void SplitGenome(const string &genome, const string &genome_splitted_dir);
    void FlushAll(const size_t lib_count);
    void BufferedOutputRead(const string &read, const string &contig_name, const size_t lib_count);
    void GetAlignedContigs(const string &read, set<string> &contigs) const;
    std::string GetLibDir(const size_t lib_count) const;
    void SplitSingleLibrary(const string &out_contigs_filename, const size_t lib_count);
    void SplitPairedLibrary(const string &all_reads, const size_t lib_count);
    void GlueSplittedContigs(string &out_contigs_filename);
    std::string RunPairedBwa(const string &left, const string &right, const size_t lib) const;
    std::string RunSingleBwa(const string &single, const size_t lib) const;
    void PrepareContigDirs(const size_t lib_count);

};
}
;
