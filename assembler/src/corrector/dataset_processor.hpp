#pragma once

#include "path_helper.hpp"

#include "io/file_reader.hpp"
#include "path_helper.hpp"
#include "io/library.hpp"

#include <string>
#include <set>
#include <vector>
#include <unordered_map>

namespace corrector {

typedef std::vector<std::pair<std:: string, io::LibraryType> > sam_files_type;

struct OneContigDescription {
    std::string input_contig_filename;
    std::string output_contig_filename;
    size_t contig_length;
    sam_files_type sam_filenames;
    std::string sam_filename;
};
typedef std::unordered_map<std::string, OneContigDescription> ContigInfoMap;

class DatasetProcessor {

    const std::string &genome_file_;
    std::string output_contig_file_;
    ContigInfoMap all_contigs_;
    sam_files_type unsplitted_sam_files_;
    const std::string &work_dir_;
    std::unordered_map<std::string, std::vector<std::string> > buffered_reads_;
    size_t nthreads_;
    size_t buffered_count_;

    const size_t kBuffSize = 100000;
    const size_t kMinContigLengthForInfo = 20000;
public:
    DatasetProcessor(const std::string &genome_file, const std::string &work_dir, const std::string &output_dir, const size_t &thread_num)
            : genome_file_(genome_file), work_dir_(work_dir), nthreads_(thread_num) {
        output_contig_file_ = path::append_path(output_dir, "corrected_contigs.fasta");
        buffered_count_ = 0;
    }

    void ProcessDataset();
private:
    void SplitGenome(const std::string &genome_splitted_dir);
    void FlushAll(const size_t lib_count);
    void BufferedOutputRead(const std::string &read, const std::string &contig_name, const size_t lib_count);
    void GetAlignedContigs(const std::string &read, std::set<std::string> &contigs) const;
    void SplitSingleLibrary(const std::string &out_contigs_filename, const size_t lib_count);
    void SplitPairedLibrary(const std::string &all_reads, const size_t lib_count);
    void GlueSplittedContigs(std::string &out_contigs_filename);
    std::string RunPairedBwa(const std::string &left, const std::string &right, const size_t lib) const;
    std::string RunSingleBwa(const std::string &single, const size_t lib) const;
    void PrepareContigDirs(const size_t lib_count);

};
}
;
