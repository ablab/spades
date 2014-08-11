#include "dataset_processor.hpp"
#include "variants_table.hpp"
#include "contig_processor.hpp"
#include "config_struct.hpp"

#include "io/file_reader.hpp"
#include "path_helper.hpp"
#include "io/osequencestream.hpp"
#include "openmp_wrapper.h"

#include <iostream>
#include <unistd.h>

using namespace std;

namespace corrector {
static std::string GetLibDir(const size_t lib_count) {
    return path::append_path(corr_cfg::get().work_dir, "lib" + to_string(lib_count));
}

void DatasetProcessor::SplitGenome(const string &genome_splitted_dir) {
    io::FileReadStream frs(genome_file_);
    while (!frs.eof()) {
        io::SingleRead cur_read;
        frs >> cur_read;
        string contig_name = cur_read.name();
        string contig_seq = cur_read.GetSequenceString();
        if (all_contigs_.find(contig_name) != all_contigs_.end()) {
            WARN("Duplicated contig names! Multiple contigs with name" << contig_name);
        }
        string full_path = path::append_path(genome_splitted_dir, contig_name + ".fasta");
        string out_full_path = path::append_path(genome_splitted_dir, contig_name + ".ref.fasta");
        string sam_filename = path::append_path(genome_splitted_dir, contig_name + ".pair.sam");
        all_contigs_[contig_name] = {full_path, out_full_path, contig_seq.length(), sam_files_type(), sam_filename};
        buffered_reads_[contig_name].clear();
        io::osequencestream oss(full_path);
        oss << io::SingleRead(contig_name, contig_seq);
        DEBUG("full_path " + full_path)
    }
}

//contigs - set of aligned contig names
void DatasetProcessor::GetAlignedContigs(const string &read, set<string> &contigs) const {
    vector<string> arr;
    boost::split(arr, read, boost::is_any_of("\t"));
    if (arr.size() > 5) {
        if (arr[2] != "*" && stoi(arr[4]) > 0) {
// here can be multuple aligned parsing if neeeded;
            contigs.insert(arr[2]);
        }
    }

}

void DatasetProcessor::SplitSingleLibrary(const string &all_reads_filename, const size_t lib_count) {
    ifstream fs(all_reads_filename);
    while (!fs.eof()) {
        set<string> contigs;
        string r1;
        getline(fs, r1);
        if (r1[0] == '@')
            continue;
        GetAlignedContigs(r1, contigs);
        for (auto &contig : contigs) {
            VERIFY_MSG(all_contigs_.find(contig) != all_contigs_.end(), "wrong contig name in SAM file header: " + contig);
            BufferedOutputRead(r1, contig, lib_count);
        }
    }
    FlushAll(lib_count);
}

void DatasetProcessor::FlushAll(const size_t lib_count) {
    for (const auto &ac : all_contigs_) {
        ofstream stream(ac.second.sam_filenames[lib_count].first.c_str(), std::ios_base::app | std::ios_base::out);
        for (const string &read : buffered_reads_[ac.first]) {
            stream << read;
            stream << '\n';
        }
        buffered_reads_[ac.first].clear();
    }
}

void DatasetProcessor::BufferedOutputRead(const string &read, const string &contig_name, const size_t lib_count) {
    buffered_reads_[contig_name].push_back(read);
    buffered_count_++;
    if (buffered_count_ % kBuffSize == 0) {
        if (buffered_count_ % (10 * kBuffSize) == 0)
            INFO("processed " << buffered_count_ << "reads, flushing");
        FlushAll(lib_count);
    }
}

void DatasetProcessor::SplitPairedLibrary(const string &all_reads_filename, const size_t lib_count) {
    ifstream fs(all_reads_filename);
    while (!fs.eof()) {
        set<string> contigs;
        string r1;
        string r2;
        getline(fs, r1);
        if (r1[0] == '@')
            continue;
        getline(fs, r2);
        GetAlignedContigs(r1, contigs);
        GetAlignedContigs(r2, contigs);
        for (const auto &contig : contigs) {
            VERIFY_MSG(all_contigs_.find(contig) != all_contigs_.end(), "wrong contig name in SAM file header: " + contig);
            if (all_contigs_.find(contig) != all_contigs_.end()) {
                BufferedOutputRead(r1, contig, lib_count);
                BufferedOutputRead(r2, contig, lib_count);
            }
        }
    }
    FlushAll(lib_count);
}

string DatasetProcessor::RunPairedBwa(const string &left, const string &right, const size_t lib) const {
    string cur_dir = GetLibDir(lib);
    if (!path::make_dir(cur_dir)) {
        WARN("Failed to create directory " << cur_dir << ", skipping sublib");
        return "";
    }
    int run_res = 0;
    string tmp1_sai_filename = path::append_path(cur_dir, "tmp1.sai");
    string tmp2_sai_filename = path::append_path(cur_dir, "tmp2.sai");
    string tmp_sam_filename = path::append_path(cur_dir, "tmp.sam");
    string isize_txt_filename = path::append_path(cur_dir, "isize.txt");
    string tmp_file = path::append_path(cur_dir, "bwa.flood");

    string index_line = corr_cfg::get().bwa + " index " + "-a " + "is " + genome_file_ + " " + " 2>" + tmp_file;
    INFO("Running bwa index ...: " << index_line);
    run_res = system(index_line.c_str());
    if (run_res != 0) {
        INFO("bwa failed, skipping sublib");
        return "";
    }
    string nthreads_str = to_string(nthreads_);
    string left_line = corr_cfg::get().bwa + " aln " + genome_file_ + " " + left + " -t " + nthreads_str + " -O 7 -E 2 -k 3 -n 0.08 -q 15 > "
            + tmp1_sai_filename + " 2>" + tmp_file;

    string right_line = corr_cfg::get().bwa + " aln " + genome_file_ + " " + right + " -t " + nthreads_str + " -O 7 -E 2 -k 3 -n 0.08 -q 15 > "
            + tmp2_sai_filename + " 2>" + tmp_file;
    INFO("Running bwa aln ...: " << left_line);
    int res1 = system(left_line.c_str());
    INFO("Running bwa aln ...: " << right_line);
    int res2 = system(right_line.c_str());
    if (res1 != 0 || res2 != 0) {
        INFO("bwa failed, skipping sublib");
        return "";
    }

    string last_line = corr_cfg::get().bwa + " sampe " + genome_file_ + " " + tmp1_sai_filename + " " + tmp2_sai_filename + " " + left + " " + right + "  > "
            + tmp_sam_filename + " 2>" + isize_txt_filename;
    INFO("Running bwa sampe ...:" << last_line);
    run_res = system(last_line.c_str());
    if (run_res != 0) {
        INFO("bwa failed, skipping sublib");
        return "";
    }
    res1 = unlink(tmp1_sai_filename.c_str());
    res2 = unlink(tmp2_sai_filename.c_str());
    if (res1 != 0 || res2 != 0) {
        INFO("Failed to delete temporary sai files");
    }
    return tmp_sam_filename;
}

string DatasetProcessor::RunSingleBwa(const string &single, const size_t lib) const {
    int run_res = 0;
    string cur_dir = GetLibDir(lib);
    if (!path::make_dir(cur_dir)) {
        WARN("Failed to create directory " << cur_dir << ", skipping sublib");
        return "";
    }
    string tmp_sai_filename = path::append_path(cur_dir, "tmp1.sai");
    string tmp_sam_filename = path::append_path(cur_dir, "tmp.sam");
    string isize_txt_filename = path::append_path(cur_dir, "isize.txt");
    string tmp_file = path::append_path(cur_dir, "bwa.flood");

    string index_line = corr_cfg::get().bwa + " index " + "-a " + "is " + genome_file_ + " 2 " + "2>" + tmp_file;
    INFO("Running bwa index ...: " << index_line);
    run_res = system(index_line.c_str());
    if (run_res != 0) {
        INFO("bwa failed, skipping sublib");
        return "";
    }
    string nthreads_str = to_string(nthreads_);
    string single_sai_line = corr_cfg::get().bwa + " aln " + genome_file_ + " " + single + " -t " + nthreads_str + " -O 7 -E 2 -k 3 -n 0.08 -q 15 > "
            + tmp_sai_filename + " 2>" + tmp_file;

    INFO("Running bwa aln ...:" + single_sai_line);

    run_res = system(single_sai_line.c_str());
    if (run_res != 0) {
        INFO("bwa failed, skipping sublib");
        return "";
    }
    string last_line = corr_cfg::get().bwa + " samse " + genome_file_ + " " + tmp_sai_filename + " " + single + "  > " + tmp_sam_filename + " 2>"
            + isize_txt_filename;
    INFO("Running bwa samse ...:" << last_line);
    run_res = system(last_line.c_str());
    if (run_res != 0) {
        INFO("bwa failed, skipping sublib");
        return "";
    }
    run_res = unlink(tmp_sai_filename.c_str());
    if (run_res != 0) {
        INFO("Failed to delete temporary sai files");
    }
    return tmp_sam_filename;
}

void DatasetProcessor::PrepareContigDirs(const size_t lib_count) {
    string out_dir = GetLibDir(lib_count);
    for (auto &ac : all_contigs_) {
        auto contig_name = ac.first;
        string header = "@SQ\tSN:" + contig_name + "\tLN:" + to_string(all_contigs_[contig_name].contig_length);
        BufferedOutputRead(header, contig_name, lib_count);
        string out_name = path::append_path(out_dir, ac.first + ".sam");
        ac.second.sam_filenames.push_back(make_pair(out_name, unsplitted_sam_files_[lib_count].second));
    }
    FlushAll(lib_count);
}

void DatasetProcessor::ProcessDataset() {
    size_t lib_num = 0;
    INFO("Splitting assembly...");
    INFO("Assembly file: " + genome_file_);
    SplitGenome(work_dir_);
    for (size_t i = 0; i < corr_cfg::get().dataset.lib_count(); ++i) {
        const auto& dataset = corr_cfg::get().dataset[i];
        auto lib_type = dataset.type();
        if (lib_type == io::LibraryType::PairedEnd || lib_type == io::LibraryType::HQMatePairs || lib_type == io::LibraryType::SingleReads) {
            for (auto iter = dataset.paired_begin(); iter != dataset.paired_end(); iter++) {
                INFO("Processing paired sublib of number " << lib_num);
                string left = iter->first;
                string right = iter->second;
                INFO(left + " " + right);
                string samf = RunPairedBwa(left, right, lib_num);
                if (samf != "") {
                    unsplitted_sam_files_.push_back(make_pair(samf, lib_type));
                    PrepareContigDirs(lib_num);
                    SplitPairedLibrary(samf, lib_num);
                    lib_num++;
                } else {
                    WARN("Failed to align paired reads " << left << " and " << right);
                }
            }
            for (auto iter = dataset.single_begin(); iter != dataset.single_end(); iter++) {
                INFO("Processing single sublib of number " << lib_num);
                string left = *iter;
                INFO(left);
                string samf = RunSingleBwa(left, lib_num);
                if (samf != "") {
                    INFO("Adding samfile " << samf);
                    unsplitted_sam_files_.push_back(make_pair(samf, io::LibraryType::SingleReads));
                    PrepareContigDirs(lib_num);
                    SplitSingleLibrary(samf, lib_num);
                    lib_num++;
                } else {
                    WARN("Failed to align single reads " << left);
                }
            }
        }
    }
    INFO("Processing contigs");
    vector<pair<size_t, string> > ordered_contigs;
    for (const auto &ac : all_contigs_) {
        ordered_contigs.push_back(make_pair(ac.second.contig_length, ac.first));
    }
    size_t cont_num = ordered_contigs.size();
    sort(ordered_contigs.begin(), ordered_contigs.end(), std::greater<pair<size_t, string> >());
    auto all_contigs_ptr = &all_contigs_;
# pragma omp parallel for shared(all_contigs_ptr, ordered_contigs) num_threads(nthreads_) schedule(dynamic,1)
    for (size_t i = 0; i < cont_num; i++) {
        bool long_enough = (*all_contigs_ptr)[ordered_contigs[i].second].contig_length > kMinContigLengthForInfo;
        ContigProcessor pc((*all_contigs_ptr)[ordered_contigs[i].second].sam_filenames, (*all_contigs_ptr)[ordered_contigs[i].second].input_contig_filename);
        size_t changes = pc.ProcessMultipleSamFiles();
        if (long_enough) {
#pragma omp critical
            {
                INFO("Contig " << ordered_contigs[i].second << " processed with " << changes << " changes in thread " << omp_get_thread_num());
            }
        }
    }
    INFO("Gluing processed contigs");
    GlueSplittedContigs(output_contig_file_);
}

void DatasetProcessor::GlueSplittedContigs(string &out_contigs_filename) {
    ofstream of_c(out_contigs_filename, std::ios_base::binary);
    for (const auto &ac : all_contigs_) {
        ifstream a_f(ac.second.output_contig_filename, std::ios_base::binary);
        of_c << a_f.rdbuf();
    }
}

}
;
