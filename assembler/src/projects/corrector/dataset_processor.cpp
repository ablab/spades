//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dataset_processor.hpp"
#include "variants_table.hpp"
#include "contig_processor.hpp"
#include "config_struct.hpp"

#include "io/reads/file_reader.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "io/reads/osequencestream.hpp"
#include "utils/parallel/openmp_wrapper.h"

#include <boost/algorithm/string.hpp>

#include <iostream>
#include <unistd.h>

using namespace std;

namespace corrector {
std::string DatasetProcessor::GetLibDir(const size_t lib_count) {
    if (lib_dirs_.find(lib_count) != lib_dirs_.end())
        return lib_dirs_[lib_count];
    std::string res = fs::make_temp_dir(corr_cfg::get().work_dir, "lib" + to_string(lib_count));
    lib_dirs_[lib_count] = res;
    return res;
}

void DatasetProcessor::SplitGenome(const string &genome_splitted_dir) {
    io::FileReadStream frs(genome_file_);
    size_t cur_id = 0;
    while (!frs.eof()) {
        io::SingleRead cur_read;
        frs >> cur_read;
        string contig_name = cur_read.name();
        string contig_seq = cur_read.GetSequenceString();
        if (all_contigs_.find(contig_name) != all_contigs_.end()) {
            WARN("Duplicated contig names! Multiple contigs with name" << contig_name);
        }
        string full_path = fs::append_path(genome_splitted_dir, contig_name + ".fasta");
        string out_full_path = fs::append_path(genome_splitted_dir, contig_name + ".ref.fasta");
        string sam_filename = fs::append_path(genome_splitted_dir, contig_name + ".pair.sam");
        all_contigs_[contig_name] = {full_path, out_full_path, contig_seq.length(), sam_files_type(), sam_filename, cur_id};
        cur_id ++;
        buffered_reads_[contig_name].clear();
        io::OFastaReadStream oss(full_path);
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

void DatasetProcessor::SplitLibrary(const string &all_reads_filename, const size_t lib_count, bool is_paired = false) {
    int reads_cnt = is_paired ? 2 : 1;
    ifstream fs(all_reads_filename);
    while (!fs.eof()) {
        set<string> contigs;
        std::vector<std::string> reads(reads_cnt);
        getline(fs, reads[0]);
        if (reads[0][0] == '@')
            continue;

        if (is_paired) {
            getline(fs, reads[1]);
        }

        for (int i = 0; i < reads_cnt; ++i) {
            GetAlignedContigs(reads[i], contigs);
        }

        for (auto &contig : contigs) {
            CHECK_FATAL_ERROR(all_contigs_.find(contig) != all_contigs_.end(),
                              "wrong contig name in SAM file header: " + contig);

            for (int i = 0; i < reads_cnt; ++i) {
                BufferedOutputRead(reads[i], contig, lib_count);
            }
        }
    }
    FlushAll(lib_count);
}

void DatasetProcessor::FlushAll(const size_t lib_count) {
    for (const auto &ac : all_contigs_) {
        if (buffered_reads_[ac.first].size() > 0) {
            ofstream stream(ac.second.sam_filenames[lib_count].first.c_str(), std::ios_base::app | std::ios_base::out);
            for (const string &read : buffered_reads_[ac.first]) {
                stream << read;
                stream << '\n';
            }
            buffered_reads_[ac.first].clear();
        }
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


int DatasetProcessor::RunBwaIndex() {
    string bwa_string = fs::screen_whitespaces(fs::screen_whitespaces(corr_cfg::get().bwa));
    string genome_screened = fs::screen_whitespaces(genome_file_);
    string index_line = bwa_string + " index " + genome_screened;
    INFO("Running bwa index ...: " << index_line);
    int run_res = system(index_line.c_str());
    if (run_res != 0) {
        INFO("bwa failed, skipping sublib");
    }

    return run_res;
}

std::string DatasetProcessor::RunBwaMem(const std::vector<std::string> &reads, const size_t lib,
    const std::string &params = "") {
    string cur_dir = GetLibDir(lib);
    string tmp_sam_filename = fs::append_path(cur_dir, "tmp.sam");
    string bwa_string = fs::screen_whitespaces(fs::screen_whitespaces(corr_cfg::get().bwa));
    string genome_screened = fs::screen_whitespaces(genome_file_);
    int run_res = 0;

    std::string reads_line = "";
    for (auto& filename : reads) {
        reads_line += fs::screen_whitespaces(filename) + " ";
    }

    string nthreads_str = to_string(nthreads_);
    string last_line = bwa_string + " mem  -v 1 -t " + nthreads_str + " " + params + " " + genome_screened + " " + reads_line  + "  > "
        + fs::screen_whitespaces(tmp_sam_filename) ;
    INFO("Running bwa mem ...:" << last_line);
    run_res = system(last_line.c_str());
    if (run_res != 0) {
        INFO("bwa failed, skipping sublib");
        return "";
    }
    return tmp_sam_filename;
}

void DatasetProcessor::PrepareContigDirs(const size_t lib_count) {
    string out_dir = GetLibDir(lib_count);
    for (auto &ac : all_contigs_) {
        auto contig_name = ac.first;
        string out_name = fs::append_path(out_dir, contig_name + ".sam");
        ac.second.sam_filenames.push_back(make_pair(out_name, unsplitted_sam_files_[lib_count].second));
        BufferedOutputRead("@SQ\tSN:" + contig_name + "\tLN:" + to_string(all_contigs_[contig_name].contig_length), contig_name, lib_count);
    }
    FlushAll(lib_count);
}

void DatasetProcessor::ProcessDataset() {
    size_t lib_num = 0;
    INFO("Splitting assembly...");
    INFO("Assembly file: " + genome_file_);
    SplitGenome(work_dir_);

    if (RunBwaIndex() != 0) {
        FATAL_ERROR("Failed to build bwa index for " << genome_file_);
    }

    auto handle_one_lib = [this, &lib_num](const std::vector<std::string>& reads,
        const std::string& type, const auto& lib_type){
        std::string reads_files_str = "";
        for (const auto& filename : reads) {
            reads_files_str += filename + " ";
        }

        INFO("Processing " + type + " sublib of number " << lib_num);
        INFO(reads_files_str);
        std::string param = "";
        if (type == "interlaced") {
            param = "-p";
        }

        string samf = RunBwaMem(reads, lib_num, param);
        if (samf != "") {
            INFO("Adding samfile " << samf);
            unsplitted_sam_files_.push_back(make_pair(samf, lib_type));
            PrepareContigDirs(lib_num);
            SplitLibrary(samf, lib_num,lib_type !=  io::LibraryType::SingleReads);
            lib_num++;
        } else {
            FATAL_ERROR("Failed to align " + type + " reads " << reads_files_str);
        }
    };

    for (size_t i = 0; i < corr_cfg::get().dataset.lib_count(); ++i) {
        const auto& dataset = corr_cfg::get().dataset[i];
        auto lib_type = dataset.type();
        if (lib_type == io::LibraryType::PairedEnd || lib_type == io::LibraryType::HQMatePairs || lib_type == io::LibraryType::SingleReads) {
            for (auto iter = dataset.paired_begin(); iter != dataset.paired_end(); iter++) {
                handle_one_lib({iter->first, iter->second}, "paired", lib_type);
            }

            for (auto iter = dataset.interlaced_begin(); iter != dataset.interlaced_end(); iter++) {
                handle_one_lib({*iter}, "interlaced", lib_type);
            }

            for (auto iter = dataset.single_begin(); iter != dataset.single_end(); iter++) {
                handle_one_lib({*iter}, "single", lib_type);
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
    vector<string> ordered_names;
    ordered_names.resize(all_contigs_.size());
    for (const auto &ac : all_contigs_) {
        ordered_names[ac.second.id] = ac.first;
    }
    for (size_t i = 0; i < ordered_names.size(); i++) {
        ifstream a_f(all_contigs_[ordered_names[i]].output_contig_filename, std::ios_base::binary);
        of_c << a_f.rdbuf();
    }
}

}
;
