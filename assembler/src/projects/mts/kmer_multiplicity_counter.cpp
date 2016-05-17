#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <io/kmers_io/mmapped_reader.hpp>
#include <libcxx/sort.hpp>
#include "getopt_pp/getopt_pp.h"
#include "kmc_api/kmc_file.h"
//#include "omp.h"
#include <data_structures/sequence/runtime_k.hpp>
#include <boost/optional/optional.hpp>
#include "dev_support/path_helper.hpp"

using std::string;
using std::vector;
using std::shared_ptr;
using std::make_shared;
using runtime_k::RtSeq;

const string KMER_PARSED_EXTENSION = ".bin";
const string KMC_EXTENSION = ".kmc";
const string KMER_SORTED_EXTENSION = ".sorted";

void PrintUsageInfo() {
    std::cout << "Usage: kmer_multiplicity_counter [options] -d (<sample_desc_X>)+" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "-k - kmer length" << std::endl;
    std::cout << "-o - output file" << std::endl;
    std::cout << "--mult - minimal occurence per file" << std::endl;
    std::cout << "--sample - minimal number of samples to contain kmer" << std::endl;
    std::cout << "sample_desc_X is a file with list of fastq files (gzipped or not) related to sample X" << std::endl;
}

void ParseKmc(const string& filename, size_t k) {
    CKMCFile kmcFile;
    kmcFile.OpenForListing(filename + KMC_EXTENSION);
    CKmerAPI kmer((unsigned int) k);
    uint32 count;

    std::ofstream output(filename + KMER_PARSED_EXTENSION, std::ios::binary);
    while (kmcFile.ReadNextKmer(kmer, count)) {
        RtSeq seq(k, kmer.to_string());
        seq.BinWrite(output);
        seq_element_type tmp = count;
        output.write((char*) &(tmp), sizeof(seq_element_type));
    }
    output.close();
}

void SortKmersCountFile(const string& filename, const size_t k) {
    MMappedRecordArrayReader<seq_element_type> ins(filename + KMER_PARSED_EXTENSION, RtSeq::GetDataSize(k) + 1, false);
    libcxx::sort(ins.begin(), ins.end(), array_less<seq_element_type>());
    std::ofstream out(filename + KMER_SORTED_EXTENSION);
    out.write((char*) ins.data(), ins.data_size());
    out.close();
}

string CountSample(const string& filename, size_t min_mult, size_t k) {
    stringstream cmd;
    std::string tmp_dir = filename + "_tmp";
    path::make_dir(tmp_dir);
    cmd << "kmc -k" << k << " -t8 -ci" << min_mult << " @" << filename << " " << filename << KMC_EXTENSION << " " << tmp_dir;
    system(cmd.str().c_str());
    ParseKmc(filename, k);
    system(("rm -f " + filename + KMC_EXTENSION + "*").c_str());
    SortKmersCountFile(filename, k);
    system(("rm -f " + filename + KMER_PARSED_EXTENSION).c_str());
    system(("rm -rf " + tmp_dir).c_str());
    return filename + KMER_SORTED_EXTENSION;
}

bool ReadKmerWithCount(std::ifstream& infile, std::pair<RtSeq, uint32>& res) {
    RtSeq seq(res.first.size());
    if (!seq.BinRead(infile)) {
        return false;
    }
    seq_element_type tmp;
    infile.read((char*) &tmp, sizeof(seq_element_type));
    res = {seq, (uint32) tmp};
    return true;
}

void FilterKmers(const std::vector<string>& files, size_t all_min, size_t k, const string& output_file) {
    size_t n = files.size();
    vector<shared_ptr<std::ifstream>> infiles;
    for (auto fn : files) {
        infiles.push_back(std::make_shared<std::ifstream>(fn));
    }
    vector<std::pair<RtSeq, uint32>> top_kmer(n, {RtSeq(k), 0});
    vector<bool> alive(n, false);

    for (size_t i = 0; i < n; i++) {
        alive[i] = ReadKmerWithCount(*infiles[i], top_kmer[i]);
    }

    std::ofstream output_kmer(output_file + ".kmer", std::ios::binary);
    //std::ofstream output_cnt(output_file + ".mpl", std::ios::binary);
    std::ofstream output_cnt(output_file + ".mpl");
    RtSeq::less3 kmer_less;
    while (true) {
        boost::optional<RtSeq> min_kmer;
        size_t cnt_min = 0;
        for (size_t i = 0; i < n; ++i) {
            if (alive[i]) {
                RtSeq& cur_kmer = top_kmer[i].first;
                if (!min_kmer || kmer_less(cur_kmer, *min_kmer)) {
                    min_kmer = boost::optional<RtSeq>(cur_kmer);
                    cnt_min = 0;
                }
                if (cur_kmer == *min_kmer) {
                    cnt_min++;
                }
            }
        }
        if (!min_kmer) {
            break;
        }
        if (cnt_min >= all_min) {
            std::vector<uint32> cnt_vector(n, 0);
            min_kmer.get().BinWrite(output_kmer);
//            output_kmer << min_kmer.get().str() << std::endl;
            for (size_t i = 0; i < n; ++i) {
                if (alive[i] && top_kmer[i].first == *min_kmer) {
                    cnt_vector[i] += top_kmer[i].second;
                }
            }
            string delim = "";
            for (auto cnt : cnt_vector) {
                output_cnt << delim << cnt;
                delim = " ";
            //    output_cnt << cnt;
            }
             output_cnt << std::endl;
        }
        for (size_t i = 0; i < n; ++i) {
            if (alive[i] && top_kmer[i].first == *min_kmer) {
                alive[i] = ReadKmerWithCount(*infiles[i], top_kmer[i]);
            }
        }
    }
}

int main(int argc, char *argv[]) {
    using namespace GetOpt;

    size_t min_mult = -1ul;
    size_t min_sample_count = -1ul;
    size_t k = -1ul;
    string output;
    std::vector<string> input_files;

    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('m', "mult", min_mult)
            >> Option('s', "sample", min_sample_count)
            >> Option('o', output)
            >> Option('d', input_files)
        ;
    } catch(GetOptEx &ex) {
        PrintUsageInfo();
        exit(1);
    }

    vector<string> kmer_cnt_files(input_files.size());
    #pragma omp parallel for num_threads(4)
    for (size_t i = 0; i < input_files.size(); i++) {
        std::cout << "Processing " << input_files[i] << std::endl;
        kmer_cnt_files[i] = CountSample(input_files[i], min_mult, k);
    }

    FilterKmers(kmer_cnt_files, min_sample_count, k, output);

    for (const auto& file: kmer_cnt_files) {
        remove(file.c_str());
    }
    return 0;
}
