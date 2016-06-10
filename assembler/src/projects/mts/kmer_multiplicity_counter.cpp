#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <set>
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
#include "dev_support/simple_tools.hpp"

using std::string;
using std::vector;
using std::shared_ptr;
using std::make_shared;
using runtime_k::RtSeq;

const string KMER_PARSED_EXTENSION = ".bin";
const string KMER_SORTED_EXTENSION = ".sorted";

void PrintUsageInfo() {
    std::cout << "Usage: kmer_multiplicity_counter [options] -f (<sampleX.kmc>)+" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "-k - kmer length" << std::endl;
    std::cout << "-o - output file" << std::endl;
    std::cout << "--sample - minimal number of samples to contain kmer" << std::endl;
    std::cout << "sampleX.kmc is a file with kmer multiplicities of sample X" << std::endl;
}

//TODO: get rid of intermediate .bin file
string ParseKmc(const string& filename, size_t k) {
    CKMCFile kmcFile;
    kmcFile.OpenForListing(filename);
    CKmerAPI kmer((unsigned int) k);
    uint32 count;
    std::string parsed_filename = filename + KMER_PARSED_EXTENSION;
    std::ofstream output(parsed_filename, std::ios::binary);
    while (kmcFile.ReadNextKmer(kmer, count)) {
        RtSeq seq(k, kmer.to_string());
        seq.BinWrite(output);
        seq_element_type tmp = count;
        output.write((char*) &(tmp), sizeof(seq_element_type));
    }
    output.close();
    remove((filename + ".kmc_pre").c_str());
    remove((filename + ".kmc_suf").c_str());
    return parsed_filename;
}

string SortKmersCountFile(const string& filename, const size_t k) {
    MMappedRecordArrayReader<seq_element_type> ins(filename, RtSeq::GetDataSize(k) + 1, false);
    libcxx::sort(ins.begin(), ins.end(), array_less<seq_element_type>());
    std::string sorted_filename = filename + KMER_SORTED_EXTENSION;
    std::ofstream out(sorted_filename);
    out.write((char*) ins.data(), ins.data_size());
    out.close();
    remove(filename.c_str());
    return sorted_filename;
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
        auto parsed = ParseKmc(fn, k);
        auto sorted = SortKmersCountFile(parsed, k);
        infiles.push_back(std::make_shared<std::ifstream>(sorted));
    }
    vector<std::pair<RtSeq, uint32>> top_kmer(n, {RtSeq(k), 0});
    vector<bool> alive(n, false);

    for (size_t i = 0; i < n; i++) {
        alive[i] = ReadKmerWithCount(*infiles[i], top_kmer[i]);
    }

    std::ofstream output_kmer(output_file + ".kmer", std::ios::binary);
    std::ofstream output_cnt(output_file + ".mpl");
    //std::ofstream output_full(output_file + ".kmpl");
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
            //output_full << min_kmer->str() << " ";
            for (size_t i = 0; i < n; ++i) {
                if (alive[i] && top_kmer[i].first == *min_kmer) {
                    cnt_vector[i] += top_kmer[i].second;
                }
            }
            string delim = "";
            for (auto cnt : cnt_vector) {
                output_cnt << delim << cnt;
                //output_full << delim << cnt;
                delim = " ";
            }
            output_cnt << std::endl;
            //output_full << std::endl;
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

    size_t min_sample_count = -1ul;
    size_t k = -1ul;
    string output;
    std::vector<string> input_files;

    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('s', "sample", min_sample_count)
            >> Option('o', output)
            >> Option('f', input_files)
        ;
    } catch(GetOptEx &ex) {
        PrintUsageInfo();
        exit(1);
    }

    //FIXME temporary workaround!
    VERIFY(!input_files.empty());
    std::string dir = path::parent_path(input_files.front());
    std::vector<string> new_input_files;
    for (size_t i = 1 ; i <= input_files.size() ; ++i) {
        new_input_files.push_back(dir + "/sample" + ToString(i));
    }

    for (const auto& ifl : new_input_files) {
        cout << "Input file " << ifl << endl;
    }

    FilterKmers(new_input_files, min_sample_count, k, output);

    return 0;
}
