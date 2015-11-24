#include <cstdio>
#include <string.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <io/mmapped_reader.hpp>
#include <libcxx/sort.hpp>
#include "kmc_api/kmc_file.h"
#include <runtime_k.hpp>
#include <boost/optional/optional.hpp>
#include <boost/lexical_cast.hpp>

using std::string;
using std::vector;
using std::shared_ptr;
using std::make_shared;
using runtime_k::RtSeq;
using boost::lexical_cast;

const string KMER_PARSED_EXTENSION = ".bin";
const string KMC_EXTENSION = ".kmc";
const string KMER_SORTED_EXTENSION = ".sorted";

void PrintUsageInfo() {
    std::cout << "Usage: kmer_multiplicity_counter [options] <sample_desc_1> [<sample_desc_2> ...]" << std::endl;
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
    cmd << "kmc -k" << k << " -ci" << min_mult << " @" << filename << " " << filename << KMC_EXTENSION << " .";
    system(cmd.str().c_str());
    ParseKmc(filename, k);
    system(("rm -f " + filename + KMC_EXTENSION + "*").c_str());
    SortKmersCountFile(filename, k);
    system(("rm -f " + filename + KMER_PARSED_EXTENSION).c_str());
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

//    std::ofstream output_kmer(output_file + ".kmer", std::ios::binary);
//    std::ofstream output_cnt(output_file + ".mpl", std::ios::binary);
    std::ofstream output_kmer(output_file + ".kmer");
    std::ofstream output_cnt(output_file + ".mpl");
    RtSeq::less3 kmer_less;
    while (true) {
        boost::optional<RtSeq> min_kmer;
        size_t cnt_min = 0;
        for (size_t i = 0; i < n; ++i) {
            if (alive[i]) {
                RtSeq& cur_kmer = top_kmer[i].first;
//                bool t1, t2;
//                if (min_kmer) {
//                    t1 = kmer_less(cur_kmer, *min_kmer);
//                    t2 = cur_kmer == *min_kmer;
//                }
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
        std::vector<uint32> cnt_vector(n, 0);
        if (cnt_min >= all_min) {
//            min_kmer.get().BinWrite(output_kmer);
            output_kmer << min_kmer.get().str() << std::endl;
            for (size_t i = 0; i < n; ++i) {
                if (i != 0) {
                    output_cnt << " ";
                }
                if (alive[i] && top_kmer[i].first == *min_kmer) {
                    cnt_vector[i] += top_kmer[i].second;
                }
            }
        }
        string delim = "";
        for (auto cnt : cnt_vector) {
            output_cnt << delim << cnt;
            delim = " ";
        }
        output_cnt << std::endl;
        for (size_t i = 0; i < n; ++i) {
            if (alive[i] && top_kmer[i].first == *min_kmer) {
                alive[i] = ReadKmerWithCount(*infiles[i], top_kmer[i]);
            }
        }
    }
    output_kmer.close();
    output_cnt.close();
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        PrintUsageInfo();
        return 0;
    }
    size_t min_mult = -1ul;
    size_t min_sample_count = -1ul;
    size_t k = -1ul;
    string output = "";
    int i = 1;
    while (i < argc && argv[i][0] == '-') {
        VERIFY(i + 1 < argc);

        string arg = argv[i];
        if (arg == "--mult") {
            min_mult = lexical_cast<size_t>(argv[i + 1]);
        }
        if (arg == "--sample") {
            min_sample_count = lexical_cast<size_t>(argv[i + 1]);
        }
        if (arg == "-k") {
            k = lexical_cast<size_t>(argv[i + 1]);
        }
        if (arg == "-o") {
            output = argv[i + 1];
        }
        i += 2;
    }
    if (argc - i < 3 || (argc - i - 1) % 2 != 0
        || min_mult == -1ul
        || min_sample_count == -1ul
        || k == -1ul
        || output == "") {
        PrintUsageInfo();
        return 0;
    }

    std::vector<string> input_files;
    for (; i < argc; ++i) {
        input_files.push_back(argv[i]);
    }

    vector<string> kmer_cnt_files;
    for (const auto& file : input_files) {
        kmer_cnt_files.push_back(CountSample(file, min_mult, k));
    }

    FilterKmers(kmer_cnt_files, min_sample_count, k, output);

    for (const auto& file: kmer_cnt_files) {
        remove(file.c_str());
    }
    return 0;
}
