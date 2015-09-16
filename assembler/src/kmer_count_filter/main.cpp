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

using std::string;
using std::vector;
using std::shared_ptr;
using std::make_shared;
using runtime_k::RtSeq;

const string DEFAULT_ONE_MIN = "3";
const string DEFAULT_ALL_MIN = "3";
const string DEFAULT_KMER_LENGTH = "31";
const string KMER_PARSED_EXTENSION = ".bin";
const string KMER_SORTED_EXTENSION = ".sorted";

void PrintUsageInfo() {
    puts("Usage: kmer-abundance-filter [options] <read 1> <read 2> [<read 1> <read 2> ...] <output file>");
    puts("Options:");
    puts("-one-min - minimal kmer abundance in each file");
    puts("-all-min - minimal kmer abundance in all files");
    puts("-kl - kmer length");
}

void parseKmc(const string& filename, const size_t k) {
    std::ofstream output(filename + KMER_PARSED_EXTENSION, std::ios::binary);

    CKMCFile kmcFile;
    kmcFile.OpenForListing(filename);
    CKmerAPI kmer(k);
    uint32 count;
    while (kmcFile.ReadNextKmer(kmer, count)) {
        RtSeq seq(k, kmer.to_string().c_str());
        seq.BinWrite(output);
        seq_element_type tmp = count;
        output.write((char*) &(tmp), sizeof(seq_element_type));
    }

    output.close();
}

void SortKmersCountFile(const string &filename, const size_t k) {
    MMappedRecordArrayReader<seq_element_type > ins(filename + KMER_PARSED_EXTENSION, RtSeq::GetDataSize(k) + 1, false);
    libcxx::sort(ins.begin(), ins.end(), array_less<seq_element_type>());
    std::ofstream out(filename + KMER_SORTED_EXTENSION);
    out.write((char*) ins.data(), ins.data_size());
    out.close();
}

void CountKmersOne(const string& filename, const string& abundance_min, const string& kmer_length) {
    system(("~/libs/kmc/kmc -k" + kmer_length + " -ci" + abundance_min + " " + filename + " " + filename + " .").c_str());
    parseKmc(filename, (size_t) atoi(kmer_length.c_str()));
    system(("rm -f " + filename + ".kmc_pre " + filename + ".kmc_suf").c_str());
    SortKmersCountFile(filename, (const size_t) atoi(kmer_length.c_str()));
    system(("rm -f " + filename + KMER_PARSED_EXTENSION).c_str());
}

bool ReadKmerWithCount(std::ifstream &infile, std::pair<RtSeq, uint32> &res) {
    RtSeq seq(res.first.size());
    if (!seq.BinRead(infile)) {
        return false;
    }
    seq_element_type tmp;
    infile.read((char*) &tmp, sizeof(seq_element_type));
    res = {seq, (uint32) tmp};
    return true;
}

void FilterKmersAll(const std::vector<string>& files, int all_min, size_t k, const string& output_file) {
    size_t n = files.size();
    vector<shared_ptr<std::ifstream>> infiles;
    for (size_t i = 0; i < n; ++i) {
        string name = files[i] + KMER_SORTED_EXTENSION;
        infiles.push_back(std::make_shared<std::ifstream>(name));
    }
    vector<std::pair<RtSeq, uint32>> top_kmer(n, {RtSeq(k), 0});
    vector<bool> alive(n);

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
        int cnt_min = 0;
        for (size_t i = 0; i < n; ++i) {
            if (alive[i]) {
                RtSeq& cur_kmer = top_kmer[i].first;
                bool t1, t2;
                if (min_kmer) {
                    t1 = kmer_less(cur_kmer, *min_kmer);
                    t2 = cur_kmer == *min_kmer;
                }
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
        std::vector<uint32> cnt_vector(n / 2, 0);
        if (cnt_min >= all_min) {
//            min_kmer.get().BinWrite(output_kmer);
            output_kmer << min_kmer.get().str() << "\n";
            for (size_t i = 0; i < n; ++i) {
                if (i != 0) {
                    output_cnt << " ";
                }
                if (alive[i] && top_kmer[i].first == *min_kmer) {
                    cnt_vector[i / 2] += top_kmer[i].second;
                }
            }
        }
        for (size_t i = 0; i < cnt_vector.size(); ++i) {
            if (i != 0) {
                output_cnt << " ";
            }
            output_cnt << cnt_vector[i];
        }
        output_cnt << "\n";
        for (size_t i = 0; i < n; ++i) {
            if (alive[i] && top_kmer[i].first == *min_kmer) {
                alive[i] = ReadKmerWithCount(*infiles[i], top_kmer[i]);
            }
        }
    }
    output_kmer.close();
    output_cnt.close();

    for (size_t i = 0; i < n; ++i) {
        string name = files[i] + KMER_SORTED_EXTENSION;
        system(("rm -f " + name).c_str());
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        PrintUsageInfo();
        return 0;
    }
    string one_min = DEFAULT_ONE_MIN;
    string all_min = DEFAULT_ALL_MIN;
    string kmer_length = DEFAULT_KMER_LENGTH;
    int i = 1;
    while (i < argc && argv[i][0] == '-') {
        if (strncmp(argv[i], "-one-min", 8) == 0) {
            one_min = argv[i + 1];
        }
        if (strncmp(argv[i], "-all-min", 8) == 0) {
            all_min = argv[i + 1];
        }
        if (strncmp(argv[i], "-kl", 3) == 0) {
            kmer_length = argv[i + 1];
        }
        i += 2;
    }
    if (argc - i < 3 || (argc - i - 1) % 2 != 0) {
        PrintUsageInfo();
        return 0;
    }
    std::vector<string> input_files;
    for (; i + 1 < argc; ++i) {
        input_files.push_back(argv[i]);
    }
    for (const string& file : input_files) {
        CountKmersOne(file, one_min, kmer_length);
    }

    FilterKmersAll(input_files, atoi(all_min.c_str()), atoi(kmer_length.c_str()), argv[argc - 1]);
    return 0;
}
