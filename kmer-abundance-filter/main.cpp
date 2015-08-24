#include <cstdio>
#include <string.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>

using std::string;
using std::vector;
using std::shared_ptr;
using std::make_shared;

const string DEFAULT_ONE_MIN = "3";
const string DEFAULT_ALL_MIN = "3";
const char DEFAULT_KMER_LENGTH[] = "31";
const char KMC_TMP_EXTENSION[] = ".kmc";

void PrintUsageInfo() {
    puts("Usage: kmer-abundance-filter [options] <file 1> <file 2> ... <file N> <output file>");
    puts("Options:");
    puts("-one-min - minimal kmer abundance in each file");
    puts("-all-min - minimal kmer abundance in all files");
    puts("-kl - kmer length");
}

void CountKmersOne(const string& filename, const string& abundance_min, const string& kmer_length) {
    system(("~/libs/kmc/kmc -k" + kmer_length + " -ci" + abundance_min + " " + filename + " " + filename + " .").c_str());
    system(("~/libs/kmc/kmc_dump " + filename + " " + filename + KMC_TMP_EXTENSION).c_str());
    system(("rm -f " + filename + ".kmc_pre " + filename + ".kmc_suf").c_str());
}

bool ReadKmerWithAbundance(std::ifstream& infile, std::pair<string, int>& res) {
    string kmer;
    int abundance;
    if (infile >> kmer >> abundance) {
        res = {kmer, abundance};
        return true;
    }
    return false;
}

void FilterKmersAll(const std::vector<string>& files, int all_min, const string& output_file) {
    size_t n = files.size();
    vector<shared_ptr<std::ifstream>> infiles;
    for (size_t i = 0; i < n; ++i) {
        string name = files[i] + KMC_TMP_EXTENSION;
        infiles.push_back(std::make_shared<std::ifstream>(name));
    }
    vector<std::pair<string, int>> top_kmer(n);
    vector<bool> alive(n);

    for (size_t i = 0; i < n; i++) {
        alive[i] = ReadKmerWithAbundance(*infiles[i], top_kmer[i]);
    }

    std::ofstream output_kmer(output_file + ".kmer");
    std::ofstream output_cnt(output_file + ".cnt", std::ios::binary);
    output_cnt.write((char*)&n, sizeof(size_t));
    while (true) {
        string min_kmer = "";
        int cnt_min = 0;
        for (size_t i = 0; i < n; ++i) {
            if (alive[i]) {
                string& cur_kmer = top_kmer[i].first;
                if (min_kmer == "" || cur_kmer < min_kmer) {
                    min_kmer = cur_kmer;
                    cnt_min = 0;
                }
                if (cur_kmer == min_kmer) {
                    cnt_min++;
                }
            }
        }
        if (min_kmer == "") {
            break;
        }
        if (cnt_min >= all_min) {
            output_kmer << min_kmer;
            for (size_t i = 0; i < n; ++i) {
                if (alive[i] && top_kmer[i].first == min_kmer) {
                    output_cnt.write((char*)&top_kmer[i].second, sizeof(int));
                    alive[i] = ReadKmerWithAbundance(*infiles[i], top_kmer[i]);
                } else {
                    output_cnt.write(0, sizeof(int));
                }
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
    if (argc - i < 2) {
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

    FilterKmersAll(input_files, atoi(all_min.c_str()), argv[argc - 1]);
    return 0;
}
