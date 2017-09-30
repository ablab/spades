#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <algorithm>
#include <libcxx/sort.hpp>
#include <boost/optional/optional.hpp>
#include "getopt_pp/getopt_pp.h"
#include "kmc_api/kmc_file.h"
//#include "omp.h"
#include "io/kmers/mmapped_reader.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/stl_utils.hpp"
#include "utils/ph_map/perfect_hash_map_builder.hpp"
#include "utils/kmer_mph/kmer_splitters.hpp"
#include "logger.hpp"

using std::string;
using std::vector;

const string KMER_PARSED_EXTENSION = ".bin";
const string KMER_SORTED_EXTENSION = ".sorted";

class KmerMultiplicityCounter {

    size_t k_, sample_cnt_;
    std::string file_prefix_;

    //TODO: get rid of intermediate .bin file
    string ParseKmc(const string& filename) {
        CKMCFile kmcFile;
        kmcFile.OpenForListing(filename);
        CKmerAPI kmer((unsigned int) k_);
        uint32 count;
        std::string parsed_filename = filename + KMER_PARSED_EXTENSION;
        std::ofstream output(parsed_filename, std::ios::binary);
        while (kmcFile.ReadNextKmer(kmer, count)) {
            RtSeq seq(k_, kmer.to_string());
            seq.BinWrite(output);
            seq_element_type tmp = count;
            output.write((char*) &(tmp), sizeof(seq_element_type));
        }
        output.close();
        return parsed_filename;
    }

    string SortKmersCountFile(const string& filename) {
        MMappedRecordArrayReader<seq_element_type> ins(filename, RtSeq::GetDataSize(k_) + 1, false);
        libcxx::sort(ins.begin(), ins.end(), adt::array_less<seq_element_type>());
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

    fs::TmpFile FilterCombinedKmers(fs::TmpDir workdir, const std::vector<string>& files, size_t all_min) {
        size_t n = files.size();
        vector<std::unique_ptr<ifstream>> infiles;
        infiles.reserve(n);
        for (auto fn : files) {
            INFO("Processing " << fn);
            auto parsed = ParseKmc(fn);
            auto sorted = SortKmersCountFile(parsed);
            infiles.emplace_back(new std::ifstream(sorted));
        }
        vector<std::pair<RtSeq, uint32>> top_kmer(n, {RtSeq(k_), 0});
        vector<bool> alive(n, false);

        for (size_t i = 0; i < n; i++) {
            alive[i] = ReadKmerWithCount(*infiles[i], top_kmer[i]);
        }

        auto kmer_file = fs::tmp::make_temp_file("kmer", workdir);

        typedef uint16_t Mpl;
        std::ofstream output_kmer(*kmer_file, std::ios::binary);
        std::ofstream mpl_file(file_prefix_ + ".bpr", std::ios_base::binary);

        RtSeq::less3 kmer_less;
        while (true) {
            boost::optional<RtSeq> min_kmer;
            size_t cnt_min = 0;
            for (size_t i = 0; i < n; ++i) {
                if (alive[i]) {
                    RtSeq& cur_kmer = top_kmer[i].first;
                    if (!min_kmer || kmer_less(cur_kmer, *min_kmer)) {
                        min_kmer = cur_kmer;
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
                for (size_t i = 0; i < n; ++i) {
                    if (alive[i] && top_kmer[i].first == *min_kmer) {
                        cnt_vector[i] += top_kmer[i].second;
                    }
                }

                for (size_t mpl : cnt_vector) {
                    mpl_file.write(reinterpret_cast<char *>(&mpl), sizeof(Mpl));
                }
            }
            for (size_t i = 0; i < n; ++i) {
                if (alive[i] && top_kmer[i].first == *min_kmer) {
                    alive[i] = ReadKmerWithCount(*infiles[i], top_kmer[i]);
                }
            }
        }
        return kmer_file;
    }

    void BuildKmerIndex(fs::TmpDir workdir, fs::TmpFile kmer_file, size_t sample_cnt, size_t nthreads) {
        INFO("Initializing kmer profile index");

        typedef size_t Offset;
        using namespace utils;

        static const size_t read_buffer_size = 0; //FIXME some buffer size
        DeBruijnKMerKMerSplitter<StoringTypeFilter<InvertableStoring>>
            splitter(workdir, k_, k_, true, read_buffer_size);
        splitter.AddKMers(*kmer_file);

        KMerDiskCounter<RtSeq> counter(workdir, splitter);
        KeyStoringMap<RtSeq, Offset, kmer_index_traits<RtSeq>, InvertableStoring> kmer_mpl(k_);
        BuildIndex(kmer_mpl, counter, 16, nthreads);
        INFO("Built index with " << kmer_mpl.size() << " kmers");

        //Building kmer->profile offset index
        std::ifstream kmers_in(*kmer_file, std::ios::binary);
        InvertableStoring::trivial_inverter<Offset> inverter;
        RtSeq kmer(k_);
        for (Offset offset = 0; ; offset += sample_cnt) {
            kmer.BinRead(kmers_in);
            if (kmers_in.fail()) {
                break;
            }

//            VERIFY(kmer_str.length() == k_);
//            conj_graph_pack::seq_t kmer(k_, kmer_str.c_str());
//            kmer = gp_.kmer_mapper.Substitute(kmer);

            auto kwh = kmer_mpl.ConstructKWH(kmer);
            VERIFY(kmer_mpl.valid(kwh));
            kmer_mpl.put_value(kwh, offset, inverter);
        }

        std::ofstream map_file(file_prefix_ + ".kmm", std::ios_base::binary | std::ios_base::out);
        kmer_mpl.BinWrite(map_file);
        INFO("Saved kmer profile map");
    }

public:
    KmerMultiplicityCounter(size_t k, std::string file_prefix):
        k_(k), file_prefix_(std::move(file_prefix)) {
    }

    void CombineMultiplicities(const vector<string>& input_files, size_t min_samples, const string& tmpdir, size_t nthreads = 1) {
        auto workdir = fs::tmp::make_temp_dir(tmpdir, "kmidx");
        auto kmer_file = FilterCombinedKmers(workdir, input_files, min_samples);
        BuildKmerIndex(workdir, kmer_file, input_files.size(), nthreads);
    }
private:
    DECL_LOGGER("KmerMultiplicityCounter");
};

void PrintUsageInfo() {
    std::cout << "Usage: kmer_multiplicity_counter [options] -f files_dir" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "-k - kmer length" << std::endl;
    std::cout << "-n - sample count" << std::endl;
    std::cout << "-o - output file prefix" << std::endl;
    std::cout << "-t - number of threads (default: 1)" << std::endl;
    std::cout << "-s - minimal number of samples to contain kmer" << std::endl;
    std::cout << "files_dir must contain two files (.kmc_pre and .kmc_suf) with kmer multiplicities for each sample from 1 to n" << std::endl;
}

int main(int argc, char *argv[]) {
    using namespace GetOpt;
    create_console_logger();

    size_t k, sample_cnt, min_samples, nthreads;
    string output, work_dir;

    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('n', sample_cnt)
            >> Option('s', min_samples)
            >> Option('o', output)
            >> Option('t', "threads", nthreads, size_t(1))
            >> Option('f', work_dir)
        ;
    } catch(GetOptEx &ex) {
        PrintUsageInfo();
        exit(1);
    }

    std::vector<string> input_files;
    for (size_t i = 1; i <= sample_cnt; ++i) {
        input_files.push_back(work_dir + "/sample" + std::to_string(i));
    }

    KmerMultiplicityCounter kmcounter(k, output);
    kmcounter.CombineMultiplicities(input_files, min_samples, work_dir, nthreads);
    return 0;
}
