#include <array>
#include <string>
#include <iostream>
#include "getopt_pp/getopt_pp.h"
#include "dev_support/verify.hpp"
#include "dev_support/logger/log_writers.hpp"
#include "dev_support/logger/logger.hpp"
#include "formats.hpp"
#include "pipeline/graph_pack.hpp"
#include "io/reads_io/file_reader.hpp"
#include "pipeline/graphio.hpp"

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

static const size_t CONTIG_SPLIT_LENGTH = 10000;
static const size_t MAX_SAMPLE_CNT = 20;
static const uint INVALID_MPL = uint(-1);

class ContigAbundanceCounter {
private:
    typedef uint Mpl;
    typedef typename std::array<Mpl, MAX_SAMPLE_CNT> MplVector;
    typedef typename std::array<double, MAX_SAMPLE_CNT> AbundanceVector;
    typedef typename InvertableStoring::immutant_inverter<MplVector> InverterT;

    unsigned k_;
    size_t sample_cnt_;
    size_t min_length_bound_;
    InverterT inverter_;

    KeyStoringMap<conj_graph_pack::seq_t,
        MplVector,
        kmer_index_traits<conj_graph_pack::seq_t>,
        InvertableStoring> kmer_mpl_;


    void FillMplMap(const std::string& kmers_mpl_file) {
        //INFO("Fill start");
        for (auto it = kmer_mpl_.value_begin(); it != kmer_mpl_.value_end(); ++it) {
            it->fill(INVALID_MPL);
        }
        std::ifstream kmers_in(kmers_mpl_file + ".kmer", std::ios::binary);
        std::ifstream kmers_mpl_in(kmers_mpl_file + ".mpl", std::ios::binary);
        //std::ifstream kmers_in(kmers_mpl_file + ".kmer");
        //std::ifstream kmers_mpl_in(kmers_mpl_file + ".mpl");
        while (true) {
//            std::string kmer_str;
//            kmers_in >> kmer_str;
            runtime_k::RtSeq kmer(k_);
            kmer.BinRead(kmers_in);
            if (kmers_in.fail()) {
                break;
            }

//            VERIFY(kmer_str.length() == k_);
//            conj_graph_pack::seq_t kmer(k_, kmer_str.c_str());
//            kmer = gp_.kmer_mapper.Substitute(kmer);

            MplVector mpls;
            mpls.fill(INVALID_MPL);
            for (size_t i = 0; i < sample_cnt_; ++i) {
//                kmers_mpl_in.read((char*) &kmer_mpl[i], sizeof(uint));
                kmers_mpl_in >> mpls[i];
                VERIFY(!kmers_mpl_in.fail());
            }

            auto kwh = kmer_mpl_.ConstructKWH(kmer);
            VERIFY(kmer_mpl_.valid(kwh));
            kmer_mpl_.put_value(kwh, mpls, inverter_);
        }
    }

    bool AddToAbundances(std::vector<std::vector<uint32_t>>& storage, const MplVector& kmer_mpls) const {
        bool invalid = (kmer_mpls[0] == INVALID_MPL);
        for (size_t i = 0; i < sample_cnt_; ++i) {
            if (invalid) {
                VERIFY(kmer_mpls[i] == INVALID_MPL);
            } else {
                VERIFY(kmer_mpls[i] != INVALID_MPL);
                storage[i].push_back(kmer_mpls[i]);
            }
        }
        return !invalid;
    };

    boost::optional<AbundanceVector> EstimateAbundance(const Sequence& seq) const {
        std::vector<std::vector<uint>> abundance_storage(sample_cnt_);

        size_t contributing_kmer_cnt = 0;

        VERIFY(seq.size() >= k_);
        auto kwh = kmer_mpl_.ConstructKWH(runtime_k::RtSeq(k_, seq));
        kwh >>= 'A';

        for (size_t j = k_ - 1; j < seq.size(); ++j) {
            kwh <<= seq[j];
            if (kmer_mpl_.valid(kwh)) {
                if (AddToAbundances(abundance_storage, kmer_mpl_.get_value(kwh, inverter_))) {
                    contributing_kmer_cnt++;
                }
            }
        }

        //FIXME magic constant
        if (contributing_kmer_cnt * 5 < seq.size() - k_ + 1) {
            return boost::none;
        }

        AbundanceVector contig_abundance;
        contig_abundance.fill(-1.);

        for (size_t i = 0; i < sample_cnt_; ++i) {
            VERIFY(abundance_storage[i].size() == contributing_kmer_cnt);
            //set contig abundance as mean across kmer multiplicities
            contig_abundance[i] = std::accumulate(abundance_storage[i].begin(), abundance_storage[i].end(), 0, std::plus<double>());
            contig_abundance[i] /= abundance_storage[i].size();
        }

        return boost::optional<AbundanceVector>(contig_abundance);
    }

public:
    ContigAbundanceCounter(unsigned k,
                           size_t sample_cnt,
                           size_t min_length_bound,
                           const std::string& work_dir) :
                               k_(k),
                               sample_cnt_(sample_cnt),
                               min_length_bound_(min_length_bound),
                               kmer_mpl_(k_, work_dir) {

    }

    void Init(const std::string& kmer_mpl_file, size_t read_buffer_size) {
        INFO("Filling kmer multiplicities. Sample cnt " << sample_cnt_);
        DeBruijnKMerKMerSplitter<StoringTypeFilter<InvertableStoring>>
                splitter(kmer_mpl_.workdir(), k_, k_, true, read_buffer_size);

        splitter.AddKMers(kmer_mpl_file + ".kmer");

        KMerDiskCounter<runtime_k::RtSeq> counter(kmer_mpl_.workdir(), splitter);

        kmer_mpl_.BuildIndex(counter, 16, /*nthreads*/ 1);

        FillMplMap(kmer_mpl_file);
    }

    void operator()(io::SingleStream& contigs,
                          const std::string& contigs_mpl_file) const {
        INFO("Calculating contig abundancies");
        std::ofstream id_out(contigs_mpl_file + ".id");
        std::ofstream mpl_out(contigs_mpl_file + ".mpl");

        io::SingleRead full_contig;
        size_t len = CONTIG_SPLIT_LENGTH; //TODO: adaptive?
        while (!contigs.eof()) {
            contigs >> full_contig;

            for (size_t i = 0; i < full_contig.size(); i += len) {
                if (full_contig.size() - i < min_length_bound_)
                    break;
                io::SingleRead contig = full_contig.Substr(i, std::min(i + len, full_contig.size()));
                contig_id id = GetId(contig);
                DEBUG("Processing contig " << id);

                auto abundance_vec = EstimateAbundance(contig.sequence());

                if (abundance_vec) {
                    id_out << id << std::endl;

                    std::string delim = "";
                    for (size_t i = 0; i < sample_cnt_; ++i) {
                        mpl_out << delim << (*abundance_vec)[i];
                        delim = " ";
                    }
                    mpl_out << std::endl;
                } else {
                    DEBUG("Failed to estimate abundance of contig " << id);
                }
            }
        }

        id_out.close();
        mpl_out.close();
    }
private:
    DECL_LOGGER("ContigAbundanceCounter");
};

}

int main(int argc, char** argv) {
    using namespace debruijn_graph;
    using namespace GetOpt;

    size_t k, sample_cnt, min_length_bound;
    std::string work_dir, contigs_path;
    std::string kmer_mult_fn, contigs_abundance_fn;

    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('w', work_dir)
            >> Option('c', contigs_path)
            >> Option('n', sample_cnt)
            >> Option('m', kmer_mult_fn)
            >> Option('o', contigs_abundance_fn)
            >> Option('l', min_length_bound, size_t(0))
        ;
    } catch(GetOptEx &ex) {
        std::cout << "Usage: contig_abundance_counter -k <K> -w <work_dir> -c <contigs path> "
                "-n <sample cnt> -m <kmer multiplicities path> "
                "-o <contigs abundance path> [-l <contig length bound> (default: 0)]"  << std::endl;
        exit(1);
    }

    //TmpFolderFixture fixture("tmp");
    create_console_logger();

    auto contigs_stream_ptr = make_shared<io::FileReadStream>(contigs_path);

    ContigAbundanceCounter abundance_counter(k, sample_cnt, min_length_bound, work_dir);
    abundance_counter.Init(kmer_mult_fn, /*fixme some buffer size*/0);
    abundance_counter(*contigs_stream_ptr, contigs_abundance_fn);

    return 0;
}
