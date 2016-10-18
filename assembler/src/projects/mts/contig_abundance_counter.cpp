#include <array>
#include <string>
#include <iostream>
#include "getopt_pp/getopt_pp.h"
#include "io/reads/file_reader.hpp"
#include "io/reads/osequencestream.hpp"
#include "pipeline/graphio.hpp"
#include "logger.hpp"
#include "formats.hpp"
#include "contig_abundance.hpp"

using namespace debruijn_graph;

//Helper class to have scoped DEBUG()
class Runner {
public:
    static void Run(ContigAbundanceCounter& abundance_counter, size_t min_length_bound,
                    io::FileReadStream& contigs_stream, io::osequencestream& splits_os,
                    std::ofstream& id_out, std::ofstream& mpl_out) {
        static const size_t split_length = 10000;
        io::SingleRead full_contig;
        while (!contigs_stream.eof()) {
            contigs_stream >> full_contig;
            DEBUG("Analyzing contig " << GetId(full_contig));

            for (size_t i = 0; i < full_contig.size(); i += split_length) {
                if (full_contig.size() - i < min_length_bound) {
                    DEBUG("Fragment shorter than min_length_bound " << min_length_bound);
                    break;
                }

                io::SingleRead contig = full_contig.Substr(i, std::min(i + split_length, full_contig.size()));
                splits_os << contig;

                contig_id id = GetId(contig);
                DEBUG("Processing fragment # " << (i / split_length) << " with id " << id);

                auto abundance_vec = abundance_counter(contig.GetSequenceString(), contig.name());

                if (abundance_vec) {
                    stringstream ss;
                    copy(abundance_vec->begin(), abundance_vec->end(),
                         ostream_iterator<Mpl>(ss, " "));
                    DEBUG("Successfully estimated abundance of " << id << " : " << ss.str());

                    id_out << id << std::endl;
                    mpl_out << ss.str() << std::endl;
                } else {
                    DEBUG("Failed to estimate abundance of " << id);
                }
            }
        }
    }
private:
    DECL_LOGGER("ContigAbundanceCounter");
};

int main(int argc, char** argv) {
    using namespace GetOpt;

    unsigned k;
    size_t sample_cnt, min_length_bound;
    std::string work_dir, contigs_path, splits_path;
    std::string kmer_mult_fn, contigs_abundance_fn;

    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('w', work_dir)
            >> Option('c', contigs_path)
            >> Option('f', splits_path)
            >> Option('n', sample_cnt)
            >> Option('m', kmer_mult_fn)
            >> Option('o', contigs_abundance_fn)
            >> Option('l', min_length_bound, size_t(0));
    } catch(GetOptEx &ex) {
        std::cout << "Usage: contig_abundance_counter -k <K> -w <work_dir> -c <contigs path> "
                "-n <sample cnt> -m <kmer multiplicities path> -f <splits_path> "
                "-o <contigs abundance path> [-l <contig length bound> (default: 0)]"  << std::endl;
        exit(1);
    }

    //TmpFolderFixture fixture("tmp");
    create_console_logger();

    SetSampleCount(sample_cnt);
    ContigAbundanceCounter abundance_counter(k, SingleClusterAnalyzer(), work_dir);
    abundance_counter.Init(kmer_mult_fn);

    io::FileReadStream contigs_stream(contigs_path);
    io::osequencestream splits_os(splits_path);

    std::ofstream id_out(contigs_abundance_fn + ".id");
    std::ofstream mpl_out(contigs_abundance_fn + ".mpl");

    Runner::Run(abundance_counter, min_length_bound,
                contigs_stream, splits_os,
                id_out, mpl_out);
    return 0;
}
