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
                    io::FileReadStream& contigs_stream, std::ofstream& out) {
        io::SingleRead contig;
        while (!contigs_stream.eof()) {
            contigs_stream >> contig;
            contig_id id = GetId(contig);
            DEBUG("Analyzing contig " << id);

            if (contig.size() < min_length_bound) {
                DEBUG("Fragment is shorter than min_length_bound " << min_length_bound);
                break;
            }

            auto abundance_vec = abundance_counter(contig.GetSequenceString(), contig.name());

            if (abundance_vec) {
                DEBUG("Successfully estimated abundance of " << id);
                out << id;
                for (auto mpl : *abundance_vec)
                     out << "\t" << mpl;
                out << std::endl;
            } else {
                DEBUG("Failed to estimate abundance of " << id);
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
    std::string contigs_path, kmer_mult_fn, contigs_abundance_fn;

    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('c', contigs_path)
            >> Option('n', sample_cnt)
            >> Option('m', kmer_mult_fn)
            >> Option('o', contigs_abundance_fn)
            >> Option('l', min_length_bound, size_t(0));
    } catch(GetOptEx &ex) {
        std::cout << "Usage: contig_abundance_counter -k <K> -c <contigs path> "
                "-n <sample cnt> -m <kmer multiplicities path> "
                "-o <contigs abundance path> [-l <contig length bound> (default: 0)]"  << std::endl;
        exit(1);
    }

    //TmpFolderFixture fixture("tmp");
    create_console_logger();

    KmerProfileIndex::SetSampleCount(sample_cnt);
    ContigAbundanceCounter abundance_counter(k, std::unique_ptr<ClusterAnalyzer>(new TrivialClusterAnalyzer()), kmer_mult_fn);

    io::FileReadStream contigs_stream(contigs_path);

    std::ofstream out(contigs_abundance_fn);

    Runner::Run(abundance_counter, min_length_bound,
                contigs_stream, out);
    return 0;
}
