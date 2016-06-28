#include <array>
#include <string>
#include <iostream>
#include "getopt_pp/getopt_pp.h"
#include "dev_support/verify.hpp"
#include "dev_support/logger/log_writers.hpp"
#include "dev_support/logger/logger.hpp"
#include "pipeline/graph_pack.hpp"
#include "io/reads_io/file_reader.hpp"
#include "pipeline/graphio.hpp"
#include "formats.hpp"
#include "contig_abundance.hpp"

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
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
            >> Option('l', min_length_bound, size_t(0));
    } catch(GetOptEx &ex) {
        std::cout << "Usage: contig_abundance_counter -k <K> -w <work_dir> -c <contigs path> "
                "-n <sample cnt> -m <kmer multiplicities path> "
                "-o <contigs abundance path> [-l <contig length bound> (default: 0)]"  << std::endl;
        exit(1);
    }

    //TmpFolderFixture fixture("tmp");
    create_console_logger();

    size_t split_length = 10000;

    ContigAbundanceCounter abundance_counter(k, sample_cnt, SingleClusterAnalyzer(sample_cnt), work_dir);
    abundance_counter.Init(kmer_mult_fn);

    std::ofstream id_out(contigs_abundance_fn + ".id");
    std::ofstream mpl_out(contigs_abundance_fn + ".mpl");

    io::FileReadStream contigs_stream(contigs_path);

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
            contig_id id = GetId(contig);
            DEBUG("Processing fragment # " << (i / split_length) << " with id " << id);

            auto abundance_vec = abundance_counter(contig.GetSequenceString(), contig.name());

            if (abundance_vec) {
                stringstream ss;
                copy(abundance_vec->begin(), abundance_vec->begin() + sample_cnt,
                     ostream_iterator<Mpl>(ss, " "));
                DEBUG("Successfully estimated abundance of " << id << " : " << ss.str());

                id_out << id << std::endl;
                mpl_out << ss.str() << std::endl;
            } else {
                DEBUG("Failed to estimate abundance of " << id);
            }
        }
    }
    return 0;
}
