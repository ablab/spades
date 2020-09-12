//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "common.hpp"
#include "blat_output_reader.hpp"
#include "blat_output_postprocessing.hpp"

#include "common/assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "common/toolchain/utils.hpp"
#include "common/utils/memory_limit.hpp"

#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/parallel/parallel_wrapper.hpp"
#include "utils/filesystem/temporary.hpp"

#include <clipp/clipp.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <algorithm>

using namespace std;
using namespace path_extend;

#define WISHED_COLUMNS Columns::match, Columns::strand, Columns::block_count,\
          Columns::Q_name, Columns::Q_size, Columns::Q_start, Columns::Q_end,\
          Columns::T_name, Columns::T_size, Columns::T_start, Columns::T_end

struct gcfg {
    gcfg()
        : k(21)
        , output_dir("output")
        , nthreads(omp_get_max_threads() / 2 + 1)
    {}

    size_t k;
    size_t drop_alg;
    std::string canu_contigs_file;
    std::string saves_folder;
    std::string paths_save_file;
    std::string output_dir;
    std::string blat_output;
    unsigned int nthreads;
};


void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.saves_folder << value("saves folder"),
      cfg.canu_contigs_file << value("canu contigs"),
      cfg.blat_output << value("blat_output.psl"),
      (required("-k") & integer("int", cfg.k)) % "k-mer length to use",
      (required("-d") & integer("int", cfg.drop_alg)) % "0 = drop nothing, 1 = full drop, 2 = transitive drop by IDY",
      (option("-t") & integer("int", cfg.nthreads)) % "# of threads to use",
      (option("-o") & value("dir", cfg.output_dir)) % "output directory",
      (option("-s") & value("file", cfg.paths_save_file)) % "save_path_to_scaffolds"
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }
}

void ReadScaffolds(PathContainer& scaffolds, Graph const & graph, std::string const & paths_save_file) {
    std::ifstream inp(paths_save_file);
    if (!inp.is_open())
        throw "Cannot open " + paths_save_file;
    size_t amount_of_paths;
    inp >> amount_of_paths;
    for (size_t i = 0; i < amount_of_paths; ++i) {
        GappedPath gapped_path;
        gapped_path.BinRead(inp);
        auto path = make_unique<BidirectionalPath>(graph, std::move(gapped_path));
        auto conj_path = make_unique<BidirectionalPath>(graph, std::move(path->Conjugate()));
        scaffolds.AddPair(path.release(), conj_path.release());
    }
}

void WriteWithWidth(std::ostream & out, std::string const & seq, size_t width = 50) {
    for (size_t pos = 0; pos < seq.size(); pos += width)
        out << seq.substr(pos, width) << '\n';
}

template<Columns ... columns>
void MakeAllFilteredEdgesDump(Records<columns ...> const & records, MapFromContigNameToContigFragments const & contig_fragments, Graph const & graph, std::string const & output_dir) {
    auto dir = fs::append_path(output_dir, "filtered_edges_dump");
    fs::remove_if_exists(dir);
    fs::make_dir(dir);

    for (auto const & contig : contig_fragments) {
        if (contig.second.empty())
            continue;
        std::ofstream output(fs::append_path(dir, contig.first + ".fasta"));
        for (auto index : contig.second) {
            auto const & record = records[index];
            auto start_pos = record.template Get<Columns::Q_start>();
            auto end_pos = record.template Get<Columns::Q_end>();
            auto edge_id = GetEdgeId(record.template Get<Columns::T_name>(), graph);
            output << ">mapped_onto_contig_from_" << start_pos << "_to_" << end_pos << ";_real_id_" << edge_id << ";\n";
            WriteWithWidth(output, graph.EdgeNucls(edge_id).str());
        }
    }
}

template<Columns ... columns>
void DropAnother(Records<columns ...> const & records, MapFromContigNameToContigFragments & contig_fragments) {
    auto intersected =
    [](long long start, long long end) {
        auto dt = 100;
        auto range = std::make_pair(1471516-dt, 1471563+dt);
        return range.first <= start && start <= range.second || 
            range.first <= end && end <= range.second || 
            start < range.first && range.second < start;
    };

    for (auto & contig : contig_fragments) {
        vector<size_t> good_indexes;
        for (auto index : contig.second) {
            auto const & record = records[index];
            auto start_pos = record.template Get<Columns::Q_start>();
            auto end_pos = record.template Get<Columns::Q_end>() + 1;
            if (intersected(start_pos, end_pos))
                good_indexes.push_back(index);
        }
        contig.second = std::move(good_indexes);
    }
}

template<Columns ... columns>
DropAlg<columns ...> GetDropAlg(size_t num) {
    switch (num) {
        case 0: return GetNonDropper<columns ...>(); break;
        case 1: return GetFullDropper<columns ...>(); break;
        case 2: return GetTransitiveDropperByIDY<columns ...>(); break;
        default: throw std::string("unknown the value of -d option"); break;
    }
}
constexpr char BASE_NAME[] = "graph_pack";

int main(int argc, char* argv[]) {
    utils::segfault_handler sh;
    const size_t GB = 1 << 30;
    toolchain::create_console_logger(logging::L_INFO);

    gcfg cfg;

    srand(42);
    srandom(42);

    process_cmdline(argc, argv, cfg);
    try {
        auto nthreads = cfg.nthreads;
        auto k = cfg.k;
        std::string &output_dir = cfg.output_dir;

        START_BANNER("SPAdes standalone Blat paths corrector");
        utils::limit_memory(15 * GB);

        CHECK_FATAL_ERROR(runtime_k::MIN_K <= k, "k-mer size " << k << " is too low");
        CHECK_FATAL_ERROR(k < runtime_k::MAX_K, "k-mer size " << k << " is too high, recompile with larger SPADES_MAX_K option");
        CHECK_FATAL_ERROR(k & 1, "k-mer size must be odd");

        nthreads = spades_set_omp_threads(nthreads);
        INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << nthreads);

        fs::make_dir(output_dir);

        debruijn_graph::GraphPack gp(k, output_dir, 0);
        auto p = fs::append_path(cfg.saves_folder, BASE_NAME);
        io::binary::FullPackIO().Load(p, gp);

        auto contigs = ReadContigs(cfg.canu_contigs_file);
        #ifdef GOOD_NAME
        {
            ofstream contigs_output(fs::append_path(output_dir, "canu_contig.fasta"));
            VERIFY(contigs_output.is_open());
            for (auto & contig : contigs) {
                if (contig.name == GOOD_NAME) {
                    contigs_output << '>' << contig.name << " len=" << contig.seq.size() << '\n';
                    const auto d_pos = 100;
                    for (size_t pos = 0; pos < contig.seq.size(); pos += d_pos)
                        contigs_output << contig.seq.substr(pos, d_pos) << '\n';
                    break;
                }
            }
        }
        #endif

        auto& graph = gp.get<debruijn_graph::Graph>();
        PathContainer scaffolds;
        if (!cfg.paths_save_file.empty())
            ReadScaffolds(scaffolds, graph, cfg.paths_save_file);

        ifstream blat_output_stream(cfg.blat_output);
        if (!blat_output_stream.is_open())
            throw "Cannot open " + cfg.blat_output;
        INFO("Started blat output reading");
        auto res = Read(blat_output_stream, GetFilter<WISHED_COLUMNS>());
        auto fragments = GetContigFragments(res, GetDropAlg<WISHED_COLUMNS>(cfg.drop_alg), graph.k());
        // DropAnother(res, fragments);
        MakeAllFilteredEdgesDump(res, fragments, graph, output_dir);
        auto input_paths = MakePaths(res, fragments, graph);
        INFO("Scanned " << input_paths.size() << " paths");

        auto paths = Launch(gp, PathThreadingParams(), input_paths, contigs, scaffolds, nthreads);

        ContigWriter writer(graph, MakeContigNameGenerator(config::pipeline_type::base, gp));
        writer.OutputPaths(paths, fs::append_path(output_dir, "connected_paths.fasta"));
        
        ofstream contigs_output(fs::append_path(output_dir, "remapped_paths.fasta"));
        VERIFY(contigs_output.is_open());

        for (auto const & contig : contigs) {
            if (contig.seq.empty())
                continue;
            contigs_output << '>' << contig.name << " len=" << contig.seq.size() << '\n';
            WriteWithWidth(contigs_output, contig.seq);
        }
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    } catch (...) {
        std::cerr << "unknown exception caught" << std::endl;
        return EINTR;
    }

    INFO("SPAdes standalone blat paths corrector finished");

    return 0;
}
