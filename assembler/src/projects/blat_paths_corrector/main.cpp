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
    std::string canu_contigs_file;
    std::string filtered_contigs_names_file;
    std::string saves_folder;
    std::string transcript_paths_file;
    std::string paths_save_file;
    std::string output_dir;
    std::string blat_output;
    unsigned int nthreads;
};


void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.saves_folder << value("saves folder"),
      cfg.transcript_paths_file << value("blat_filtered.paths.txt"),
      cfg.canu_contigs_file << value("canu contigs"),
      cfg.filtered_contigs_names_file << value("blat_filtered.names.txt"),
      cfg.blat_output << value("blat_output.psl"),
      (required("-k") & integer("int", cfg.k)) % "k-mer length to use",
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

std::vector<std::string> ReadContigsNames(std::string const & filtered_contigs_names_file) {
    std::ifstream inp(filtered_contigs_names_file);
    if (!inp.is_open())
        throw "Cannot open " + filtered_contigs_names_file;
    std::vector<std::string> names;
    std::string name;
    while (inp >> name) {
        if (name.empty())
            continue;
        names.push_back(std::move(name));
    }
    return names;
}

void WriteWithWidth(std::ostream & out, std::string const & seq, size_t width = 100) {
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
            auto edge_id = GetEdgeId(record.template Get<Columns::T_name>(), graph);
            output << ">mapped_onto_contig_from_" << record.template Get<Columns::Q_start>() << "_to_" << record.template Get<Columns::Q_end>() << ";_real_id_" << edge_id << ";\n";
            WriteWithWidth(output, graph.EdgeNucls(edge_id).str());
        }
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


        #ifdef USE_R_SCRIPT
            auto input_paths = ParseInputPaths(cfg.transcript_paths_file, graph);
            auto paths_names = ReadContigsNames(cfg.filtered_contigs_names_file);
            VERIFY(input_paths.size() == paths_names.size());
        #else
            ifstream blat_output_stream(cfg.blat_output);
            if (!blat_output_stream.is_open())
                throw "Cannot open " + cfg.blat_output;
            INFO("Started blat output reading");
            auto res = Read(blat_output_stream, GetFilter<WISHED_COLUMNS>());
            auto fragments = GetContigFragments(res, GetNonDropper<WISHED_COLUMNS>(), graph.k());
            MakeAllFilteredEdgesDump(res, fragments, graph, output_dir);
            auto paths_data = MakePaths(res, fragments, graph);
            auto const & input_paths = paths_data.first;
            auto const & paths_names = paths_data.second;
        #endif
        INFO("Scanned " << input_paths.size() << " paths");

        auto paths = Launch(gp, PathThreadingParams(), input_paths, contigs, paths_names, scaffolds, nthreads);

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
