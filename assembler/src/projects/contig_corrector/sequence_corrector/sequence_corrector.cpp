//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sequence_corrector.hpp"
#include "helpers/common.hpp"
#include "helpers/minimap_output_reader.hpp"
#include "helpers/blat_output_reader.hpp"
#include "helpers/aligner_output_postprocessing.hpp"

#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"
#include "common/assembly_graph/graph_support/coverage_uniformity_analyzer.hpp"
#include "common/assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "common/toolchain/utils.hpp"
#include "common/utils/memory_limit.hpp"

#include "utils/parallel/openmp_wrapper.h"
#include "utils/parallel/parallel_wrapper.hpp"
#include "utils/filesystem/temporary.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <algorithm>

namespace sequence_corrector {

using namespace std;
using namespace path_extend;
using namespace helpers;

// #define WISHED_COLUMNS BColumns, BColumns::match, BColumns::strand, BColumns::block_count,\
//           BColumns::Q_name, BColumns::Q_size, BColumns::Q_start, BColumns::Q_end,\
//           BColumns::T_name, BColumns::T_size, BColumns::T_start, BColumns::T_end

#define WISHED_COLUMNS MColumns, MColumns::match, MColumns::strand,\
          MColumns::Q_name, MColumns::Q_size, MColumns::Q_start, MColumns::Q_end,\
          MColumns::T_name, MColumns::T_size, MColumns::T_start, MColumns::T_end


struct gcfg {
    gcfg()
        : k(-1)
        , output_dir("output")
        , nthreads(omp_get_max_threads() / 2 + 1)
    {}

    size_t k;
    size_t drop_alg;
    std::string canu_contigs_file;
    std::string saves_folder;
    std::string paths_save_file;
    std::string output_dir;
    std::string aligner_output;
    unsigned int nthreads;
} cfg;


clipp::group GetCLI() {
  using namespace clipp;

  auto cli = (
      cfg.saves_folder << value("saves folder"),
      cfg.canu_contigs_file << value("canu contigs"),
      cfg.aligner_output << value("aligner_output"),
      (required("-k") & integer("int", cfg.k)) % "k-mer length to use",
      (required("-d") & integer("int", cfg.drop_alg)) % "0 = drop nothing, 1 = full drop, 2 = transitive drop by IDY",
      (option("-t") & integer("int", cfg.nthreads)) % "# of threads to use",
      (option("-o") & value("dir", cfg.output_dir)) % "output directory",
      (option("-s") & value("file", cfg.paths_save_file)) % "save_path_to_scaffolds"
  );

  return cli;
}

namespace {

void ReadScaffolds(PathContainer& scaffolds, Graph const & graph, std::string const & paths_save_file) {
    std::ifstream inp(paths_save_file);
    if (!inp.is_open())
        throw "Cannot open " + paths_save_file;
    size_t amount_of_paths;
    inp >> amount_of_paths;
    for (size_t i = 0; i < amount_of_paths; ++i) {
        SimpleBidirectionalPath gapped_path;
        gapped_path.BinRead(inp);
        auto path = BidirectionalPath::create(graph, std::move(gapped_path));
        auto conj_path = BidirectionalPath::clone_conjugate(path);
        scaffolds.AddPair(move(path), move(conj_path));
    }
}

void WriteWithWidth(std::ostream & out, std::string const & seq, size_t width = 50) {
    for (size_t pos = 0; pos < seq.size(); pos += width)
        out << seq.substr(pos, width) << '\n';
}

template<class Columns, Columns ... columns>
void MakeAllFilteredEdgesDump(Records<Columns, columns ...> const & records, MapFromContigNameToContigFragments const & contig_fragments, Graph const & graph, std::string const & output_dir) {
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

ScaffoldingUniqueEdgeStorage GetUniqueEdgeStorage(GraphPack const & gp) {
    auto const & graph = gp.get<Graph>();
    constexpr auto min_unique_length = 2000;
    constexpr auto uniformity_fraction_threshold = 0.8;
    constexpr auto nonuniform_coverage_variation = 50;
    auto unique_variation = 0.5;

    CoverageUniformityAnalyzer coverage_analyzer(graph, min_unique_length);
    double median_coverage = coverage_analyzer.CountMedianCoverage();
    double uniformity_fraction = coverage_analyzer.UniformityFraction(unique_variation, median_coverage);
    INFO ("median coverage for edges longer than " << min_unique_length << " is " << median_coverage <<
        " uniformity " << size_t(uniformity_fraction * 100) << "%");

    bool uniform_coverage = false;
    if (math::gr(uniformity_fraction, uniformity_fraction_threshold)) {
        uniform_coverage = true;
    }
    if (!uniform_coverage) {
        unique_variation = nonuniform_coverage_variation;
        INFO("Coverage is not uniform, we do not rely on coverage for long edge uniqueness");
    }

    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp, min_unique_length, unique_variation);
    ScaffoldingUniqueEdgeStorage unique_edge_storage;
    unique_edge_analyzer.FillUniqueEdgeStorage(unique_edge_storage);
    return unique_edge_storage;
}

template<class Columns, Columns ... columns>
void DropAnother(Records<Columns, columns ...> const & records, MapFromContigNameToContigFragments & contig_fragments) {
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

template<class Columns, Columns ... columns>
DropAlg<Columns, columns ...> GetDropAlg(size_t num) {
    switch (num) {
        case 0: return GetNonDropper<Columns, columns ...>(); break;
        case 1: return GetFullDropper<Columns, columns ...>(); break;
        case 2: return GetTransitiveDropperByIDY<Columns, columns ...>(); break;
        default: throw std::string("unknown the value of -d option"); break;
    }
}

PathContainer GetUnique(PathContainer const & paths) {
    PathContainer ans;
    vector<unique_ptr<BidirectionalPath>> forward_paths;
    for (auto const & path : paths)
        forward_paths.push_back(BidirectionalPath::clone_conjugate(path.second));
    std::sort(forward_paths.begin(), forward_paths.end(), [](unique_ptr<BidirectionalPath> const & lhs, unique_ptr<BidirectionalPath> const & rhs) {
        if (lhs->Size() != rhs->Size())
            return lhs->Size() < rhs->Size();
        if (lhs->Length() != rhs->Length())
            return lhs->Length() < rhs->Length();

        for (size_t i = 0; i < lhs->Size(); ++i) {
            if (lhs->At(i) != rhs->At(i))
                return lhs->At(i) < rhs->At(i);
        }
        return false;
    });

    forward_paths.erase(std::unique(forward_paths.begin(), forward_paths.end(), [](auto const & lhs, auto const & rhs) { return *lhs == *rhs; }), forward_paths.end());

    for (auto & path : forward_paths) {
        auto cpath = BidirectionalPath::clone_conjugate(path);
        ans.AddPair(move(path), move(cpath));
    }
    return ans;
}

constexpr char BASE_NAME[] = "graph_pack";

} //namespace

int main() {
    const size_t GB = 1 << 30;

    auto nthreads = cfg.nthreads;
    auto k = cfg.k;
    std::string &output_dir = cfg.output_dir;

    START_BANNER("SPAdes standalone contig corrector");
    utils::limit_memory(30 * GB);

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


    ifstream aligner_output_stream(cfg.aligner_output);
    if (!aligner_output_stream.is_open())
        throw "Cannot open " + cfg.aligner_output;

    auto unique_edge_storage = GetUniqueEdgeStorage(gp);
    auto unique_edge_checker = [&unique_edge_storage](EdgeId id) { return unique_edge_storage.IsUnique(id); };
    auto filter = GetFilter<WISHED_COLUMNS>(std::move(unique_edge_checker), graph);

    INFO("Started aligner output reading");
    auto res = ReadMinimapOutput(aligner_output_stream, filter);
    auto fragments = GetContigFragments(res, GetDropAlg<WISHED_COLUMNS>(cfg.drop_alg), graph.k());
    // DropAnother(res, fragments);
    MakeAllFilteredEdgesDump(res, fragments, graph, output_dir);
    auto input_paths = MakePaths(res, fragments, graph);
    INFO("Scanned " << input_paths.size() << " paths");

    auto paths = Launch(gp, PathThreadingParams(), input_paths, contigs, scaffolds, nthreads);

    ContigWriter writer(graph, MakeContigNameGenerator(config::pipeline_type::base, gp));
    writer.OutputPaths(GetUnique(paths), fs::append_path(output_dir, "connected_paths.fasta"));

    ofstream contigs_output(fs::append_path(output_dir, "remapped_paths.fasta"));
    VERIFY(contigs_output.is_open());

    for (auto const & contig : contigs) {
        if (contig.seq.empty())
            continue;
        contigs_output << '>' << contig.name << " len=" << contig.seq.size() << '\n';
        WriteWithWidth(contigs_output, contig.seq);
    }

    INFO("SPAdes standalone contig corrector finished");

    return 0;
}

} // namespace sequence_corrector