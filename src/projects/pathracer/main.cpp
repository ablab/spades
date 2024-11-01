//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "fees.hpp"
#include "find_best_path.hpp"
#include "debruijn_graph_cursor.hpp"
#include "cursor_neighborhood.hpp"
#include "cursor_conn_comps.hpp"
#include "io/binary/binary.hpp"
#include "path_utils.hpp"
#include "cached_cursor.hpp"
#include "superpath_index.hpp"
#include "hmm_path_info.hpp"
#include "fasta_reader.hpp"

#include "stack_limit.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

#include "sequence/aa.hpp"

#include "visualization/visualization.hpp"
#include "io/graph/gfa_reader.hpp"
#include "io/reads/io_helper.hpp"
#include "io/reads/osequencestream.hpp"
#include "io/reads/file_reader.hpp"

#include "io/binary/graph.hpp"

#include "hmm/hmmfile.hpp"
#include "hmm/hmmmatcher.hpp"

#include "utils/parallel/openmp_wrapper.h"
#include "utils/segfault_handler.hpp"
#include "utils/memory_limit.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/verify.hpp"

#include "version.hpp"

#include <clipp/clipp.h>
#include <debug_assert/debug_assert.hpp>

#include <filesystem>
#include <string>
#include <functional>

#include <type_traits>

#include <unistd.h>  // getpid()
#include <sys/types.h>
#include <sys/stat.h>

extern "C" {
    #include "easel.h"
    #include "esl_sqio.h"
    #include "esl_vectorops.h"
}

struct main_assert : debug_assert::default_handler,
                     debug_assert::set_level<1> {};

enum class Mode {
    hmm,
    nucl,
    aa
};

enum class SeedMode {
    edges,
    scaffolds,
    scaffolds_one_by_one,
    edges_one_by_one,
    edges_scaffolds,
    exhaustive
};

struct PathracerConfig {
    std::string load_from = "";
    std::string hmmfile = "";
    std::filesystem::path output_dir = "";
    enum Mode mode = Mode::hmm;
    enum SeedMode seed_mode = SeedMode::edges_scaffolds;
    int threads = 16;
    size_t top = 10000;
    std::vector<std::string> edges;
    std::vector<std::string> queries;
    size_t max_size = size_t(-1);
    bool debug = false;
    bool draw = false;
    bool rescore = false;
    bool annotate_graph = true;
    double expand_coef = 2.;
    int expand_const = 20;
    size_t state_limits_coef = 1;
    bool local = false;
    bool parallel_component_processing = false;
    bool disable_depth_filter = false;
    size_t memory = 100;  // 100GB
    int use_experimental_i_loop_processing = true;
    std::string known_sequences = "";
    bool export_event_graph = false;
    double minimal_match_length = 0.9;
    size_t max_insertion_length = 30;
    double indel_rate = 0;
    bool quiet = false;
    bool ignore_names = false;

    hmmer::hmmer_cfg hcfg;
};

extern "C" {
#include "p7_config.h"

#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"

#include "hmmer.h"
}

void process_cmdline(int argc, char **argv, PathracerConfig &cfg) {
  using namespace clipp;

  std::string output_dir;
  auto cli = (
      cfg.hmmfile    << value("input file (HMM or sequence)"),
      cfg.load_from  << value("load from"),
      required("--output", "-o") & value("output directory", output_dir)    % "output directory",
      (option("--global").set(cfg.local, false) % "perform global-local (aka glocal) HMM matching [default]") |
      (cfg.local << option("--local") % "perform local-local HMM matching"),
      (option("--length", "-l") & value("value", cfg.minimal_match_length)) % "minimal length of resultant matched sequence; if <=1 then to be multiplied on the length of the pHMM [default: 0.9]",
      (option("--indel-rate", "-r") & value("value", cfg.indel_rate)) % "expected rate of nucleotides indels in graph edges [default: 0]. Used for AA pHMM alignment with frame shifts",
      (option("--top") & integer("N", cfg.top)) % "extract top N paths [default: 10000]",
      (option("--threads", "-t") & integer("NTHREADS", cfg.threads)) % "the number of parallel threads [default: 16]",
      (option("--memory", "-m") & integer("MEMORY", cfg.memory)) % "RAM limit for PathRacer in GB (terminates if exceeded) [default: 100]",
      (option("--max-size") & integer("SIZE", cfg.max_size)) % "maximal component size to consider [default: INF]",
      (cfg.quiet << option("--quiet", "-q")) % "be quiet, do not output anything to the console",
      (cfg.ignore_names << option("--ignore-names")) % "ignore HMM/sequence names (in case of dups, etc)",
      "Query type:" %
      one_of(option("--hmm").set(cfg.mode, Mode::hmm) % "match against HMM(s) [default]",
             option("--nt").set(cfg.mode, Mode::nucl) % "match against nucleotide string(s)",
             option("--aa").set(cfg.mode, Mode::aa) % "match agains amino acid string(s)"),
      (option("--queries") & values("queries", cfg.queries)) % "queries names to lookup [default: all queries from input query file]",
      (option("--edges") & values("edges", cfg.edges)) % "match around particular edges [default: all graph edges]",
      "Seeding options:" % (
          "Seeding mode:" %
          one_of(option("--seed-edges").set(cfg.seed_mode, SeedMode::edges) % "use graph edges as seeds",
                 option("--seed-scaffolds").set(cfg.seed_mode, SeedMode::scaffolds) % "use scaffolds paths as seeds",
                 option("--seed-edges-scaffolds").set(cfg.seed_mode, SeedMode::edges_scaffolds) % "use edges AND scaffolds paths as seeds [default]",
                 option("--seed-exhaustive").set(cfg.seed_mode, SeedMode::exhaustive) % "exhaustive mode, use ALL edges",
                 option("--seed-edges-1-by-1").set(cfg.seed_mode, SeedMode::edges_one_by_one) % "use edges as seeds (1 by 1)",
                 option("--seed-scaffolds-1-by-1").set(cfg.seed_mode, SeedMode::scaffolds_one_by_one) % "use scaffolds paths as seeds (1 by 1)")
      ),
      "Control of the output:" % (
          cfg.debug << option("--debug") % "enable extensive debug output",
          cfg.draw  << option("--draw")  % "draw pictures around the interesting edges",
          cfg.rescore  << option("--rescore")  % "rescore paths via HMMer",
          cfg.annotate_graph << option("--annotate-graph") % "emit paths in GFA graph"
      ),
      "HMMER options (used for seeding and rescoring):" % (
          cfg.hcfg.acc     << option("--acc")          % "prefer accessions over names in output",
          cfg.hcfg.noali   << option("--noali")        % "don't output alignments, so output is smaller",
          // Control of reporting thresholds
          (option("-E") & number("value", cfg.hcfg.E))        % "report sequences <= this E-value threshold in output",
          (option("-T") & number("value", cfg.hcfg.T))        % "report sequences >= this score threshold in output",
          (option("--domE") & number("value", cfg.hcfg.domE)) % "report domains <= this E-value threshold in output",
          (option("--domT") & number("value", cfg.hcfg.domT)) % "report domains >= this score cutoff in output",
          // Control of inclusion (significance) thresholds
          (option("-incE") & number("value", cfg.hcfg.incE))       % "consider sequences <= this E-value threshold as significant",
          (option("-incT") & number("value", cfg.hcfg.incT))       % "consider sequences >= this score threshold as significant",
          (option("-incdomE") & number("value", cfg.hcfg.incdomE)) % "consider domains <= this E-value threshold as significant",
          (option("-incdomT") & number("value", cfg.hcfg.incdomT)) % "consider domains >= this score threshold as significant",
          // Model-specific thresholding for both reporting and inclusion
          cfg.hcfg.cut_ga  << option("--cut_ga")       % "use profile's GA gathering cutoffs to set all thresholding",
          cfg.hcfg.cut_nc  << option("--cut_nc")       % "use profile's NC noise cutoffs to set all thresholding",
          cfg.hcfg.cut_tc  << option("--cut_tc")       % "use profile's TC trusted cutoffs to set all thresholding",
          // Control of acceleration pipeline
          cfg.hcfg.max     << option("--max")             % "Turn all heuristic filters off (less speed, more power)",
          (option("--F1") & number("value", cfg.hcfg.F1)) % "Stage 1 (MSV) threshold: promote hits w/ P <= F1",
          (option("--F2") & number("value", cfg.hcfg.F2)) % "Stage 2 (Vit) threshold: promote hits w/ P <= F2",
          (option("--F3") & number("value", cfg.hcfg.F3)) % "Stage 3 (Fwd) threshold: promote hits w/ P <= F3"
      ),
      "Developer options:" % (
          cfg.parallel_component_processing << option("--parallel-components") % "process connected components of neighborhood subgraph in parallel [default: false]",
          (option("--max-insertion-length") & integer("x", cfg.max_insertion_length)) % "maximal allowed number of successive I-emissions [default: 30]",
          (option("--expand-coef") & number("value", cfg.expand_coef)) % "overhang expansion coefficient for neighborhood search [default: 2]",
          (option("--expand-const") & integer("value", cfg.expand_const)) % "const addition to overhang values for neighborhood search [default: 20]",
          (option("--no-top-score-filter").set(cfg.state_limits_coef, size_t(100500))) % "disable top score Event Graph vertices filter [default: false]",
          option("--no-fast-forward").set(cfg.use_experimental_i_loop_processing, 0) % "disable fast forward in I-loops processing [default: false]",
          // cfg.disable_depth_filter << option("--disable-depth-filter") % "disable depth filter",  // TODO restore this option
          (option("--known-sequences") & value("filename", cfg.known_sequences)) % "FASTA file with known sequnces that should be definitely found",
          cfg.export_event_graph << option("--export-event-graph") % "export event graph in binary format"
      )
  );

  if (!parse(argc, argv, cli)) {
    std::cout << make_man_page(cli, argv[0]);
    exit(1);
  }

  cfg.output_dir = output_dir;
}

void DrawComponent(const omnigraph::GraphComponent<debruijn_graph::ConjugateDeBruijnGraph> &component,
                   const debruijn_graph::ConjugateDeBruijnGraph &graph,
                   const std::string &prefix,
                   const std::vector<debruijn_graph::EdgeId> &match_edges) {
    using namespace visualization;
    using namespace visualization::visualization_utils;
    using namespace debruijn_graph;

    // FIXME: This madness needs to be refactored
    graph_labeler::StrGraphLabeler<ConjugateDeBruijnGraph> tmp_labeler1(graph);
    graph_labeler::CoverageGraphLabeler<ConjugateDeBruijnGraph> tmp_labeler2(graph);
    graph_labeler::CompositeLabeler<ConjugateDeBruijnGraph> labeler{tmp_labeler1, tmp_labeler2};

    auto colorer = graph_colorer::DefaultColorer(graph);
    auto edge_colorer = std::make_shared<graph_colorer::CompositeEdgeColorer<ConjugateDeBruijnGraph>>("black");
    edge_colorer->AddColorer(colorer);
    edge_colorer->AddColorer(std::make_shared<graph_colorer::SetColorer<ConjugateDeBruijnGraph>>(graph, match_edges, "green"));
    std::shared_ptr<graph_colorer::GraphColorer<ConjugateDeBruijnGraph>>
            resulting_colorer = std::make_shared<graph_colorer::CompositeGraphColorer<Graph>>(colorer, edge_colorer);

    WriteComponent(component,
                   prefix + ".dot",
                   resulting_colorer,
                   labeler);
}

using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;
using debruijn_graph::ConjugateDeBruijnGraph;

using SuperpathIndex = superpath_index::SuperpathIndex<EdgeId>;


template <typename T>
class PseudoVector {
public:
    size_t size() const { return size_; }
    bool empty() const { return size() == 0; }
    T operator[](size_t i) const { return function_(i); }
    PseudoVector(size_t size, const std::function<T(size_t)> function) : size_{size}, function_{function} {}

private:
    size_t size_;
    std::function<T(size_t)> function_;
};

template <typename StringArray>
auto ScoreSequences(const StringArray &seqs,
                    const std::vector<std::string> &refs,
                    const hmmer::HMM &hmm, const PathracerConfig &cfg) {
    DEBUG("ScoreSequences started");
    bool hmm_in_aas = hmm.abc()->K == 20;
    hmmer::HMMMatcher matcher(hmm, cfg.hcfg);

    if (!hmm_in_aas) {
        DEBUG("HMM in nucleotides");
    } else {
        DEBUG("HMM in amino acids");
    }

    for (size_t i = 0; i < seqs.size(); ++i) {
        std::string seq = seqs[i];
        if (seq.size() < 20)
            continue;
        std::string ref = refs.size() > i ? refs[i] : std::to_string(i);
        if (!hmm_in_aas) {
            matcher.match(ref.c_str(), seq.c_str());
        } else {
            for (size_t shift = 0; shift < 3; ++shift) {
                std::string ref_shift = ref + "/" + std::to_string(shift);
                std::string seq_aas = aa::translate(seq.c_str() + shift);
                matcher.match(ref_shift.c_str(), seq_aas.c_str());
            }
        }
    }

    matcher.summarize();
    return matcher;
}

using PathAlnInfo = std::vector<std::pair<size_t, std::pair<int, int>>>;

template <typename StringArray>
PathAlnInfo GetOverhangs(const hmmer::HMMMatcher &matcher,
                         const StringArray &seqs,
                         const hmmer::HMM &hmm) {
    // TODO Move this logic to ScoreSequences()
    // we need only alphabet size (actually aa/nt flag) from hmm
    // and only lengths from the initial seqs
    // Lets store all this stuff into matcher object or its wrapper
    bool hmm_in_aas = hmm.abc()->K == 20;

    PathAlnInfo matches;
    for (const auto &hit : matcher.hits()) {
        if (!hit.reported() || !hit.included())
            continue;

        const std::string name = hit.name();
        size_t id = std::stoull(name);  // Slash and everything after is ignored automatically
        size_t slash_pos = name.find('/');
        VERIFY((slash_pos == std::string::npos ) ^ (hmm_in_aas));
        int shift = hmm_in_aas ? static_cast<int>(std::strtol(name.c_str() + slash_pos + 1, nullptr, 10)) : 0;
        VERIFY(0 <= shift && shift < 3);  // shift should be 0, 1, or 3
        size_t seqlen = seqs[id].length();
        VERIFY(seqlen >= 3);

        for (const auto &domain : hit.domains()) {
            // Calculate HMM overhang
            std::pair<int, int> seqpos = domain.seqpos();
            std::pair<int, int> hmmpos = domain.hmmpos();

            int roverhang = static_cast<int>(domain.M() - hmmpos.second) - static_cast<int>(domain.L() - seqpos.second);
            int loverhang = static_cast<int>(hmmpos.first) - static_cast<int>(seqpos.first);

            if (hmm_in_aas) {
                loverhang *= 3;
                roverhang *= 3;
                // Take shift and sequence length into account
                int extra_left = shift;
                int extra_right = static_cast<int>(seqlen - shift) % 3;
                loverhang -= extra_left;
                roverhang -= extra_right;
            }

            matches.push_back({id, std::make_pair(loverhang, roverhang)});
        }
    }
    DEBUG("Total matched sequences: " << matches.size());

    return matches;
}

std::string PathToString(const std::vector<EdgeId>& path,
                         const ConjugateDeBruijnGraph &graph) {
    if (path.size() == 0) {
        return "";
    }

    std::string res = graph.EdgeNucls(path[0]).str();
    for (size_t i = 1; i < path.size(); ++i) {
        const auto &e = path[i];
        res = res + graph.EdgeNucls(e).Last(graph.length(e)).str();
        size_t k = graph.k();
        VERIFY(graph.EdgeNucls(path[i - 1]).Last(k).str() == graph.EdgeNucls(path[i]).First(k).str());
    }

    return res;
}

PathAlnInfo MatchedPaths(const std::vector<std::vector<EdgeId>> &paths,
                         const ConjugateDeBruijnGraph &graph,
                         const hmmer::HMM &hmm, const PathracerConfig &cfg) {
    DEBUG("MatchedPaths started");
    auto get = [&](size_t i) -> std::string {
        return PathToString(paths[i], graph);
    };

    PseudoVector<std::string> seqs(paths.size(), get);
    auto matcher = ScoreSequences(seqs, {}, hmm, cfg);

    auto matched = GetOverhangs(matcher, seqs, hmm);
    INFO(matched.size() << " matched edges found");

    if (matched.size() && cfg.debug) {
        int textw = 120;
        #pragma omp critical(console)
        {
            p7_tophits_Targets(stdout, matcher.top_hits(), matcher.pipeline(), textw);
            p7_tophits_Domains(stdout, matcher.top_hits(), matcher.pipeline(), textw);
            p7_pli_Statistics(stdout, matcher.pipeline(), nullptr);
        }
    }

    return matched;
}

using EdgeAlnInfo = std::vector<std::pair<EdgeId, std::pair<int, int>>>;
using graph_t = debruijn_graph::ConjugateDeBruijnGraph;

size_t path_length(const graph_t &graph, const std::vector<EdgeId> &path) {
    if (path.size() == 0)
        return 0;

    size_t sum = 0;
    for (const auto &e : path)
        sum += graph.length(e);

    return sum + graph.k();
}

EdgeAlnInfo PathAlignments2EdgeAlignments(const PathAlnInfo &painfo,
                                          const std::vector<std::vector<EdgeId>> &paths,
                                          const graph_t &graph) {
    DEBUG("PathAlignments2EdgeAlignments started");
    EdgeAlnInfo result;
    size_t k = graph.k();

    for (const auto &aln : painfo) {
        const auto &path = paths[aln.first];
        size_t path_len = path_length(graph, path);
        size_t position = 0;
        for (const auto &e : path) {
            size_t edge_len = graph.length(e) + k;
            int loverhang = aln.second.first + static_cast<int>(position);
            int roverhang = aln.second.second + static_cast<int>(path_len - position - edge_len);

            if (-loverhang >= static_cast<int>(edge_len) || -roverhang >= static_cast<int>(edge_len)) {
                // Do nothing, edge lays outside the matched region
            } else {
                result.push_back({e, {loverhang, roverhang}});
            }

            position += graph.length(e);
        }
        VERIFY(position + k == path_len);
    }

    return result;
}

void OutputMatches(const hmmer::HMM &hmm, const hmmer::HMMMatcher &matcher, const std::string &filename,
                   const std::string &format = "tblout") {
    P7_HMM *p7hmm = hmm.get();
    FILE *fp = fopen(filename.c_str(), "w");
    if (format == "domtblout") {
        p7_tophits_TabularDomains(fp, p7hmm->name, p7hmm->acc, matcher.top_hits(), matcher.pipeline(), true);
        // TODO Output tail
    } else if (format == "tblout") {
        p7_tophits_TabularTargets(fp, p7hmm->name, p7hmm->acc, matcher.top_hits(), matcher.pipeline(), true);
        // TODO Output tail
    } else if (format == "pfamtblout") {
        p7_tophits_TabularXfam(fp, p7hmm->name, p7hmm->acc, matcher.top_hits(), matcher.pipeline());
    } else {
        FATAL_ERROR("unknown output format");
    }
    fclose(fp);
}

std::vector<hmmer::HMM> ParseHMMFile(const std::string &filename) {
    /* Open the query profile HMM file */
    hmmer::HMMFile hmmfile(filename);
    if (!hmmfile.valid()) {
        FATAL_ERROR("Error opening HMM file " << filename);
    }

    std::vector<hmmer::HMM> hmms;

    while (auto hmmw = hmmfile.read())
        hmms.emplace_back(std::move(hmmw.get()));

    if (hmms.empty()) {
        FATAL_ERROR("Error reading HMM file " << filename);
    }

    return hmms;
}

std::vector<hmmer::HMM> ParseFASTAFile(const std::string &filename, enum Mode mode,
                                       bool ignore_name) {
    std::vector<hmmer::HMM> res;
    hmmer::HMMSequenceBuilder builder(mode == Mode::nucl ? hmmer::Alphabet::DNA : hmmer::Alphabet::AMINO,
                                      hmmer::ScoreSystem::Default);

    ESL_ALPHABET   *abc  = esl_alphabet_Create(mode == Mode::nucl ? eslDNA : eslAMINO);
    ESL_SQ         *qsq  = esl_sq_CreateDigital(abc);
    ESL_SQFILE     *qfp  = NULL;
    const char *qfile = filename.c_str();

    // Open the query sequence file in FASTA format
    int status = esl_sqfile_Open(qfile, eslSQFILE_FASTA, NULL, &qfp);
    if (status == eslENOTFOUND) {
        FATAL_ERROR("No such file " << filename);
    } else if (status == eslEFORMAT) {
        FATAL_ERROR("Format of " << filename << " unrecognized.");
    } else if (status == eslEINVAL) {
        FATAL_ERROR("Can't autodetect stdin or .gz.");
    } else if (status != eslOK) {
        FATAL_ERROR("Open of " << filename << " failed, code " << status);
    }

    // For each sequence, build a model and save it.
    size_t id = 0;
    while ((status = esl_sqio_Read(qfp, qsq)) == eslOK) {
        INFO("Converting " << qsq->name << ", len: " << qsq->n);
        if (ignore_name)
            esl_sq_SetName(qsq, std::to_string(id++).c_str());
        res.push_back(builder.from_string(qsq));
        esl_sq_Reuse(qsq);
    }
    if (status != eslEOF) {
        FATAL_ERROR("Unexpected error " << status << " reading sequence file " << filename);
    }

    esl_sq_Destroy(qsq);
    esl_sqfile_Close(qfp);
    esl_alphabet_Destroy(abc);

    return res;
}

using MappingF = std::function<std::string(EdgeId)>;

template <typename Path>
std::string edgepath2str(const Path &path,
                         const MappingF &mapping_f) {
    std::vector<std::string> mapped_ids;
    for (const auto &id : path) {
        mapped_ids.push_back(mapping_f(id));
    }
    return join_sep(mapped_ids, "_");
}

template <typename Container>
auto EdgesToSequences(const Container &entries,
                      const debruijn_graph::ConjugateDeBruijnGraph &graph,
                      const MappingF &mapping_f) {
    std::vector<std::pair<std::string, std::string>> ids_n_seqs;
    for (const auto &entry : entries) {
        std::string id = edgepath2str(entry, mapping_f);
        std::string seq = debruijn_graph::MergeSequences(graph, entry).str();
        ids_n_seqs.push_back({id, seq});
    }

    std::sort(ids_n_seqs.begin(), ids_n_seqs.end());
    return ids_n_seqs;
}

std::string SuperPathInfo(const std::vector<EdgeId> &path,
                          const SuperpathIndex &index,
                          const MappingF &mapping_f) {
    std::vector<std::string> results;
    for (const auto &p : index.query(path)) {
        std::stringstream super_path_info;
        super_path_info << edgepath2str(index[p.first], mapping_f) << ":" << p.first << "/" << p.second;
        results.push_back(super_path_info.str());
    }

    return join_sep(results, ",");
}

template <typename Container>
void ExportEdges(const Container &entries,
                 const debruijn_graph::ConjugateDeBruijnGraph &graph,
                 const SuperpathIndex &scaffold_paths,
                 const std::string &filename,
                 const MappingF &mapping_f) {
    if (entries.size() == 0)
        return;

    std::ofstream o(filename, std::ios::out);

    for (const auto &path : entries) {
        std::string id = edgepath2str(path, mapping_f);
        std::string seq = debruijn_graph::MergeSequences(graph, path).str();
        // FIXME return sorting like in EdgesToSequences
        o << ">" << id << "|ScaffoldSuperpaths=" << SuperPathInfo(path, scaffold_paths, mapping_f) << "\n";
        io::WriteWrapped(seq, o);
    }
}

std::vector<EdgeId> conjugate_path(const std::vector<EdgeId> &path,
                                   const debruijn_graph::ConjugateDeBruijnGraph &g) {
    std::vector<EdgeId> cpath;
    for (auto it = path.crbegin(), e = path.crend(); it != e; ++it) {
        cpath.push_back(g.conjugate(*it));
    }
    return cpath;
}

void SaveResults(const hmmer::HMM &hmm, const ConjugateDeBruijnGraph & /* graph */,
                 const PathracerConfig &cfg,
                 const std::vector<HMMPathInfo> &results,
                 const SuperpathIndex &scaffold_paths,
                 const MappingF &mapping_f) {
    const P7_HMM *p7hmm = hmm.get();  // TODO We use only hmm name from this object, may be we should just pass the name itself
    bool hmm_in_aas = hmm.abc()->K == 20;

    INFO("Total " << results.size() << " resultant paths extracted");

    if (results.empty())
        return;

#pragma omp critical(save_results)
    {
        std::ofstream o_seqs(cfg.output_dir / (std::string{p7hmm->name} + ".seqs.fa"), std::ios::out);
        std::ofstream o_nucs;
        if (hmm_in_aas)
            o_nucs.open(cfg.output_dir / (std::string{p7hmm->name} + ".nucs.fa"), std::ios::out);

        for (const auto &result : results) {
            if (result.seq.size() == 0)
                continue;
            auto scaffold_path_info = SuperPathInfo(result.path, scaffold_paths, mapping_f);
            std::stringstream component_info;
            const std::string edge_prefix = "edge_";
            if (result.label.size() && result.label.substr(0, edge_prefix.size()) != edge_prefix) {
                component_info << scaffold_paths[std::stoll(result.label)] << ":" << result.label;
            } else {
                component_info << result.label;
            }
            std::stringstream header;
            std::string path_id = edgepath2str(result.path, mapping_f);
            header << ">Score=" << result.score << "|Edges=" << path_id << "|Position=" << result.pos << "|Alignment=" << result.alignment << "|ScaffoldSuperpaths=" << scaffold_path_info << "|OriginScaffoldPath=" << component_info.str() << '\n';

            o_seqs << header.str();
            io::WriteWrapped(result.seq, o_seqs);

            if (hmm_in_aas) {
                o_nucs << header.str();
                io::WriteWrapped(result.nuc_seq, o_nucs);
            }
        }
    }
}

void Rescore(const hmmer::HMM &hmm, const ConjugateDeBruijnGraph &graph,
             const PathracerConfig &cfg,
             const std::vector<HMMPathInfo> &results,
             const SuperpathIndex &scaffold_paths,
             const MappingF &mapping_f) {
    P7_HMM *p7hmm = hmm.get();

    std::unordered_set<std::vector<EdgeId>> to_rescore;

    for (const auto &result : results) {
        if (result.path.size() == 0)
            continue;

        to_rescore.insert(result.path);
    }

    INFO("Total " << to_rescore.size() << " local paths to rescore");

    std::unordered_map<std::vector<EdgeId>, double> path2score;
    for (const auto &result : results) {
        // results are already ordered by score
        const auto &path = result.path;
        if (!path2score.count(path)) {
            path2score[path] = result.score;
        }
    }

    std::vector<std::vector<EdgeId>> to_rescore_ordered(to_rescore.cbegin(), to_rescore.cend());
    to_rescore.clear();
    sort_by(to_rescore_ordered.begin(), to_rescore_ordered.end(),
            [&](const auto &path){ return std::make_tuple(-path2score[path], path); });

    // TODO export paths along with their best score
#pragma omp critical(export_edges)
    ExportEdges(to_rescore_ordered, graph, scaffold_paths,
                cfg.output_dir / (std::string(p7hmm->name) + ".edges.fa"),
                mapping_f);

    std::vector<std::string> seqs_to_rescore;
    std::vector<std::string> refs_to_rescore;
    for (const auto &kv : EdgesToSequences(to_rescore_ordered, graph, mapping_f)) {
        refs_to_rescore.push_back(kv.first);
        seqs_to_rescore.push_back(kv.second);
    }

    if (cfg.rescore) {
        INFO("Rescore edges using HMMER");
        auto matcher = ScoreSequences(seqs_to_rescore, refs_to_rescore, hmm, cfg);
        INFO("Edges rescored, output");
#pragma omp critical(rescore)
        {
            OutputMatches(hmm, matcher, cfg.output_dir / (std::string(p7hmm->name) + ".tblout"), "tblout");
            OutputMatches(hmm, matcher, cfg.output_dir / (std::string(p7hmm->name) + ".domtblout"), "domtblout");
            OutputMatches(hmm, matcher, cfg.output_dir / (std::string(p7hmm->name) + ".pfamtblout"), "pfamtblout");
        }
    }
}

using GraphCursor = DebruijnGraphCursor;

auto ConnCompsFromEdgesMatches(const EdgeAlnInfo &matched_edges, const graph_t &graph,
                               double expand_coef, int expand_const,
                               bool parallel_component_processing) {
    DEBUG("ConnCompsFromEdgesMatches started");
    using GraphCursor = DebruijnGraphCursor;
    std::vector<std::pair<GraphCursor, size_t>> left_queries, right_queries;
    std::unordered_set<GraphCursor> cursors;
    for (const auto &kv : matched_edges) {
        EdgeId e = kv.first;
        int loverhang = kv.second.first + expand_const;
        int roverhang = kv.second.second + expand_const;

        if (loverhang > 0) {
            for (const auto &start : GraphCursor::get_cursors(graph, e, 0)) {
                left_queries.push_back({start, static_cast<size_t>(loverhang * expand_coef)});
            }
        }

        size_t len = graph.length(e) + graph.k();
        DEBUG("Edge length: " << len <<"; edge id " << e.int_id() << " edge overhangs: " << loverhang << " " << roverhang);
        if (roverhang > 0) {
            for (const auto &end : GraphCursor::get_cursors(graph, e, len - 1)) {
                right_queries.push_back({end, static_cast<size_t>(roverhang * expand_coef)});
            }
        }

        for (size_t i = std::max(0, -loverhang); i < len - std::max(0, -roverhang); ++i) {
            auto position_cursors = GraphCursor::get_cursors(graph, e, i);
            cursors.insert(std::make_move_iterator(position_cursors.begin()), std::make_move_iterator(position_cursors.end()));
        }
    }

    INFO("Getting leftside neighbourhoods");
    auto left_cursors = neighborhood(left_queries, &graph, /* forward */ false);
    INFO("Getting rightside neighbourhoods");
    auto right_cursors = neighborhood(right_queries, &graph, /* forward */ true);

    INFO("Exporting cursors");
    cursors.insert(left_cursors.cbegin(), left_cursors.cend());
    cursors.insert(right_cursors.cbegin(), right_cursors.cend());

    std::vector<GraphCursor> cursors_vector(cursors.cbegin(), cursors.cend());
    auto cursor_conn_comps = parallel_component_processing ? cursor_connected_components(cursors_vector, &graph) : fake_cursor_connected_components(cursors_vector, &graph);
    std::stable_sort(cursor_conn_comps.begin(), cursor_conn_comps.end(),
                     [](const auto &c1, const auto &c2) { return c1.size() > c2.size(); });

    return cursor_conn_comps;
}

void TraceHMM(const hmmer::HMM &hmm,
              const debruijn_graph::ConjugateDeBruijnGraph &graph, const std::vector<EdgeId> &edges,
              const SuperpathIndex &scaffold_paths,
              const PathracerConfig &cfg,
              std::vector<HMMPathInfo> &results) {
    const P7_HMM *p7hmm = hmm.get();

    INFO("Query:         " << p7hmm->name << "  [M=" << p7hmm->M << "]");
    if (p7hmm->acc) {
        INFO("Accession:     " << p7hmm->acc);
    }
    if (p7hmm->desc) {
        INFO("Description:   " << p7hmm->desc);
    }

    auto fees = hmm::fees_from_hmm(p7hmm, hmm.abc());
    fees.state_limits.l25 = 1000000 * cfg.state_limits_coef;
    fees.state_limits.l100 = 50000 * cfg.state_limits_coef;
    fees.state_limits.l500 = 10000 * cfg.state_limits_coef;
    if (cfg.minimal_match_length <= 1.0) {
        fees.minimal_match_length = static_cast<size_t>(cfg.minimal_match_length * static_cast<double>(fees.M));
    } else {
        fees.minimal_match_length = static_cast<size_t>(cfg.minimal_match_length);
    }

    fees.max_insertion_length = cfg.max_insertion_length;
    fees.local = cfg.local;
    fees.use_experimental_i_loop_processing = cfg.use_experimental_i_loop_processing;
    fees.frame_shift_cost = fees.all_matches_score() / static_cast<double>(fees.M) / cfg.indel_rate / 3;
    if (cfg.indel_rate > 0.0 && fees.minimal_match_length) {
        WARN("Current depth filter implementation does not support alignment with frame shifts! minimal_match_length will be set to 0");
        fees.minimal_match_length = 0;
        // FIXME Implement proper filter or at least disable it in more gentle way (leaving post-processing filter)
    }

    INFO("HMM consensus: " << fees.consensus);
    INFO("HMM " << p7hmm->name << " has " << fees.count_negative_loops() << " positive-score I-loops over " << fees.ins.size());
    INFO("All-matches consensus sequence score: " << fees.all_matches_score());
    INFO("Empty sequence score: " << fees.empty_sequence_score());

    std::vector<std::vector<GraphCursor>> cursor_conn_comps;
    std::vector<std::string> component_names;

    if (cfg.seed_mode == SeedMode::scaffolds_one_by_one) {
        for (size_t idx = 0; idx < scaffold_paths.size(); ++idx) {
            const auto &path = scaffold_paths[idx];
            std::vector<std::vector<EdgeId>> paths = {path};
            auto matched_paths = MatchedPaths(paths, graph, hmm, cfg);
            if (!matched_paths.size()) {
                // path not matched
                continue;
            }
            auto matched_edges = PathAlignments2EdgeAlignments(matched_paths, paths, graph);
            auto cursor_conn_comps_local = ConnCompsFromEdgesMatches(matched_edges, graph,
                                                                     cfg.expand_coef, cfg.expand_const,
                                                                     cfg.parallel_component_processing);
            cursor_conn_comps.insert(cursor_conn_comps.end(), cursor_conn_comps_local.cbegin(), cursor_conn_comps_local.cend());
            for (size_t cmp_idx = 0; cmp_idx < cursor_conn_comps_local.size(); ++cmp_idx) {
                // TODO add cmp_idx? (it could not be trivial!!!)
                component_names.push_back(std::to_string(idx));
            }
            // TODO be more verbose
        }
    } else if (cfg.seed_mode == SeedMode::edges_one_by_one) {
        for (const auto &e : edges) {
            std::vector<std::vector<EdgeId>> paths = {{e}};
            auto matched_paths = MatchedPaths(paths, graph, hmm, cfg);
            if (!matched_paths.size()) {
                // path not matched
                continue;
            }
            auto matched_edges = PathAlignments2EdgeAlignments(matched_paths, paths, graph);
            auto cursor_conn_comps_local = ConnCompsFromEdgesMatches(matched_edges, graph,
                                                                     cfg.expand_coef, cfg.expand_const,
                                                                     cfg.parallel_component_processing);
            cursor_conn_comps.insert(cursor_conn_comps.end(), cursor_conn_comps_local.cbegin(), cursor_conn_comps_local.cend());
            for (size_t cmp_idx = 0; cmp_idx < cursor_conn_comps_local.size(); ++cmp_idx) {
                // TODO add cmp_idx? (it could not be trivial!!!)
                component_names.push_back("edge_" + std::to_string(e.int_id()));
            }
        }
    } else if (cfg.seed_mode == SeedMode::edges) {
        std::vector<std::vector<EdgeId>> paths;
        // Fill paths by single edges
        for (const auto &e : edges) {
            paths.push_back(std::vector<EdgeId>({e}));
        }
        auto matched_paths = MatchedPaths(paths, graph, hmm, cfg);
        auto matched_edges = PathAlignments2EdgeAlignments(matched_paths, paths, graph);
        cursor_conn_comps = ConnCompsFromEdgesMatches(matched_edges, graph,
                                                      cfg.expand_coef, cfg.expand_const,
                                                      cfg.parallel_component_processing);
    } else if (cfg.seed_mode == SeedMode::scaffolds) {
        std::vector<std::vector<EdgeId>> paths;
        // Fill paths by paths read from GFA
        paths.insert(paths.end(), scaffold_paths.cbegin(), scaffold_paths.cend());
        auto matched_paths = MatchedPaths(paths, graph, hmm, cfg);
        auto matched_edges = PathAlignments2EdgeAlignments(matched_paths, paths, graph);
        cursor_conn_comps = ConnCompsFromEdgesMatches(matched_edges, graph,
                                                      cfg.expand_coef, cfg.expand_const,
                                                      cfg.parallel_component_processing);
    } else if (cfg.seed_mode == SeedMode::edges_scaffolds) {
        std::vector<std::vector<EdgeId>> paths;
        // Fill paths by single edges
        for (const auto &e : edges) {
            paths.push_back(std::vector<EdgeId>({e}));
        }
        // Fill paths by paths read from GFA
        paths.insert(paths.end(), scaffold_paths.cbegin(), scaffold_paths.cend());
        auto matched_paths = MatchedPaths(paths, graph, hmm, cfg);
        auto matched_edges = PathAlignments2EdgeAlignments(matched_paths, paths, graph);
        cursor_conn_comps = ConnCompsFromEdgesMatches(matched_edges, graph,
                                                      cfg.expand_coef, cfg.expand_const,
                                                      cfg.parallel_component_processing);
    } else if (cfg.seed_mode == SeedMode::exhaustive) {
        cursor_conn_comps.resize(1);
        auto &cursors = cursor_conn_comps[0];

        for (EdgeId edge : graph.edges()) {
            size_t len = graph.length(edge) + graph.k();
            for (size_t i = 0; i < len; ++i) {
                auto position_cursors = GraphCursor::get_cursors(graph, edge, i);
                cursors.insert(cursors.end(),
                               std::make_move_iterator(position_cursors.begin()), std::make_move_iterator(position_cursors.end()));
            }
        }
    }

    if (!cursor_conn_comps.size()) {
        WARN("No components to process!");
        return;
    }
    VERIFY(cursor_conn_comps.size());

    INFO("The number of connected components: " << cursor_conn_comps.size());
    std::vector<size_t> cursor_conn_comps_sizes;
    for (const auto &comp : cursor_conn_comps) {
        cursor_conn_comps_sizes.push_back(comp.size());
    }
    INFO("Connected component sizes: " << cursor_conn_comps_sizes);

    auto run_search = [&fees, &p7hmm, &cfg, &graph](const auto &cached_context,
                                                    const auto &cursors,
                                                    const auto context,
                                                    size_t top,
                                                    std::vector<HMMPathInfo> &local_results,
                                                    const std::string &component_name = "") -> void {
        auto cached_cursors = cached_context.Cursors();
        for (const auto &cursor : cached_cursors) {
            DEBUG_ASSERT(check_cursor_symmetry(cursor, &cached_context), main_assert{}, debug_assert::level<2>{});
            // VERIFY(check_cursor_symmetry(cursor, &cached_context));
        }
        auto result = find_best_path(fees, cached_cursors, &cached_context);
        INFO("Collapsing event graph");
        size_t collapsed_count = result.pathlink_mutable()->collapse_all();
        INFO(collapsed_count << " event graph vertices modified");
        VERIFY(collapsed_count == 0);
        INFO("Event graph depth " << result.pathlink()->max_prefix_size());

        if (!cfg.known_sequences.empty()) {
            auto seqs = read_fasta(cfg.known_sequences);
            for (const auto &kv : seqs) {
                const std::string seq = kv.second;
                if (result.pathlink()->has_sequence(seq, &cached_context)) {
                    INFO("Sequence " << kv.first << " found");
                } else {
                    auto has_prefix = [&](size_t l) -> bool { return result.pathlink()->has_sequence(seq.substr(0, l), &cached_context); };
                    auto has_suffix = [&](size_t l) -> bool { return result.pathlink()->has_sequence(seq.substr(seq.size() - l), &cached_context); };
                    auto has_infix = [&](size_t l) -> bool { return result.pathlink()->has_sequence(seq.substr((seq.size() - l) / 2, l), &cached_context); };
                    size_t max_prefix_size = int_max_binsearch<size_t>(has_prefix, 0, seq.length() + 1);
                    size_t max_suffix_size = int_max_binsearch<size_t>(has_suffix, 0, seq.length() + 1);
                    size_t max_infix_size = int_max_binsearch<size_t>(has_infix, 0, seq.length() + 1);
                    INFO("Sequence " << kv.first << " not found, max prefix " << max_prefix_size << ", max suffix " << max_suffix_size << " max infix " << max_infix_size << " over " << seq.size());
                }
            }
        }

        if (cfg.export_event_graph) {
            std::ofstream of(cfg.output_dir /
                             (std::string{"event_graph_"} + p7hmm->name +
                              "_component_" + int_to_hex(hash_value(cursors)) +
                              "_size_" + std::to_string(cursors.size()) +
                              ".bin"));
            // FIXME: We cannot save cursors as we need to have some kind of context to serialize / deserialize shared pointers
            io::binary::BinWrite(of, cursors, cached_context /*, result */);
            INFO("Event graph exported");
        }

        INFO("Extracting top paths");
        auto top_paths = result.top_k(&cached_context, top);
        INFO("Total " << top_paths.size() << " paths extracted");
        bool x_as_m_in_alignment = fees.is_proteomic();
        if (!top_paths.empty()) {
            INFO("Best score in the current component: " << result.best_score());
            INFO("Best sequence in the current component");
            INFO(top_paths.str(0, &cached_context));
            INFO("Alignment: " << compress_alignment(top_paths.alignment(0, fees, &cached_context), x_as_m_in_alignment));
        }

        std::unordered_set<std::tuple<std::vector<EdgeId>, size_t, size_t>> extracted_paths;

        for (const auto& annotated_path : top_paths) {
            VERIFY(annotated_path.path.size());
            std::string seq = annotated_path.str(&cached_context);
            if (seq.length() < fees.minimal_match_length) {
                continue;
            }
            auto unpacked_path = cached_context.UnpackPath(annotated_path.path, cursors);
            auto alignment = compress_alignment(annotated_path.alignment(fees, &cached_context), x_as_m_in_alignment);
            auto nucl_path = to_nucl_path(unpacked_path);
            std::string nucl_seq = pathtree::path2string(nucl_path, context);
            auto edge_path = to_path(nucl_path);
            // DEBUG_ASSERT(check_path_continuity(nucl_path, context), main_assert{}, debug_assert::level<2>{});
            VERIFY(check_path_continuity(nucl_path, context));
            // DEBUG_ASSERT(!edge_path.empty(), main_assert{});
            // auto edge_path_aas = to_path(unpacked_path);
            // if (edge_path != edge_path_aas) {
            //     ERROR("NT: " << edge_path);
            //     ERROR("AA: " << edge_path_aas);
            // }
            // DEBUG_ASSERT(edge_path == edge_path_aas, main_assert{}, debug_assert::level<2>{});
            size_t pos = nucl_path[0].position();
            HMMPathInfo info(p7hmm->name, annotated_path.score, seq, nucl_seq, std::move(edge_path), std::move(alignment),
                             component_name, pos);
            info.trim_first_edges(graph);
            auto tpl = std::make_tuple(info.path, info.pos, info.nuc_seq.length());
            if (!extracted_paths.count(tpl)) {
                local_results.push_back(std::move(info));
                extracted_paths.insert(tpl);
            }
        }
    };

    std::vector<EdgeId> match_edges;
    for (const auto &comp : cursor_conn_comps) {
        for (const auto &cursor : comp)
        match_edges.push_back(cursor.edge());
    }
    remove_duplicates(match_edges);

    auto process_component = [&hmm, &run_search, &cfg, &results, &graph](const auto &component_cursors,
                                                                         const std::string &component_name = "") -> std::unordered_set<std::vector<EdgeId>> {
        assert(!component_cursors.empty());
        INFO("Component size " << component_cursors.size());

        if (component_cursors.size() > cfg.max_size) {
            WARN("The component is too large, skipping");
            return {};
        }

        DEBUG_EXPR(
            for (const auto &cursor : component_cursors) {
                DEBUG_ASSERT(check_cursor_symmetry(cursor, &graph), main_assert{}, debug_assert::level<2>{});
            }
         );

        std::unordered_set<EdgeId> edges;
        for (const auto& cursor : component_cursors) {
            edges.insert(cursor.edge());
        }

        INFO("# edges in the component: " << edges.size());
        DEBUG("Edges: " << edges);

        INFO("Running path search");
        std::vector<HMMPathInfo> local_results;
        std::unordered_set<GraphCursor> component_set(component_cursors.cbegin(), component_cursors.cend());
        auto restricted_context = make_optimized_restricted_cursor_context(component_set, &graph);
        auto restricted_component_cursors = make_optimized_restricted_cursors(component_cursors);

        bool hmm_in_aas = hmm.abc()->K == 20;
        if (hmm_in_aas) {
            if (cfg.indel_rate == 0.0) {
                auto aa_cursors = make_aa_cursors(restricted_component_cursors, &restricted_context);
                CachedCursorContext ccc(aa_cursors, &restricted_context);
                run_search(ccc, aa_cursors, &restricted_context, cfg.top, local_results, component_name);
            } else {
                CachedAACursorContext caacc(restricted_component_cursors, &restricted_context);
                run_search(caacc, restricted_component_cursors, &restricted_context, cfg.top, local_results, component_name);
            }
            // run_search(make_aa_cursors(restricted_component_cursors, &restricted_context), cfg.top, local_results, &restricted_context, component_name);
        } else {
            // run_search(restricted_component_cursors, cfg.top, local_results, &restricted_context, component_name);
            CachedCursorContext ccc(restricted_component_cursors, &restricted_context);
            run_search(ccc, restricted_component_cursors, &restricted_context, cfg.top, local_results, component_name);
        }

        results.insert(results.end(), local_results.begin(), local_results.end());

        std::unordered_set<std::vector<EdgeId>> paths;
        for (const auto& entry : local_results) {
            paths.insert(entry.path);
        }
        return paths;
    };


    //#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < cursor_conn_comps.size(); ++i) {
        const auto &component_cursors = cursor_conn_comps[i];
        const std::string &component_name = component_names.size() ? component_names[i] : "";
        auto paths = process_component(component_cursors, component_name);

        INFO("Total " << paths.size() << " unique edge paths extracted");
        // size_t count = 0;  // FIXME this ad-hoc
        // for (const auto &path : paths) {
        //     INFO("Path length : " << path.size() << " edges " << path);  // FIXME use id_mapper here as well
        //     ++count;
        //     if (count > 1000) break;
        // }

        if (cfg.draw) {
            auto component_name_with_hash = component_name + int_to_hex(hash_value(component_cursors));
            INFO("Construct component as omnigraph-component" << component_name_with_hash);
            auto component = omnigraph::GraphComponent<ConjugateDeBruijnGraph>::FromEdges(graph, edges, true, component_name_with_hash);

            INFO("Writing component " << component_name_with_hash);
            DrawComponent(component, graph, cfg.output_dir / component_name_with_hash, match_edges);

            size_t idx = 0;
            for (const auto &path : paths) {
                INFO("Writing component around path " << idx);
                DrawComponent(component, graph, cfg.output_dir / (component_name_with_hash + "_" + std::to_string(idx)), path);
                ++idx;
            }
        }
    }
}

void hmm_main(const PathracerConfig &cfg,
              const debruijn_graph::ConjugateDeBruijnGraph &graph,
              const std::vector<EdgeId> &edges,
              const std::vector<std::vector<EdgeId>> scaffold_paths,
              std::unordered_set<std::vector<EdgeId>> &to_rescore,
              std::set<std::pair<std::string, std::vector<EdgeId>>> &gfa_paths,
              const std::function<std::string(EdgeId)> &mapping_f) {
    std::vector<hmmer::HMM> hmms;
    if (cfg.mode == Mode::hmm)
        hmms = ParseHMMFile(cfg.hmmfile);
    else
        hmms = ParseFASTAFile(cfg.hmmfile, cfg.mode, cfg.ignore_names);

    SuperpathIndex scaffold_path_index(scaffold_paths);

    // Filter input hmms
    if (!cfg.queries.empty()) {
        std::unordered_set<std::string> queries(cfg.queries.cbegin(), cfg.queries.cend());
        hmms.erase(std::remove_if(hmms.begin(), hmms.end(),
                                  [&queries](const auto &hmm) -> bool { return !queries.count(hmm.get()->name); }),
                   hmms.end());
    }

    // Outer loop: over each query HMM in <hmmfile>.
    omp_set_num_threads(cfg.threads);
    #pragma omp parallel for schedule(guided)
    for (size_t _i = 0; _i < hmms.size(); ++_i) {
        const auto &hmm = hmms[_i];

        std::vector<HMMPathInfo> results;

        TraceHMM(hmm, graph, edges, scaffold_path_index,
                 cfg, results);

        std::sort(results.begin(), results.end());
        unique_hmm_path_info(results, scaffold_path_index);
        SaveResults(hmm, graph, cfg, results, scaffold_path_index, mapping_f);

        if (cfg.annotate_graph) {
            size_t idx = 0;
            for (const auto &result : results) {
                #pragma omp critical
                {
                    gfa_paths.insert({ std::string(hmm.get()->name) + "_" + std::to_string(idx++) + "_score_" + std::to_string(result.score), result.path });
                }
            }
        }

        Rescore(hmm, graph, cfg, results, scaffold_paths, mapping_f);

        for (const auto &result : results) {
#pragma omp critical
            {
                to_rescore.insert(result.path);
            }
        }
    } // end outer loop over query HMMs
}

int pathracer_main(int argc, char* argv[]) {
    utils::segfault_handler sh;
    utils::perf_counter pc;

    srand(42);
    srandom(42);

    PathracerConfig cfg;
    process_cmdline(argc, argv, cfg);

    int status = mkdir(cfg.output_dir.c_str(), 0775);
    create_console_logger(cfg.output_dir / "pathracer.log", cfg.quiet);

    if (status != 0) {
        if (errno == EEXIST) {
            WARN("Output directory exists: " << cfg.output_dir);
        } else {
            ERROR("Cannot create output directory: " << cfg.output_dir);
            std::exit(1);
        }
    }

    START_BANNER("Graph HMM aligning engine");
    std::string cmd_line = join_sep(adt::make_range(argv, argv + argc), " ");
    INFO("Command line: " << cmd_line);

    // Set memory limit
    const size_t GB = 1 << 30;
    utils::limit_memory(cfg.memory * GB);

    // Stack limit
    size_t stacklimit = stack_limit();
    INFO("Soft stack limit: " << (stacklimit == size_t(-1) ? "UNLIMITED" : std::to_string(stacklimit)));
    INFO("Process ID: " << getpid());

    // Some additional defaults
    if (cfg.local && cfg.minimal_match_length > 0.75 && cfg.minimal_match_length <= 1.0) {
        cfg.minimal_match_length = 0.1;
        INFO("Setting minimal match length threshold in local mode to: " << cfg.minimal_match_length);
    }
    
    using namespace debruijn_graph;

    std::vector<std::vector<EdgeId>> scaffold_paths;
    std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());
    using namespace debruijn_graph;

    debruijn_graph::ConjugateDeBruijnGraph graph(0);
    gfa::GFAReader gfa(cfg.load_from);
    unsigned gfa_k = gfa.to_graph(graph, id_mapper.get());
    INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
    if (gfa_k != -1U) {
        INFO("Detected k: " << gfa_k);
        VERIFY_MSG(gfa_k == 0 || gfa_k % 2 == 1, "k-mer length must be odd");
        VERIFY(graph.k() == gfa_k);
    } else
        FATAL_ERROR("Failed to determine k-mer length");

    scaffold_paths.reserve(gfa.num_paths());
    for (const auto &path : gfa.paths()) {
        VERIFY_MSG(!path.edges.empty(), "empty scaffold path in GFA");
        // We do not support jumps as we need proper graph path,
        // split scaffold into contigs
        std::vector<EdgeId> contig_paths{path.edges.front()};
        for (size_t i = 1; i < path.edges.size(); ++i) {
            EdgeId last = contig_paths.back();
            EdgeId next = path.edges[i];
            if (graph.EdgeEnd(last) == graph.EdgeStart(next)) {
                contig_paths.push_back(next);
            } else {
                scaffold_paths.push_back(contig_paths);
                scaffold_paths.push_back(conjugate_path(contig_paths, graph));
                contig_paths.clear();
                contig_paths.push_back(next);
            }
        }

        scaffold_paths.push_back(contig_paths);
        scaffold_paths.push_back(conjugate_path(contig_paths, graph));
    }

    size_t letters = 0;
    for (EdgeId edge : graph.canonical_edges())
        letters += (graph.length(edge) + graph.k()) * (graph.conjugate(edge) == edge ? 1 : 2);
    INFO("Graph loaded. Total vertices: " << graph.size() << ", edges: " << graph.e_size() << ", letters: " << letters);

    INFO("Total graph paths " << gfa.num_paths() << ", to score: " << scaffold_paths.size() << " paths");

    // Collect all the edges
    std::vector<EdgeId> edges;
    for (std::string &edge : cfg.edges) {
        std::replace(edge.begin(), edge.end(), '^', '\'');
    }
    std::unordered_set<std::string> allowed_edges(cfg.edges.cbegin(), cfg.edges.cend());
    for (EdgeId edge : graph.edges()) {
        if (allowed_edges.empty() || allowed_edges.count((*id_mapper)[edge.int_id()])) {
            edges.push_back(edge);
        }
    }

    std::unordered_set<std::vector<EdgeId>> to_rescore;
    std::set<std::pair<std::string, std::vector<EdgeId>>> gfa_paths;

    const auto mapping_f = [&id_mapper, &graph](EdgeId id) -> std::string { return (*id_mapper)[graph.int_id(id)]; };
    hmm_main(cfg, graph, edges, scaffold_paths, to_rescore, gfa_paths,
             mapping_f);

    if (cfg.rescore) {
        INFO("Total " << to_rescore.size() << " paths to rescore");
        ExportEdges(to_rescore, graph, scaffold_paths,
                    cfg.output_dir / "all.edges.fa",
                    mapping_f);
    }

    if (cfg.annotate_graph) {
        std::string fname = cfg.output_dir / "graph_with_hmm_paths.gfa";
        INFO("Saving annotated graph to " << fname)
        std::ofstream os(fname);
        path_extend::GFAPathWriter gfa_writer(graph, os,
                                              io::MapNamingF<debruijn_graph::ConjugateDeBruijnGraph>(*id_mapper));
        gfa_writer.WriteSegmentsAndLinks();

        for (const auto& entry : gfa_paths) {
            gfa_writer.WritePaths(entry.second, entry.first);
        }
    }

    INFO("Pathracer successfully finished! Thanks for flying us!");
    return 0;
}


#include "cached_aa_cursor.hpp"
struct PathracerSeqFsConfig {
    // std::string load_from = "";
    std::string sequence_file = "";
    std::string hmmfile = "";
    std::string output_dir = "";
    // enum Mode mode = Mode::hmm;
    // enum SeedMode seed_mode = SeedMode::edges_scaffolds;
    // size_t k = 0;
    int threads = 16;
    size_t top = 100;
    size_t max_fs = 10;
    // std::vector<std::string> edges;
    std::vector<std::string> sequences;
    std::vector<std::string> queries;
    size_t max_size = size_t(-1);
    bool debug = false;
    bool draw = false;
    bool rescore = false;
    bool annotate_graph = true;
    double expand_coef = 2.;
    int expand_const = 50;
    size_t state_limits_coef = 1;
    bool local = false;
    bool parallel_component_processing = false;
    bool disable_depth_filter = false;
    size_t memory = 100;  // 100GB
    int use_experimental_i_loop_processing = true;
    std::string known_sequences = "";
    bool export_event_graph = false;
    // double minimal_match_length = 0;
    size_t max_insertion_length = 30;
    double indel_rate = 0.05;
    double cutoff = 0.7;

    hmmer::hmmer_cfg hcfg;
    bool exhaustive = false;
    bool quiet = false;
};

void process_cmdline_seq_fs(int argc, char **argv, PathracerSeqFsConfig &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.hmmfile    << value("input HMM file"),  // TODO support sequence mode
      cfg.sequence_file  << value("input sequence file"),
      required("--output", "-o") & value("output directory", cfg.output_dir)    % "output directory",
      (option("--global").set(cfg.local, false) % "perform global-local (aka glocal) HMM matching [default]") |
      (cfg.local << option("--local") % "perform local-local HMM matching"),
      (option("--indel-rate", "-r") & value("value", cfg.indel_rate)) % "expected rate of nucleotides indels in graph edges [default: 0.05]",
      (option("--top") & integer("N", cfg.top)) % "extract top N paths [default: 100]",
      (option("--max-fs") & integer("MAX #FS", cfg.max_fs)) % "maximal allowed number of frameshifts in a reported sequence [default: 10]",
      (option("--threads", "-t") & integer("NTHREADS", cfg.threads)) % "the number of parallel threads [default: 16]",
      (option("--memory", "-m") & integer("MEMORY", cfg.memory)) % "RAM limit for PathRacer in GB (terminates if exceeded) [default: 100]",
      (option("--cutoff") & value("CUTOFF", cfg.cutoff)) % "bitscore cutoff for reported match; if <= 1 then to be multiplied on GA HMM cutoff [default: 0.7]",
      (option("--queries") & values("queries", cfg.queries)) % "queries names to lookup [default: all queries from input query file]",
      cfg.exhaustive << option("--exhaustive") % "run in exhaustive mode, disable HMM filter",
      (option("--sequences") & values("sequences", cfg.sequences)) % "sequence IDs to process [default: all input sequences]",
      (cfg.quiet << option("--quiet", "-q")) % "be quiet, do not output anything to the console",
      "HMMER options (used for seeding and rescoring):" % (
          cfg.hcfg.acc     << option("--acc")          % "prefer accessions over names in output",
          cfg.hcfg.noali   << option("--noali")        % "don't output alignments, so output is smaller",
          // Control of reporting thresholds
          (option("-E") & number("value", cfg.hcfg.E))        % "report sequences <= this E-value threshold in output",
          (option("-T") & number("value", cfg.hcfg.T))        % "report sequences >= this score threshold in output",
          (option("--domE") & number("value", cfg.hcfg.domE)) % "report domains <= this E-value threshold in output",
          (option("--domT") & number("value", cfg.hcfg.domT)) % "report domains >= this score cutoff in output",
          // Control of inclusion (significance) thresholds
          (option("-incE") & number("value", cfg.hcfg.incE))       % "consider sequences <= this E-value threshold as significant",
          (option("-incT") & number("value", cfg.hcfg.incT))       % "consider sequences >= this score threshold as significant",
          (option("-incdomE") & number("value", cfg.hcfg.incdomE)) % "consider domains <= this E-value threshold as significant",
          (option("-incdomT") & number("value", cfg.hcfg.incdomT)) % "consider domains >= this score threshold as significant",
          // Model-specific thresholding for both reporting and inclusion
          cfg.hcfg.cut_ga  << option("--cut_ga")       % "use profile's GA gathering cutoffs to set all thresholding",
          cfg.hcfg.cut_nc  << option("--cut_nc")       % "use profile's NC noise cutoffs to set all thresholding",
          cfg.hcfg.cut_tc  << option("--cut_tc")       % "use profile's TC trusted cutoffs to set all thresholding",
          // Control of acceleration pipeline
          cfg.hcfg.max     << option("--max")             % "Turn all heuristic filters off (less speed, more power)",
          (option("--F1") & number("value", cfg.hcfg.F1)) % "Stage 1 (MSV) threshold: promote hits w/ P <= F1",
          (option("--F2") & number("value", cfg.hcfg.F2)) % "Stage 2 (Vit) threshold: promote hits w/ P <= F2",
          (option("--F3") & number("value", cfg.hcfg.F3)) % "Stage 3 (Fwd) threshold: promote hits w/ P <= F3"
      ),
      "Developer options:" % (
          (option("--max-insertion-length") & integer("x", cfg.max_insertion_length)) % "maximal allowed number of successive I-emissions [default: 30]",
          (option("--expand-coef") & number("value", cfg.expand_coef)) % "overhang expansion coefficient for neighborhood search [default: 2]",
          (option("--expand-const") & integer("value", cfg.expand_const)) % "const addition to overhang values for neighborhood search [default: 20]",
          (option("--no-top-score-filter").set(cfg.state_limits_coef, size_t(100500))) % "disable top score Event Graph vertices filter [default: false]",
          (option("--no-fast-forward").set(cfg.use_experimental_i_loop_processing, 0)) % "disable fast forward in I-loops processing [default: false]"//,
      )
  );

  if (!parse(argc, argv, cli)) {
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }
}


float max_bitscore(const std::string &seq, hmmer::HMMMatcher &matcher) {
    matcher.reset();
    matcher.match("seq", seq.c_str());
    matcher.summarize();
    float bitscore = -std::numeric_limits<float>::infinity();
    for (const auto &hit : matcher.hits()) {
        if (!hit.reported() || !hit.included())
            continue;

        for (const auto &domain : hit.domains()) {
            bitscore = std::max(bitscore, domain.bitscore());
        }
    }

    return bitscore;
}


std::vector<std::string> split_seq(const std::string &s, const std::unordered_set<char> &delims) {
    std::string current = "";
    std::vector<std::string> result;
    for (char ch : s) {
        if (delims.count(ch) && !current.empty()) {
            result.push_back(std::move(current));
            current = "";
        } else {
            current += ch;
        }
    }
    if (!current.empty()) {
        result.push_back(std::move(current));
    }
    return result;
}

std::vector<float> max_bitscores(const std::vector<std::string> &seqs, hmmer::HMMMatcher &matcher) {
    std::vector<float> result;
    for (const std::string &seq : seqs) {
        float score = max_bitscore(seq, matcher);
        result.push_back(score);
    }
    return result;
}


int aling_fs(int argc, char* argv[]) {
    PathracerSeqFsConfig cfg;
    process_cmdline_seq_fs(argc, argv, cfg);

    utils::segfault_handler sh;
    utils::perf_counter pc;

    srand(42);
    srandom(42);

    int status = mkdir(cfg.output_dir.c_str(), 0775);
    create_console_logger(cfg.output_dir + "/pathracer-seq-fs.log", cfg.quiet);

    if (status != 0) {
        if (errno == EEXIST) {
            WARN("Output directory exists: " << cfg.output_dir);
        } else {
            ERROR("Cannot create output directory: " << cfg.output_dir);
            std::exit(1);
        }
    }

    START_BANNER("Sequence HMM aligning engine");
    std::string cmd_line = join_sep(adt::make_range(argv, argv + argc), " ");
    INFO("Command line: " << cmd_line);

    // Set memory limit
    const size_t GB = 1 << 30;
    utils::limit_memory(cfg.memory * GB);

    // Stack limit
    size_t stacklimit = stack_limit();
    INFO("Soft stack limit: " << (stacklimit == size_t(-1) ? "UNLIMITED" : std::to_string(stacklimit)));
    INFO("Process ID: " << getpid());

    auto hmms = ParseHMMFile(cfg.hmmfile);
    auto all_seqs = read_fasta(cfg.sequence_file);
    INFO("Total " << all_seqs.size() << " input sequences");
    std::vector<std::pair<std::string, std::string>> seqs;

    std::unordered_set<std::string> allowed_seqs(cfg.sequences.cbegin(), cfg.sequences.cend());
    for (const auto &kv : all_seqs) {
        if (!allowed_seqs.empty() && !allowed_seqs.count(kv.first))
            continue;

        seqs.push_back(kv);
        std::string rc = rev_comp(kv.second);
        if (rc != kv.second)
            seqs.push_back({kv.first + "'", std::move(rc)});
    }

    std::unordered_set<std::string> queries(cfg.queries.cbegin(), cfg.queries.cend());

    omp_set_num_threads(cfg.threads);
    #pragma omp parallel for schedule(guided)
    for (size_t i = 0; i < hmms.size(); ++i) {
        const auto &hmm = hmms[i];
        const P7_HMM *p7hmm = hmm.get();

        if (!queries.empty() && !queries.count(p7hmm->name))
            continue;

        auto fees = hmm::fees_from_hmm(p7hmm, hmm.abc());
        VERIFY(fees.is_proteomic());
        fees.state_limits.l25 = 1000000 * cfg.state_limits_coef;
        fees.state_limits.l100 = 50000 * cfg.state_limits_coef;
        fees.state_limits.l500 = 10000 * cfg.state_limits_coef;
        fees.max_insertion_length = cfg.max_insertion_length;
        fees.local = cfg.local;
        fees.use_experimental_i_loop_processing = cfg.use_experimental_i_loop_processing;
        fees.frame_shift_cost = fees.all_matches_score() / static_cast<double>(fees.M) / cfg.indel_rate / 3;
        fees.minimal_match_length = 0;  // FIXME fix depth filter for frame shifts

        INFO("Query:         " << p7hmm->name << "  [M=" << p7hmm->M << "]");
        if (p7hmm->acc) {
            INFO("Accession:     " << p7hmm->acc);
        }
        if (p7hmm->desc) {
            INFO("Description:   " << p7hmm->desc);
        }

        std::ofstream o_seqs(cfg.output_dir + "/" + hmm.get()->name + ".seqs.fa");
        std::ofstream o_nucs(cfg.output_dir + "/" + hmm.get()->name + ".nucs.fa");
        bool at_least_one_match_found = false;

        // Set proper search space
        cfg.hcfg.Z = 3 * seqs.size();
        hmmer::HMMMatcher matcher(hmm, cfg.hcfg);
        for (size_t j = 0; j < seqs.size(); ++j) {
            const auto &id = seqs[j].first;
            const auto &seq = seqs[j].second;

            matcher.reset_top_hits();

            for (size_t shift = 0; shift < 3; ++shift) {
                std::string ref = std::to_string(j) + std::string("/") + std::to_string(shift);
                std::string seq_aas = aa::translate(seq.c_str() + shift);
                matcher.match(ref.c_str(), seq_aas.c_str());
            }

            matcher.summarize();
            std::unordered_set<size_t> indices;
            auto get = [&](size_t i) -> const std::string& {
                return seqs[i].second;
            };

            PseudoVector<std::string> local_seqs(seqs.size(), get);
            auto overs = GetOverhangs(matcher, local_seqs, hmm);
            for (const auto &over : overs) {
                int loverhang = over.second.first + cfg.expand_const;
                int roverhang = over.second.second + cfg.expand_const;
                size_t start = -std::min(loverhang, 0);
                size_t finish = int(seq.length()) + std::min(roverhang, 0);
                INFO("Overhangs: " << loverhang << " " << roverhang);
                INFO("START " << start << " FINISH " << finish);
                for (size_t i = start; i < finish; ++i) {
                    indices.insert(i);
                }
            }

            std::vector<StringCursor> cursors;
            for (size_t i : indices) {
                cursors.push_back(StringCursor(i));
            }

            if (cfg.exhaustive) {
                for (size_t i = 0; i < seq.size(); ++i)
                    cursors.push_back(StringCursor(i));
            }

            if (!cursors.size())
                continue;

            INFO("Sequence: " << id << ", seed hits: " << overs.size());

            std::unordered_set<StringCursor> cursor_set(cursors.cbegin(), cursors.cend());

            auto restricted_context = make_optimized_restricted_cursor_context(cursor_set, &seq);
            auto restricted_component_cursors = make_optimized_restricted_cursors(cursors);

            CachedAACursorContext caacc(restricted_component_cursors, &restricted_context);
            auto all_cursors = caacc.Cursors();
            for (const auto &cursor : all_cursors) {
                VERIFY(check_cursor_symmetry(cursor, &caacc));
            }
            auto result = find_best_path(fees, all_cursors, &caacc);

            // INFO("Collapsing event graph");
            // size_t collapsed_count = result.pathlink_mutable()->collapse_all();
            // INFO(collapsed_count << " event graph vertices modified");
            // VERIFY(collapsed_count == 0);
            // INFO("Event graph depth " << result.pathlink()->max_prefix_size());

            INFO("Extracting top paths");
            auto top_paths = result.top_k(&caacc, cfg.top);
            INFO("Total " << top_paths.size() << " paths extracted");
            bool x_as_m_in_alignment = fees.is_proteomic();
            if (!top_paths.empty()) {
                INFO("Best score in the current component: " << result.best_score());
                INFO("Best sequence in the current component");
                const auto top_string = top_paths.str(0, &caacc);
                INFO(top_string);
                const auto alignment = compress_alignment(top_paths.alignment(0, fees, &caacc), x_as_m_in_alignment);
                INFO("Alignment: " << alignment);
            }
            auto &context = restricted_context;
            for (const auto& annotated_path : top_paths) {
                VERIFY(annotated_path.path.size());
                std::string seq = annotated_path.str(&caacc);
                if (seq.length() < fees.minimal_match_length) {
                    continue;
                }
                auto unpacked_path = caacc.UnpackPath(annotated_path.path, restricted_component_cursors);
                auto alignment = compress_alignment(annotated_path.alignment(fees, &caacc), x_as_m_in_alignment);
                auto nucl_path = to_nucl_path(unpacked_path);
                VERIFY(check_path_continuity(nucl_path, &context));
                std::string nucl_seq = pathtree::path2string(nucl_path, &context);
                size_t pos = nucl_path[0].position();
                HMMPathInfo info(p7hmm->name, annotated_path.score, seq, nucl_seq, {}, std::move(alignment),
                                 "NA", pos);

                auto subseqs = split_seq(seq, {'=', '-'});
                auto subscores = max_bitscores(subseqs, matcher);
                float max_subscore = subscores.empty() ? -std::numeric_limits<float>::infinity() : *std::max_element(subscores.cbegin(), subscores.cend());

                // std::string seq_without_gaps = seq;
                // seq_without_gaps.erase(std::remove_if(seq_without_gaps.begin(), seq_without_gaps.end(), [](char ch) {return ch == '-' || ch == '=';}),
                //                        seq_without_gaps.end());
                //
                std::string seq_without_gaps = join_sep(subseqs, "");
                // VERIFY(seq_without_gaps == seq_without_gaps2);
                float bitscore = max_bitscore(seq_without_gaps, matcher);

                double cutoff = cfg.cutoff <= 1.0 ? cfg.cutoff * p7hmm->cutoff[cfg.local ? p7_GA2 : p7_GA1] : cfg.cutoff;  // FIXME check global and local
                int nindels = int(subseqs.size()) - 1;
                VERIFY(nindels >= 0);
                if (nindels > int(cfg.max_fs))
                    continue;

                if (cfg.cutoff == -1.0) {  // Auto cutoff
                    const double coef = 15.0;
                    double alpha = static_cast<double>(nindels) / p7hmm->M * coef;
                    cutoff = (1 - alpha) * p7hmm->cutoff[cfg.local ? p7_GA2 : p7_GA1];  // FIXME check global and local
                }

                if (bitscore < cutoff)
                    continue;

                std::stringstream header;
                // FIXME report PartialScore correspondent to the currunt much rather than just maximal partial score
                header << ">Score=" << info.score << "|Bitscore=" << bitscore << "|PartialBitscore=" << max_subscore << "|Seq=" << id << "|Position=" << info.pos << "|Frameshifts=" << nindels << "|Alignment=" << info.alignment << '\n';
                o_seqs << header.str();
                io::WriteWrapped(info.seq, o_seqs);
                o_nucs << header.str();
                io::WriteWrapped(info.nuc_seq, o_nucs);
                at_least_one_match_found = true;
            }
        }

        o_nucs.close();
        o_seqs.close();
        if (!at_least_one_match_found) {
            remove((cfg.output_dir + "/" + hmm.get()->name + ".seqs.fa").c_str());
            remove((cfg.output_dir + "/" + hmm.get()->name + ".nucs.fa").c_str());
        }
    }

    INFO("Pathracer successfully finished! Thanks for flying us!");
    return 0;
}

int aling_kmers_main(int argc, char* argv[]) {
    create_console_logger("");
    using namespace clipp;

    std::string hmm_file;
    std::string sequence_file;
    std::string output_file;
    size_t k;

    auto cli =
        (sequence_file << value("input sequence file"),
         hmm_file << value("HMM file"),
         k << integer("k-mer size"),
         required("--output", "-o") & value("output file", output_file) % "output file"
         );

    if (!parse(argc, argv, cli)) {
        std::cout << make_man_page(cli, argv[0]);
        exit(1);
    }

    auto hmms = ParseHMMFile(hmm_file);
    auto seqs = read_fasta(sequence_file);

    hmmer::hmmer_cfg hcfg;
    for (const auto &hmm : hmms) {
    const P7_HMM *p7hmm = hmm.get();
    auto fees = hmm::fees_from_hmm(p7hmm, hmm.abc());
    const size_t state_limits_coef = 100500;
    fees.state_limits.l25 = 1000000 * state_limits_coef;
    fees.state_limits.l100 = 50000 * state_limits_coef;
    fees.state_limits.l500 = 10000 * state_limits_coef;
    fees.minimal_match_length = static_cast<size_t>(0.9 * static_cast<double>(fees.M));
    fees.use_experimental_i_loop_processing = true;

    std::ofstream of(output_file + hmm.get()->name);
    std::vector<double> best_scores;
    for (size_t seq_id = 0; seq_id < seqs.size(); ++seq_id) {
        std::string seq = seqs[seq_id].second;
        std::string id = seqs[seq_id].first;

        if (seq.find('X') != std::string::npos) {
            INFO("Stop-codon found, skipping the sequence");
            continue;
        }
        // remove -
        seq.erase(std::remove(seq.begin(), seq.end(), '-'), seq.end());
        // double score = score_sequence(fees, seq);  //FIXME initial cursor in score_sequence
        /*
        double score = score_subsequence(fees, seq);
        if (score < 0) {
            continue;
        }
        of << ">" << id << "|Score=" << score <<  "\n";
        io::WriteWrapped(seq, of);
        */

        hmmer::HMMMatcher matcher(hmm, hcfg);
        // Split into k-mers
        if (k > seq.size()) {
            continue;
        }
        for (size_t i = 0; i < seq.size() - k + 1; ++i) {
            std::string kmer = seq.substr(i, k);
            std::string ref = std::string("seq_") + std::to_string(seq_id) + std::string("_kmer_") + std::to_string(i);
            matcher.match(ref.c_str(), kmer.c_str());
            // if (!hmm_in_aas) {
            //     matcher.match(ref.c_str(), kmer.c_str());
            // } else {
            //     VERIFY(seq.size() >= 2);
            //     for (size_t shift = 0; shift < 3; ++shift) {
            //         std::string ref_shift = ref + "/" + std::to_string(shift);
            //         std::string seq_aas = aa::translate(kmer.c_str() + shift);
            //         matcher.match(ref_shift.c_str(), seq_aas.c_str());
            //     }
            // }
        }

        size_t hit_count = 0;

        matcher.summarize();
        for (const auto &hit : matcher.hits()) {
            if (!hit.reported() || !hit.included())
                continue;

            const std::string name = hit.name();
            ++hit_count;
            // size_t id = std::stoull(name);  // Slash and everything after is ignored automatically
        }
        INFO("Hits: " << hit_count << " over " << seq.size() - k + 1);
    }
    INFO("Best scores: " << best_scores);
    }

    return 0;
}

// vim: set ts=4 sw=4 et :
