//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/components/graph_component.hpp"

#include "visualization/visualization.hpp"
#include "pipeline/graphio.hpp"
#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"
#include "io/reads/io_helper.hpp"

#include "version.hpp"

#include "hmmfile.hpp"
#include "hmmmatcher.hpp"
#include "fees.hpp"
#include "omnigraph_wrapper.hpp"

#include <clipp/clipp.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <string>

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

struct cfg {
    std::string load_from;
    std::string hmmfile;
    size_t k;
    uint64_t int_id;

    hmmer::hmmer_cfg hcfg;
    cfg()
            : load_from(""), hmmfile(""), k(0), int_id(0)
    {}
};

extern "C" {
#include "p7_config.h"

#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
}

void process_cmdline(int argc, char **argv, cfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.hmmfile    << value("hmm file"),
      cfg.load_from  << value("load from"),
      cfg.k          << value("k-mer size"),
      (option("--edge_id") & number("value", cfg.int_id)) % "match around edge",
      // Control of output
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
  );

  if (!parse(argc, argv, cli)) {
    std::cout << make_man_page(cli, argv[0]);
    exit(1);
  }
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

template<class GraphCursor>
std::vector<debruijn_graph::EdgeId> to_path(const std::vector<GraphCursor> &cpath) {
    std::vector<debruijn_graph::EdgeId> path;

    auto it = cpath.begin();
    while (it->is_empty())
        ++it;

    for (; it != cpath.end(); ++it) {
        if (it->is_empty())
            continue;

        if (path.size() == 0 || it->edge() != path.back())
            path.push_back(it->edge());
    }

    return path;
}


int main(int argc, char* argv[]) {
    utils::perf_counter pc;
    int textw = 120;

    srand(42);
    srandom(42);

    cfg cfg;
    process_cmdline(argc, argv, cfg);

    create_console_logger();
    INFO("Starting Graph HMM aligning engine, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

    /* Open the query profile HMM file */
    hmmer::HMMFile hmmfile(cfg.hmmfile);
    if (!hmmfile.valid())
        FATAL_ERROR("Error reading HMM file "<< cfg.hmmfile);

    using namespace debruijn_graph;
    ConjugateDeBruijnGraph graph(cfg.k);
    graphio::ScanBasicGraph(cfg.load_from, graph);
    INFO("Graph loaded. Total vertices: " << graph.size());

    // Collect all the edges
    std::vector<EdgeId> edges;
    for (auto it = graph.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        EdgeId edge = *it;
        if (cfg.int_id == 0 ||
            edge.int_id() == cfg.int_id)
        edges.push_back(edge);
    }

    auto hmmw = hmmfile.read();
    ESL_STOPWATCH *w = esl_stopwatch_Create();

    // Outer loop: over each query HMM in <hmmfile>.
    while (hmmw) {
        P7_HMM *hmm = hmmw->get();

        hmmer::HMMMatcher matcher(hmmw.get(), cfg.hcfg);

        if (fprintf(stderr, "Query:       %s  [M=%d]\n", hmm->name, hmm->M) < 0) FATAL_ERROR("write failed");
        if (hmm->acc)  { if (fprintf(stderr, "Accession:   %s\n", hmm->acc)  < 0) FATAL_ERROR("write failed"); }
        if (hmm->desc) { if (fprintf(stderr, "Description: %s\n", hmm->desc) < 0) FATAL_ERROR("write failed"); }

        esl_stopwatch_Start(w);
        for (size_t i = 0; i < edges.size(); ++i) {
            // FIXME: this conversion is pointless
            std::string ref = std::to_string(i);
            std::string seq = graph.EdgeNucls(edges[i]).str();
            //  INFO("EdgeId: " << edge << ", length: " << graph.length(edge) << ", seq: " << graph.EdgeNucls(edge));
            matcher.match(ref.c_str(), seq.c_str());
        }

        matcher.summarize();
        esl_stopwatch_Stop(w);

        auto th = matcher.hits();
        std::vector<EdgeId> match_edges;
        for (size_t h = 0; h < th->N; h++) {
            if (!(th->hit[h]->flags & p7_IS_REPORTED))
                continue;
            if (!(th->hit[h]->flags & p7_IS_INCLUDED))
                continue;

            match_edges.push_back(edges[std::stoull(th->hit[h]->name)]);
        }
        INFO("Total matched edges: " << match_edges.size());

        // Collect the neighbourhood of the matched edges
        auto fdijkstra = omnigraph::CreateBoundedDijkstra(graph, 2 * hmm->M);
        auto bdijkstra = omnigraph::CreateBackwardBoundedDijkstra(graph, 2 * hmm->M);

        for (EdgeId e : match_edges) {
            INFO("Extracting neighbourhood of edge " << e);
            fdijkstra.Run(graph.EdgeEnd(e));
            bdijkstra.Run(graph.EdgeStart(e));
            std::vector<VertexId> vertices = fdijkstra.ReachedVertices();
            std::vector<VertexId> bvertices = bdijkstra.ReachedVertices();

            vertices.insert(vertices.end(), bvertices.begin(), bvertices.end());
            bvertices.clear();

            auto component = omnigraph::GraphComponent<ConjugateDeBruijnGraph>::FromVertices(graph,
                                                                                             vertices.begin(), vertices.end(),
                                                                                             true);
            INFO("Neighbourhood vertices: " << component.v_size() << ", edges: " << component.e_size());

            if (1) {
                INFO("Writing component around edge " << e);
                DrawComponent(component, graph, std::to_string(graph.int_id(e)), match_edges);
            }

            auto fees = hmm::fees_from_hmm(hmm, hmmw->abc());

            auto initial = all(component);
            auto result = find_best_path(fees, initial);

            INFO("Best score: " << result.best_score());
            INFO("Best of the best");
            INFO(result.best_path_string());
            INFO("Extracting top paths");
            auto top_paths = result.top_k(5);
            size_t idx = 0;
            for (const auto& kv : top_paths) {
                // INFO("" << kv.second << ":" << top_paths.str(kv.first));
                auto path = to_path(kv.first);
                INFO("Path length : " << path.size() << " edges");
                for (EdgeId e : path)
                    INFO("" << e.int_id());
                DrawComponent(component, graph, std::to_string(graph.int_id(e)) + "_" + std::to_string(idx), path);
                idx += 1;
            }
        }

        if (0) {
            p7_tophits_Targets(stderr, matcher.hits(), matcher.pipeline(), textw); if (fprintf(stderr, "\n\n") < 0) FATAL_ERROR("write failed");
            p7_tophits_Domains(stderr, matcher.hits(), matcher.pipeline(), textw); if (fprintf(stderr, "\n\n") < 0) FATAL_ERROR("write failed");
            p7_pli_Statistics(stderr, matcher.pipeline(), w); if (fprintf(stderr, "//\n") < 0) FATAL_ERROR("write failed");
        }

        hmmw = hmmfile.read();
    } // end outer loop over query HMMs

    esl_stopwatch_Destroy(w);

    return 0;
}
