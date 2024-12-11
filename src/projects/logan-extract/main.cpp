//***************************************************************************
//* Copyright (c) 2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/construction_helper.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "assembly_graph/construction/debruijn_graph_constructor.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "kmer_index/extension_index/kmer_extension_index_builder.hpp"

#include "toolchain/utils.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "io/reads/io_helper.hpp"
#include "io/reads/osequencestream.hpp"
#include "io/utils/id_mapper.hpp"
#include "io/graph/gfa_writer.hpp"
#include "utils/filesystem/temporary.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/parallel/openmp_wrapper.h"

#include <clipp/clipp.h>
#include <filesystem>
#include <iostream>
#include <string>

using namespace debruijn_graph;

enum class output_type {
    unitigs, gfa
};

struct gcfg {
    unsigned k = 29;
    std::filesystem::path tmpdir;
    std::filesystem::path unitigs;
    std::filesystem::path outfile;
    std::filesystem::path list;
    unsigned nthreads;
    unsigned dist = 5000;
    size_t buff_size = 512ULL << 20;
    enum output_type mode = output_type::gfa;

    gcfg()
        : nthreads(omp_get_max_threads() / 2 + 1)
    {}
};


static void process_cmdline(int argc, char** argv, gcfg& cfg) {
    using namespace clipp;

    std::string outfile;
    std::string unitigs;
    std::string tmpdir;
    std::string list;

    auto cli = (
        unitigs << value("Logan unitigs (in FASTA)"),
        list << value("File with list of edges to extract"),
        outfile << value("output file"),
        (option("-d") & integer("value", cfg.dist)) % "neighbourhood size",
        (option("-k") & integer("value", cfg.k)) % "k-mer length to use",
        (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
        (option("-tmp-dir") & value("dir", tmpdir)) % "scratch directory to use",
        (option("-b") & integer("value", cfg.buff_size)) % "sorting buffer size, per thread",
        one_of(option("--unitigs").set(cfg.mode, output_type::unitigs) % "produce unitigs",
               option("--gfa").set(cfg.mode, output_type::gfa) % "produce graph in GFA1 format (default)")
    );

    auto result = parse(argc, argv, cli);
    if (!result) {
        std::cout << make_man_page(cli, argv[0]);
        exit(1);
    }

    cfg.unitigs = unitigs;
    cfg.outfile = outfile;
    cfg.tmpdir = tmpdir.empty() ? "tmp" : tmpdir;
    cfg.list = list;
}

typedef std::vector<std::tuple<EdgeId, std::string, bool>> Links;

static void HandleLink(Links &links,
                       EdgeId e1,
                       std::string_view link,
                       const ConjugateDeBruijnGraph &g) {
    VERIFY_MSG(utils::starts_with(link, "L"), "invalid link record: " << link);
    auto record = utils::split(link, ":");
    VERIFY_MSG(record.size() == 4, "invalid link record: " << record);

    if (record[1] == "-")
        e1 = g.conjugate(e1);

    links.emplace_back(e1, record[2], record[3] == "-");
}

static void HandleSequence(const io::SingleRead &read,
                           Links &links,
                           io::IdMapper<std::string> &mapper,
                           ConjugateDeBruijnGraph &g,
                           ConjugateDeBruijnGraph::HelperT &helper) {
    auto meta = utils::split(read.comment(), " ", /* compress */ true);
    VERIFY_MSG(meta.size() >= 1, "invalid comment string for read: " << read);

    unsigned cov = 0;
    {
        VERIFY_MSG(utils::starts_with(meta.front(), "ka"), "invalid `ka` tag: " << meta.front());
        auto kaTag = utils::split(meta.front(), ":");
        VERIFY_MSG(kaTag.size() == 3, "invalid `ka` tag: " << meta.front());
        float ka = std::stof(std::string(kaTag.back()));
        cov = unsigned(ka * float(read.size()));
    }

    EdgeId e = helper.AddEdge(DeBruijnEdgeData(Sequence{read.GetSequenceString()}));

    g.coverage_index().SetRawCoverage(e, cov);
    g.coverage_index().SetRawCoverage(g.conjugate(e), cov);

    // Logan contigs format: ACCESSION_ID, we need to chop out the accession
    auto accAndName = utils::split(read.name(), "_");
    VERIFY_MSG(accAndName.size() > 0, "unexpected name format: " << read.name());
    std::string name(accAndName.back());

    DEBUG("Map ids: " << e.int_id() << ":" << name);
    mapper.map(name, e.int_id());

    EdgeId ce = g.conjugate(e);
    if (e != ce) {
        DEBUG("Map ids: " << ce.int_id() << ":" << name << "'");
        mapper.map(name + '\'', ce.int_id());
    }
    std::vector<LinkId> empty_links;
    VertexId v1 = helper.CreateVertex(DeBruijnVertexData(empty_links));
    helper.LinkIncomingEdge(v1, e);

    if (e != ce) {
        VertexId v2 = helper.CreateVertex(DeBruijnVertexData(empty_links));
        helper.LinkIncomingEdge(v2, ce);
    }

    auto first_link = std::next(meta.begin());
    while (first_link != meta.end())
        HandleLink(links, e, *first_link++, g);
}

static void ProcessLinks(DeBruijnGraph &g,
                         ConjugateDeBruijnGraph::HelperT &helper,
                         const io::IdMapper<std::string> &mapper,
                         const Links &links) {
    for (const auto &link : links) {
        EdgeId e1 = std::get<0>(link), e2 = mapper[std::get<1>(link)];
        if (std::get<2>(link))
            e2 = g.conjugate(e2);
        VertexId v1 = g.EdgeEnd(e1);
        VertexId v2 = g.EdgeStart(e2);

        g.set_overlap(v1, 0);
        g.set_overlap(v2, 0);
        g.set_overlap(g.conjugate(v1), 0);
        g.set_overlap(g.conjugate(v2), 0);

        helper.LinkEdges(e1, e2);
    }
}

int main(int argc, char** argv) {
  utils::segfault_handler sh;
  gcfg cfg;

  process_cmdline(argc, argv, cfg);

  toolchain::create_console_logger();

  START_BANNER("Logan graph neighborhood extractor");

  unsigned nthreads = cfg.nthreads;
  unsigned k = cfg.k;
  std::filesystem::path tmpdir = cfg.tmpdir;
  size_t buff_size = cfg.buff_size;

  if (k < runtime_k::MIN_K)
      FATAL_ERROR("k-mer size " << k << " is too low");
  if (k >= runtime_k::MAX_K)
      FATAL_ERROR("k-mer size " << k << " is too high, recompile with larger SPADES_MAX_K option");
  if (k % 2 == 0)
      FATAL_ERROR("k-mer size must be odd");


  INFO("K-mer length set to " << k);
  switch (cfg.mode) {
      case output_type::unitigs:
          INFO("Producing unitigs only");
          break;
      case output_type::gfa:
          INFO("Producing graph in GFA1 format");
          break;
  }

  nthreads = spades_set_omp_threads(nthreads);
  INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << nthreads);

  try {
      auto id_mapper = std::make_unique<io::IdMapper<std::string>>();
      auto graph = std::make_unique<ConjugateDeBruijnGraph>(0);

      {
          auto helper = graph->GetConstructionHelper();
          auto unitigs = io::EasyStream(cfg.unitigs,
                                        false, false,
                                        io::FileReadFlags::names_and_comments());
          Links links;
          while (!unitigs.eof()) {
              io::SingleRead read;
              unitigs >> read;
              HandleSequence(read, links, *id_mapper, *graph, helper);
          }

          INFO("Processing links");
          ProcessLinks(*graph, helper, *id_mapper, links);

          INFO("Finalizing the graph");
          // Add "point tips" of edges
          for (EdgeId e: graph->edges()) {
              if (graph->EdgeEnd(e))
                  continue;

              helper.LinkIncomingEdge(helper.CreateVertex(DeBruijnVertexData(0)), e);
          }

          // "Filtering dangling vertices"
          for (VertexId v : graph->vertices()) {
              if (graph->OutgoingEdgeCount(v) > 0 || graph->IncomingEdgeCount(v) > 0)
                  continue;
              graph->DeleteVertex(v);
          }
      }

      INFO("Graph loaded. Total vertices: " << graph->size() << ", total edges: " << graph->e_size());

      std::vector<VertexId> vertices;
      {
          std::ifstream list(cfg.list);
          VERIFY_MSG(list.is_open(), "failed to open file: " << cfg.list);
          for (std::string line; std::getline(list, line); ) {
              auto entries = utils::split(line, "_");
              // Expecting diamond format, e.g. DRR472539_15_ka_f_1875.342_L_8_L_24_4_1_177_
              VERIFY_MSG(entries.size() > 2, "invalid list line: " << line);

              EdgeId e = (*id_mapper)[std::string(entries[1])];
              vertices.push_back(graph->EdgeStart(e));
              vertices.push_back(graph->EdgeEnd(e));
          }
      }

      INFO("Total vertices to explore: " << vertices.size());

      std::unordered_set<VertexId> reached;
      {
          std::vector<std::unordered_set<VertexId>> reached_by_thread(omp_get_max_threads());
#         pragma omp parallel for schedule(guided)
          for (size_t i = 0; i < vertices.size(); ++i) {
              auto dijkstra = omnigraph::CreateUnorientedEdgeBoundedDijkstra(*graph, cfg.dist);
              dijkstra.Run(vertices[i]);
              for (const auto &entry : dijkstra.reached())
                  reached_by_thread[omp_get_thread_num()].insert(entry.first);
          }

          for (const auto &entry : reached_by_thread)
              reached.insert(entry.begin(), entry.end());
      }

      auto component = omnigraph::GraphComponent<ConjugateDeBruijnGraph>::FromVertices(*graph,
                                                                                       reached.begin(), reached.end(),
                                                                                       true);
      INFO("Neighbourhood vertices: " << component.v_size() << ", edges: " << component.e_size());

      if (0) {
          std::ofstream f(cfg.outfile.string() + ".initial");
          gfa::GFAWriter gfa_writer(*graph, f);
          gfa_writer.WriteSegmentsAndLinks(component);
      }

      create_directory(tmpdir);
      auto workdir = fs::tmp::make_temp_dir(tmpdir, "construction");

      auto edges_dir = tmpdir / "edges";
      create_directory(edges_dir);
      io::ReadConverter::ConvertComponentEdgeSequencesToBinary(component, edges_dir, nthreads);

      io::SequencingLibraryT seq_lib;
      seq_lib.set_type(io::LibraryType::TrustedContigs);
      seq_lib.set_orientation(io::LibraryOrientation::Undefined);
      seq_lib.data().lib_index = size_t(-1);
      auto& bin_info = seq_lib.data().binary_reads_info;
      bin_info.single_read_prefix = edges_dir / "edges";
      bin_info.bin_reads_info_file = edges_dir / "edges_info";
      bin_info.binary_converted = true;
      bin_info.chunk_num = nthreads;

      auto read_streams = io::single_binary_readers(seq_lib, true, false);

      // Step 1: build extension index
      VERIFY_MSG(read_streams.size(), "No input streams specified");
      kmers::DeBruijnExtensionIndex<> ext_index(k);

      auto kmers = kmers::DeBruijnExtensionIndexBuilder().BuildExtensionIndexFromStream(workdir, ext_index,
                                                                                        read_streams, buff_size);
      // Step 2: extract unbranching paths
      bool keep_perfect_loops = true;
      std::vector<Sequence> edge_sequences;
      unsigned nchunks = 16 * omp_get_max_threads();
      if (keep_perfect_loops)
          edge_sequences = debruijn_graph::UnbranchingPathExtractor(ext_index, k).ExtractUnbranchingPathsAndLoops(nchunks);
      else
          edge_sequences = debruijn_graph::UnbranchingPathExtractor(ext_index, k).ExtractUnbranchingPaths(nchunks);

      if (cfg.mode == output_type::unitigs) {
          // Step 3: output stuff
          INFO("Saving unitigs to " << cfg.outfile);
          size_t idx = 1;
          std::ofstream f(cfg.outfile);
          for (const auto &edge: edge_sequences) {
              f << std::string(">") << io::MakeContigId(idx++, edge.size(), "EDGE") << std::endl;
              io::WriteWrapped(edge.str(), f);
          }
      } else if (cfg.mode == output_type::gfa) {
          // Step 3: build the graph
          INFO("Building graph");
          debruijn_graph::DeBruijnGraph g(k);
          debruijn_graph::FastGraphFromSequencesConstructor<debruijn_graph::DeBruijnGraph>(k, ext_index).ConstructGraph(g, edge_sequences);

          INFO("Saving graph to " << cfg.outfile);
          std::ofstream f(cfg.outfile);
          gfa::GFAWriter gfa_writer(g, f);
          gfa_writer.WriteSegmentsAndLinks();
      } else
          FATAL_ERROR("Invalid mode");
  } catch (const std::string& s) {
      std::cerr << s << std::endl;
      return EINTR;
  } catch (const std::exception& e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return EINTR;
  }

  return 0;
}
