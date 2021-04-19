//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

#include "io/graph/gfa_writer.hpp"
#include "toolchain/utils.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/parallel/openmp_wrapper.h"

#include "adt/concurrent_dsu.hpp"

#include <clipp/clipp.h>
#include <filesystem>
#include <iostream>
#include <string>

using namespace debruijn_graph;

struct gcfg {
    std::filesystem::path graph;
    std::filesystem::path output_base;
};


static void process_cmdline(int argc, char** argv, gcfg& cfg) {
    using namespace clipp;

    std::string output_base;
    std::string graph;
  
  auto cli = (
      graph << value("graph (in GFA)"),
      output_base << value("output base")
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
    std::cout << make_man_page(cli, argv[0]);
    exit(1);
  }

  cfg.graph = graph;
  cfg.output_base = output_base;
}


// Finds weakly connected components of the given graph
class WeaklyConnectedComponentsFinder {
    const Graph &graph_;

  public:
    explicit WeaklyConnectedComponentsFinder(const Graph &graph)
            : graph_(graph) {};

    WeaklyConnectedComponentsFinder() = delete;

    size_t Run(dsu::ConcurrentDSU &components, const size_t max_id) {
        const unsigned nthreads = omp_get_max_threads();

        omnigraph::IterationHelper<Graph, EdgeId> edges(graph_);
        const auto &ranges = edges.Ranges(nthreads);

#pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < ranges.size(); ++i) {
            for (EdgeId e : ranges[i]) {
                VertexId v = graph_.EdgeStart(e);
                VertexId u = graph_.EdgeEnd(e);
                
                components.unite(v.int_id(), u.int_id());
                components.unite(graph_.conjugate(v).int_id(), u.int_id());
                components.unite(graph_.conjugate(u).int_id(), v.int_id());
            }
        }
        const size_t nonexistent_vertices = max_id + 1 - graph_.size();
        return components.num_sets() - nonexistent_vertices;
    }
};

int main(int argc, char** argv) {
  utils::segfault_handler sh;
  gcfg cfg;

  process_cmdline(argc, argv, cfg);

  toolchain::create_console_logger();

  START_BANNER("GFA splitter");

  try {
      std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());
      std::unique_ptr<ConjugateDeBruijnGraph> graph;

      {
          gfa::GFAReader gfa(cfg.graph);
          INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links() << ", paths: " << gfa.num_paths());
          INFO("Detected k: " << gfa.k());
          VERIFY_MSG(gfa.k() != -1U, "Failed to determine k-mer length");
          VERIFY_MSG(gfa.k() == 0 || gfa.k() % 2 == 1, "k-mer length must be odd");

          graph.reset(new ConjugateDeBruijnGraph(gfa.k()));
      
          gfa.to_graph(*graph, id_mapper.get());

          INFO("Graph loaded. Total vertices: " << graph->size() << ", total edges: " << graph->e_size());
      }

      INFO("Looking for weakly connected components");
      std::vector<std::vector<VertexId>> component_to_vertices;
      std::vector<std::vector<gfa::GFAReader::GFAPath>> component_to_paths;

      size_t cnt = 0;

      {
          dsu::ConcurrentDSU weakly_connected_components(graph->max_vid() + 1);
          cnt = WeaklyConnectedComponentsFinder(*graph).Run(weakly_connected_components, graph->max_vid());
          INFO("Done, total components: " << cnt);

          phmap::flat_hash_map<size_t, size_t> root_id_to_component_id;
          component_to_vertices.resize(cnt);
          component_to_paths.resize(cnt);

          INFO("Assigning roots");
          {
              size_t cur_component = 0;
              for (VertexId v : *graph) {
                  size_t root = weakly_connected_components.find_set(v.int_id());
                  if (root_id_to_component_id.count(root))
                      continue;

                  root_id_to_component_id[root] = cur_component;
                  component_to_vertices[cur_component].reserve(weakly_connected_components.set_size(root));
                  ++cur_component;
              }
          }

          INFO("Assigning vertices to roots");
          for (VertexId v : *graph) {
              size_t root = weakly_connected_components.find_set(v.int_id());
              size_t cur_component_id = root_id_to_component_id[root];
              component_to_vertices[cur_component_id].push_back(v);
          }

          INFO("Splitting paths");
          for (auto &path : paths) {
              VertexId v = graph->EdgeEnd(path.edges.front());
              size_t root = weakly_connected_components.find_set(v.int_id());
              size_t cur_component_id = root_id_to_component_id[root];
              component_to_paths[cur_component_id].emplace_back(std::move(path));
          }
      }

      INFO("Writing components");
      create_directories(cfg.output_base);
      size_t total_v = 0;
      for (size_t i = 0; i < component_to_vertices.size(); ++i) {
          const auto& component = component_to_vertices[i];
          total_v += component.size();
          auto subgraph =
                  omnigraph::GraphComponent<ConjugateDeBruijnGraph>::FromVertices(*graph,
                                                                                  component.begin(),  component.end(),
                                                                                  true);
          std::ofstream os(cfg.output_base / ("subgraph_" + std::to_string(i) + ".gfa"));
          path_extend::GFAPathWriter writer(*graph, os,
                                            io::MapNamingF<debruijn_graph::ConjugateDeBruijnGraph>(*id_mapper));
          writer.WriteSegmentsAndLinks(subgraph);
          for (const auto& path : component_to_paths[i])
              writer.WritePaths(path);
          VERBOSE_POWER_T2(i, 0, "Written " << i << " components, total vertices: " << total_v);
      }
      INFO("Written " << component_to_vertices.size() << " components, total vertices: " << total_v);
  } catch (const std::string& s) {
      std::cerr << s << std::endl;
      return EINTR;
  } catch (const std::exception& e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return EINTR;
  }

  return 0;
}
