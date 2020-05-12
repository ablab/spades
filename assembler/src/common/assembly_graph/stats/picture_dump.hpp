//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "statistics.hpp"
#include "assembly_graph/core/graph.hpp"

#include "pipeline/graph_pack.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "io/binary/graph_pack.hpp"
#include "io/binary/paired_index.hpp"
//FIXME awful dependency to get write_lib_data
#include "pipeline/config_struct.hpp"
#include "visualization/position_filler.hpp"
#include "visualization/visualization.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "io/reads/rc_reader_wrapper.hpp"
#include "io/reads/delegating_reader_wrapper.hpp"
#include "io/reads/io_helper.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "io/reads/osequencestream.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "utils/filesystem/copy_file.hpp"
#include "sequence/genome_storage.hpp"
#include "paired_info/paired_info.hpp"

#include <boost/algorithm/string.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

namespace debruijn_graph {

namespace stats {

template<class Graph, class Index>
MappingPath<typename Graph::EdgeId>
FindGenomeMappingPath(const Sequence& genome, const Graph& g,
                      const Index& index,
                      const KmerMapper<Graph>& kmer_mapper) {
    BasicSequenceMapper<Graph, Index> srt(g, index, kmer_mapper);
    return srt.MapSequence(genome);
}

inline MappingPath<Graph::EdgeId> FindGenomeMappingPath(const Sequence& genome, const GraphPack &gp) {
    const auto &graph = gp.get<Graph>();
    const auto &index = gp.get<EdgeIndex<Graph>>();
    const auto &kmer_mapper = gp.get<KmerMapper<Graph>>();
    return FindGenomeMappingPath(genome, graph, index, kmer_mapper);
}

inline std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> DefaultColorer(const GraphPack &gp) {
    const auto &graph = gp.get<Graph>();
    const auto &index = gp.get<EdgeIndex<Graph>>();
    const auto &kmer_mapper = gp.get<KmerMapper<Graph>>();
    const auto &genome = gp.get<GenomeStorage>();
    return visualization::graph_colorer::DefaultColorer(graph,
        FindGenomeMappingPath(genome.GetSequence(), graph, index, kmer_mapper).path(),
        FindGenomeMappingPath(!genome.GetSequence(), graph, index, kmer_mapper).path());
}

inline void CollectContigPositions(GraphPack &gp) {
    if (!cfg::get().pos.contigs_for_threading.empty() &&
        fs::FileExists(cfg::get().pos.contigs_for_threading))
      visualization::position_filler::FillPos(gp, cfg::get().pos.contigs_for_threading, "thr_", true);

    if (!cfg::get().pos.contigs_to_analyze.empty() &&
        fs::FileExists(cfg::get().pos.contigs_to_analyze))
      visualization::position_filler::FillPos(gp, cfg::get().pos.contigs_to_analyze, "anlz_", true);
}

template<class Graph, class Index>
class GenomeMappingStat: public AbstractStatCounter {
  private:
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    const Index& index_;
    Sequence genome_;
    size_t k_;
  public:
    GenomeMappingStat(const Graph &graph, const Index &index, const GenomeStorage &genome, size_t k) :
            graph_(graph), index_(index), genome_(genome.GetSequence()), k_(k) {}

    virtual ~GenomeMappingStat() {}

    virtual void Count() {
        INFO("Mapping genome");
        size_t break_number = 0;
        size_t covered_kp1mers = 0;
        size_t fail = 0;
        if (genome_.size() <= k_)
            return;

        RtSeq cur = genome_.start<RtSeq>(k_ + 1);
        cur >>= 0;
        bool breaked = true;
        std::pair<EdgeId, size_t> cur_position;
        for (size_t cur_nucl = k_; cur_nucl < genome_.size(); cur_nucl++) {
            cur <<= genome_[cur_nucl];
            if (index_.contains(cur)) {
                auto next = index_.get(cur);
                if (!breaked
                    && cur_position.second + 1
                    < graph_.length(cur_position.first)) {
                    if (next.first != cur_position.first
                        || cur_position.second + 1 != next.second) {
                        fail++;
                    }
                }
                cur_position = next;
                covered_kp1mers++;
                breaked = false;
            } else {
                if (!breaked) {
                    breaked = true;
                    break_number++;
                }
            }
        }
        INFO("Genome mapped");
        INFO("Genome mapping results:");
        INFO("Covered k+1-mers:" << covered_kp1mers << " of " << (genome_.size() - k_) << " which is "
             << (100.0 * (double) covered_kp1mers / (double) (genome_.size() - k_)) << "%");
        INFO("Covered k+1-mers form " << break_number + 1 << " contiguous parts");
        INFO("Continuity failtures " << fail);
    }
};

template<class Graph>
void WriteErrorLoc(const Graph &g,
                   const std::string& folder_name,
                   std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> genome_colorer,
                   const visualization::graph_labeler::GraphLabeler<Graph>& labeler) {
    INFO("Writing error localities for graph to folder " << folder_name);
    auto all = GraphComponent<Graph>::WholeGraph(g);
    auto edges = genome_colorer->ColoredWith(all.edges().begin(), all.edges().end(), "black");
    std::set<typename Graph::VertexId> to_draw;
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        to_draw.insert(g.EdgeEnd(*it));
        to_draw.insert(g.EdgeStart(*it));
    }
    auto splitter = StandardSplitter(g, to_draw);
    visualization::visualization_utils::WriteComponents(g, folder_name, splitter, genome_colorer, labeler);
    INFO("Error localities written written to folder " << folder_name);
}

inline void CountStats(const GraphPack &gp) {
    typedef typename Graph::EdgeId EdgeId;
    INFO("Counting stats");
    StatList stats;
    const auto &graph = gp.get<Graph>();
    const auto &index = gp.get<EdgeIndex<Graph>>();
    const auto &kmer_mapper = gp.get<KmerMapper<Graph>>();
    const auto &genome = gp.get<GenomeStorage>();

    Path<EdgeId> path1 = FindGenomeMappingPath(genome.GetSequence(), graph, index, kmer_mapper).path();
    Path<EdgeId> path2 = FindGenomeMappingPath(!genome.GetSequence(), graph, index, kmer_mapper).path();
    stats.AddStat(new VertexEdgeStat<Graph>(graph));
    stats.AddStat(new BlackEdgesStat<Graph>(graph, path1, path2));
    stats.AddStat(new NStat<Graph>(graph, path1, 50));
    stats.AddStat(new SelfComplementStat<Graph>(graph));
    stats.AddStat(new GenomeMappingStat<Graph, EdgeIndex<Graph>>(graph, index, genome, gp.k()));
    stats.AddStat(new IsolatedEdgesStat<Graph>(graph, path1, path2));
    stats.Count();
    INFO("Stats counted");
}

template<class Graph>
void WriteGraphComponentsAlongGenome(const Graph &g,
                                     const visualization::graph_labeler::GraphLabeler<Graph> &labeler,
                                     const std::string &folder,
                                     const Path<typename Graph::EdgeId> &path1,
                                     const Path<typename Graph::EdgeId> &path2) {
    INFO("Writing graph components along genome");

    fs::make_dir(folder);
    visualization::visualization_utils::WriteComponentsAlongPath(g, path1, folder,
                                                                 visualization::graph_colorer::DefaultColorer(g, path1, path2),
                                                                 labeler);

    INFO("Writing graph components along genome finished");
}

//todo refactoring needed: use graph pack instead!!!
template<class Graph, class Mapper>
void WriteGraphComponentsAlongContigs(const Graph &g,
                                      Mapper &mapper,
                                      const std::string &folder,
                                      std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer,
                                      const visualization::graph_labeler::GraphLabeler<Graph> &labeler) {
    INFO("Writing graph components along contigs");
    auto contigs_to_thread = io::EasyStream(cfg::get().pos.contigs_to_analyze, false);
    contigs_to_thread.reset();
    io::SingleRead read;
    while (!contigs_to_thread.eof()) {
        contigs_to_thread >> read;
        fs::make_dir(folder + read.name());
        visualization::visualization_utils::WriteComponentsAlongPath(g, mapper.MapSequence(read.sequence()).simple_path(),
                                                                     folder + read.name() + "/", colorer, labeler);
    }
    INFO("Writing graph components along contigs finished");
}

inline void WriteKmerComponent(const GraphPack &gp, const RtSeq &kp1mer, const std::string &file,
                        std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer,
                        const visualization::graph_labeler::GraphLabeler<Graph> &labeler) {
    const auto &graph = gp.get<Graph>();
    const auto &index = gp.get<EdgeIndex<Graph>>();
    if(!index.contains(kp1mer)) {
        WARN("no such kmer in the graph");
        return;
    }
    VERIFY(index.contains(kp1mer));
    auto pos = index.get(kp1mer);
    typename Graph::VertexId v = pos.second * 2 < graph.length(pos.first)
            ? graph.EdgeStart(pos.first)
            : graph.EdgeEnd(pos.first);
    GraphComponent<Graph> component = omnigraph::VertexNeighborhood<Graph>(graph, v);
    visualization::visualization_utils::WriteComponent<Graph>(component, file, colorer, labeler);
}

inline boost::optional<RtSeq> FindCloseKP1mer(const GraphPack &gp, size_t genome_pos, size_t k) {
    const auto &genome = gp.get<GenomeStorage>();
    const auto &kmer_mapper = gp.get<KmerMapper<Graph>>();
    const auto &index = gp.get<EdgeIndex<Graph>>();
    VERIFY(genome.size() > 0);
    VERIFY(genome_pos < genome.size());
    static const size_t magic_const = 200;
    for (size_t diff = 0; diff < magic_const; diff++) {
        for (int dir = -1; dir <= 1; dir += 2) {
            size_t pos = (genome.size() - k + genome_pos + dir * diff) % (genome.size() - k);
            RtSeq kp1mer = kmer_mapper.Substitute(RtSeq (k + 1, genome.GetSequence(), pos));
            if (index.contains(kp1mer))
                return boost::optional<RtSeq>(kp1mer);
        }
    }
    return boost::none;
}

inline void PrepareForDrawing(GraphPack &gp) {
    gp.EnsureDebugInfo();
    CollectContigPositions(gp);
}


struct detail_info_printer {
    detail_info_printer(GraphPack &gp,
                        const visualization::graph_labeler::GraphLabeler<Graph> &labeler,
                        const std::string &folder)
            :  gp_(gp),
               labeler_(labeler),
               folder_(folder) {
    }

    void operator()(config::info_printer_pos pos, const std::string &folder_suffix = "") {
        auto pos_name = ModeName(pos, config::InfoPrinterPosNames());

        ProduceDetailedInfo(pos_name + folder_suffix, pos);
    }

  private:

    template<typename T>
    std::string ToString(const T& t, size_t length) {
        std::ostringstream ss;
        ss << t;
        std::string result = ss.str();
        while (result.size() < length)
            result = "0" + result;
        return result;
    }


    void ProduceDetailedInfo(const std::string &pos_name,
                             config::info_printer_pos pos) {
        using namespace visualization;
        using namespace io::binary;

        const auto &graph = gp_.get<Graph>();
        const auto &index = gp_.get<EdgeIndex<Graph>>();
        const auto &kmer_mapper = gp_.get<KmerMapper<Graph>>();

        static size_t call_cnt = 0;

        auto it = cfg::get().info_printers.find(pos);
        VERIFY(it != cfg::get().info_printers.end());
    
        const config::debruijn_config::info_printer & config = it->second;
    
        if (config.basic_stats) {
            VertexEdgeStat<Graph> stats(graph);
            INFO("Number of vertices : " << stats.vertices() << ", number of edges : "
                  << stats.edges() << ", sum length of edges : " << stats.edge_length());
        }

        if (config.save_graph_pack) {
            auto saves_folder = fs::append_path(fs::append_path(folder_, "saves/"),
                                                ToString(call_cnt++, 2) + "_" + pos_name + "/");
            fs::make_dirs(saves_folder);
            BasePackIO().Save(saves_folder + "graph_pack", gp_);
            //TODO: separate
            const auto &indices = gp_.get<omnigraph::de::PairedInfoIndicesT<Graph>>("clustered_indices");
            Save(saves_folder + "graph_pack", indices);
        }

        if (config.save_all) {
            auto saves_folder = fs::append_path(fs::append_path(folder_, "saves/"),
                                                ToString(call_cnt++, 2) + "_" + pos_name);
            fs::make_dirs(saves_folder);
            auto p = saves_folder + "/saves";
            INFO("Saving current state to " << p);

            FullPackIO().Save(saves_folder + "graph_pack", gp_);
            debruijn_graph::config::write_lib_data(p);
        }

        if (config.save_full_graph) {
            auto saves_folder = fs::append_path(fs::append_path(folder_, "saves/"),
                                                ToString(call_cnt++, 2) + "_" + pos_name + "/");
            fs::make_dirs(saves_folder);

            BasicGraphIO<Graph>().Save(saves_folder + "graph", graph);
        }

        if (config.lib_info) {
            auto saves_folder = fs::append_path(fs::append_path(folder_, "saves/"),
                                                ToString(call_cnt++, 2) + "_" + pos_name + "/");
            fs::make_dirs(saves_folder);
            config::write_lib_data(saves_folder + "lib_info");
        }

        if (config.extended_stats) {
            VERIFY(cfg::get().developer_mode);
            CountStats(gp_);
        }

        if (!(config.write_error_loc ||
            config.write_full_graph ||
            config.write_full_nc_graph ||
            config.write_components ||
            !config.components_for_kmer.empty() ||
            config.write_components_along_genome ||
            config.write_components_along_contigs ||
            !config.components_for_genome_pos.empty())) {
            return;
        } 

        VERIFY(cfg::get().developer_mode);
        auto pics_folder = fs::append_path(fs::append_path(folder_, "pictures/"),
                                           ToString(call_cnt++, 2) + "_" + pos_name + "/");
        fs::make_dirs(pics_folder);
        PrepareForDrawing(gp_);
    
        auto path1 = FindGenomeMappingPath(gp_.get<GenomeStorage>().GetSequence(), graph, index, kmer_mapper).path();
    
        auto colorer = DefaultColorer(gp_);
    
        if (config.write_error_loc) {
            fs::make_dir(pics_folder + "error_loc/");
            WriteErrorLoc(graph, pics_folder + "error_loc/", colorer, labeler_);
        }
    
        if (config.write_full_graph) {
            visualization_utils::WriteComponent(GraphComponent<Graph>::WholeGraph(graph),
                                                pics_folder + "full_graph.dot", colorer, labeler_);
        }
    
        if (config.write_full_nc_graph) {
            visualization_utils::WriteSimpleComponent(GraphComponent<Graph>::WholeGraph(graph),
                                                      pics_folder + "nc_full_graph.dot", colorer, labeler_);
        }
    
        if (config.write_components) {
            fs::make_dir(pics_folder + "components/");
            visualization_utils::WriteComponents(graph, pics_folder + "components/",
                                                 omnigraph::ReliableSplitter<Graph>(graph), colorer, labeler_);
        }
    
        if (!config.components_for_kmer.empty()) {
            auto kmer_folder = fs::append_path(pics_folder, "kmer_loc/");
            fs::make_dir(kmer_folder);
            auto kmer = RtSeq(gp_.k() + 1, config.components_for_kmer.substr(0, gp_.k() + 1).c_str());
            auto file_name = fs::append_path(kmer_folder, pos_name + ".dot");
            WriteKmerComponent(gp_, kmer, file_name, colorer, labeler_);
        }
    
        if (config.write_components_along_genome) {
            fs::make_dir(pics_folder + "along_genome/");
            visualization_utils::WriteComponentsAlongPath
                    (graph, path1.sequence(), pics_folder + "along_genome/", colorer, labeler_);
        }
    
        if (config.write_components_along_contigs) {
            fs::make_dir(pics_folder + "along_contigs/");
            BasicSequenceMapper<Graph, EdgeIndex<Graph>> mapper(graph, index, kmer_mapper);
            WriteGraphComponentsAlongContigs(graph, mapper, pics_folder + "along_contigs/", colorer, labeler_);
        }

        if (!config.components_for_genome_pos.empty()) {
            auto pos_loc_folder = fs::append_path(pics_folder, "pos_loc/");
            fs::make_dir(pos_loc_folder);
            std::vector<std::string> positions;
            boost::split(positions, config.components_for_genome_pos,
                         boost::is_any_of(" ,"), boost::token_compress_on);
            const auto &genome = gp_.get<GenomeStorage>();
            for (auto it = positions.begin(); it != positions.end(); ++it) {
                auto close_kp1mer = FindCloseKP1mer(gp_, std::stoi(*it), gp_.k());
                if (close_kp1mer) {
                    auto locality_folder = fs::append_path(pos_loc_folder, *it + "/");
                    fs::make_dir(locality_folder);
                    WriteKmerComponent(gp_, *close_kp1mer, fs::append_path(locality_folder, pos_name + ".dot"), colorer, labeler_);
                } else {
                    WARN("Failed to find genome kp1mer close to the one at position "
                        << *it << " in the graph. Which is "
                        << RtSeq (gp_.k() + 1, genome.GetSequence(), std::stoi(*it)));
                }
            }
        }
    }

    GraphPack& gp_;
    const visualization::graph_labeler::GraphLabeler<Graph>& labeler_;
    std::string folder_;
};

inline
std::string ConstructComponentName(const std::string &file_name, size_t cnt) {
    return file_name + std::to_string(cnt);
}


}
}
