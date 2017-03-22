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
#include "pipeline/graphio.hpp"
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

template<class graph_pack>
MappingPath<typename graph_pack::graph_t::EdgeId>
FindGenomeMappingPath(const Sequence& genome, const graph_pack& gp) {
    return FindGenomeMappingPath(genome, gp.g, gp.index, gp.kmer_mapper);
}

template <class graph_pack>
shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> DefaultColorer(const graph_pack& gp) {
    return visualization::graph_colorer::DefaultColorer(gp.g,
        FindGenomeMappingPath(gp.genome.GetSequence(), gp.g, gp.index, gp.kmer_mapper).path(),
        FindGenomeMappingPath(!gp.genome.GetSequence(), gp.g, gp.index, gp.kmer_mapper).path());
}

template <class graph_pack>
void CollectContigPositions(graph_pack &gp) {
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
    GenomeMappingStat(const Graph &graph, const Index &index, GenomeStorage genome, size_t k) :
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
        pair<EdgeId, size_t> cur_position;
        for (size_t cur_nucl = k_; cur_nucl < genome_.size(); cur_nucl++) {
            cur <<= genome_[cur_nucl];
            if (index_.contains(cur)) {
                pair<EdgeId, size_t> next = index_.get(cur);
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
        INFO("Covered k+1-mers form " << break_number + 1 << " contigious parts");
        INFO("Continuity failtures " << fail);
    }
};

template<class Graph>
void WriteErrorLoc(const Graph &g,
                   const string& folder_name,
                   std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> genome_colorer,
                   const visualization::graph_labeler::GraphLabeler<Graph>& labeler) {
    INFO("Writing error localities for graph to folder " << folder_name);
    auto all = GraphComponent<Graph>::WholeGraph(g);
    set<typename Graph::EdgeId> edges = genome_colorer->ColoredWith(all.edges().begin(),
                                                    all.edges().end(), "black");
    set<typename Graph::VertexId> to_draw;
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        to_draw.insert(g.EdgeEnd(*it));
        to_draw.insert(g.EdgeStart(*it));
    }
    shared_ptr<GraphSplitter<Graph>> splitter = StandardSplitter(g, to_draw);
    visualization::visualization_utils::WriteComponents(g, folder_name, splitter, genome_colorer, labeler);
    INFO("Error localities written written to folder " << folder_name);
}

template<class graph_pack>
void CountStats(const graph_pack& gp) {
    typedef typename graph_pack::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;
    INFO("Counting stats");
    StatList stats;
    Path<EdgeId> path1 = FindGenomeMappingPath(gp.genome.GetSequence(), gp.g, gp.index,
                                      gp.kmer_mapper).path();
    Path<EdgeId> path2 = FindGenomeMappingPath(!gp.genome.GetSequence(), gp.g, gp.index,
                                      gp.kmer_mapper).path();
    stats.AddStat(new VertexEdgeStat<Graph>(gp.g));
    stats.AddStat(new BlackEdgesStat<Graph>(gp.g, path1, path2));
    stats.AddStat(new NStat<Graph>(gp.g, path1, 50));
    stats.AddStat(new SelfComplementStat<Graph>(gp.g));
    stats.AddStat(
        new GenomeMappingStat<Graph, Index>(gp.g, gp.index,
                                            gp.genome, gp.k_value));
    stats.AddStat(new IsolatedEdgesStat<Graph>(gp.g, path1, path2));
    stats.Count();
    INFO("Stats counted");
}

template<class Graph>
void WriteGraphComponentsAlongGenome(const Graph& g,
                                     const visualization::graph_labeler::GraphLabeler<Graph>& labeler,
                                     const string& folder,
                                     const Path<typename Graph::EdgeId>& path1,
                                     const Path<typename Graph::EdgeId>& path2) {
    INFO("Writing graph components along genome");

    make_dir(folder);
    visualization::visualization_utils::WriteComponentsAlongPath(g, path1, folder,
                                                                 visualization::graph_colorer::DefaultColorer(g, path1, path2),
                                                                 labeler);

    INFO("Writing graph components along genome finished");
}

//todo refactoring needed: use graph pack instead!!!
template<class Graph, class Mapper>
void WriteGraphComponentsAlongContigs(const Graph& g,
                                      Mapper &mapper,
                                      const std::string& folder,
                                      std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer,
                                      const visualization::graph_labeler::GraphLabeler<Graph>& labeler) {
    INFO("Writing graph components along contigs");
    auto contigs_to_thread = io::EasyStream(cfg::get().pos.contigs_to_analyze, false);
    contigs_to_thread->reset();
    io::SingleRead read;
    while (!contigs_to_thread->eof()) {
        (*contigs_to_thread) >> read;
        make_dir(folder + read.name());
        visualization::visualization_utils::WriteComponentsAlongPath(g, mapper.MapSequence(read.sequence()).simple_path(),
                                                                     folder + read.name() + "/", colorer, labeler);
    }
    INFO("Writing graph components along contigs finished");
}

template<class Graph>
void WriteKmerComponent(conj_graph_pack &gp, RtSeq const& kp1mer, const std::string& file,
                        std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer,
                        const visualization::graph_labeler::GraphLabeler<Graph>& labeler) {
    if(!gp.index.contains(kp1mer)) {
        WARN("no such kmer in the graph");
        return;
    }
    VERIFY(gp.index.contains(kp1mer));
    auto pos = gp.index.get(kp1mer);
    typename Graph::VertexId v = pos.second * 2 < gp.g.length(pos.first) ? gp.g.EdgeStart(pos.first) : gp.g.EdgeEnd(pos.first);
    GraphComponent<Graph> component = omnigraph::VertexNeighborhood<Graph>(gp.g, v);
    visualization::visualization_utils::WriteComponent<Graph>(component, file, colorer, labeler);
}

inline
optional<RtSeq> FindCloseKP1mer(const conj_graph_pack &gp,
                                           size_t genome_pos, size_t k) {
    VERIFY(gp.genome.size() > 0);
    VERIFY(genome_pos < gp.genome.size());
    static const size_t magic_const = 200;
    for (size_t diff = 0; diff < magic_const; diff++) {
        for (int dir = -1; dir <= 1; dir += 2) {
            size_t pos = (gp.genome.size() - k + genome_pos + dir * diff) % (gp.genome.size() - k);
            RtSeq kp1mer = gp.kmer_mapper.Substitute(
                RtSeq (k + 1, gp.genome.GetSequence(), pos));
            if (gp.index.contains(kp1mer))
                return optional<RtSeq>(kp1mer);
        }
    }
    return boost::none;
}

inline
void PrepareForDrawing(conj_graph_pack &gp) {
    gp.EnsureDebugInfo();
    CollectContigPositions(gp);
}


struct detail_info_printer {
    detail_info_printer(conj_graph_pack &gp,
                        const visualization::graph_labeler::GraphLabeler<Graph>& labeler,
                        const string& folder)
            :  gp_(gp),
               labeler_(labeler),
               folder_(folder) {
    }

    void operator() (config::info_printer_pos pos,
                    const string& folder_suffix = "") {
        string pos_name = ModeName(pos, config::InfoPrinterPosNames());

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


    void ProduceDetailedInfo(const string &pos_name,
                             config::info_printer_pos pos) {
        using namespace visualization;

        static size_t call_cnt = 0;

        auto it = cfg::get().info_printers.find(pos);
        VERIFY(it != cfg::get().info_printers.end());
    
        const config::debruijn_config::info_printer & config = it->second;
    
        if (config.basic_stats) {
            VertexEdgeStat<conj_graph_pack::graph_t> stats(gp_.g);
            INFO("Number of vertices : " << stats.vertices() << ", number of edges : "
                  << stats.edges() << ", sum length of edges : " << stats.edge_length());
        }

        if (config.save_graph_pack) {
            string saves_folder = fs::append_path(fs::append_path(folder_, "saves/"),
                                              ToString(call_cnt++, 2) + "_" + pos_name + "/");
            fs::make_dirs(saves_folder);
            graphio::ConjugateDataPrinter<conj_graph_pack::graph_t> printer(gp_.g);
            graphio::PrintGraphPack(saves_folder + "graph_pack", printer, gp_);
            //TODO: separate
            graphio::PrintClusteredIndices(saves_folder + "graph_pack", printer, gp_.clustered_indices);
        }

        if (config.save_all) {
            string saves_folder = fs::append_path(fs::append_path(folder_, "saves/"),
                                                          ToString(call_cnt++, 2) + "_" + pos_name);
            fs::make_dirs(saves_folder);
            string p = saves_folder + "/saves";
            INFO("Saving current state to " << p);

            debruijn_graph::graphio::PrintAll(p, gp_);
            debruijn_graph::config::write_lib_data(p);
        }

        if (config.save_full_graph) {
            string saves_folder = fs::append_path(fs::append_path(folder_, "saves/"),
                                              ToString(call_cnt++, 2) + "_" + pos_name + "/");
            fs::make_dirs(saves_folder);
            graphio::ConjugateDataPrinter<conj_graph_pack::graph_t> printer(gp_.g);
            graphio::PrintBasicGraph(saves_folder + "graph", printer);
        }

        if (config.lib_info) {
            string saves_folder = fs::append_path(fs::append_path(folder_, "saves/"),
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
        string pics_folder = fs::append_path(fs::append_path(folder_, "pictures/"),
                                          ToString(call_cnt++, 2) + "_" + pos_name + "/");
        fs::make_dirs(pics_folder);
        PrepareForDrawing(gp_);
    
        auto path1 = FindGenomeMappingPath(gp_.genome.GetSequence(), gp_.g, gp_.index,
                                           gp_.kmer_mapper).path();
    
        auto colorer = DefaultColorer(gp_);
    
        if (config.write_error_loc) {
            make_dir(pics_folder + "error_loc/");
            WriteErrorLoc(gp_.g, pics_folder + "error_loc/", colorer, labeler_);
        }
    
        if (config.write_full_graph) {
            visualization_utils::WriteComponent(GraphComponent<Graph>::WholeGraph(gp_.g),
                                                pics_folder + "full_graph.dot", colorer, labeler_);
        }
    
        if (config.write_full_nc_graph) {
            visualization_utils::WriteSimpleComponent(GraphComponent<Graph>::WholeGraph(gp_.g),
                                                      pics_folder + "nc_full_graph.dot", colorer, labeler_);
        }
    
        if (config.write_components) {
            make_dir(pics_folder + "components/");
            visualization_utils::WriteComponents(gp_.g, pics_folder + "components/",
                                                 omnigraph::ReliableSplitter<Graph>(gp_.g), colorer, labeler_);
        }
    
        if (!config.components_for_kmer.empty()) {
            string kmer_folder = fs::append_path(pics_folder, "kmer_loc/");
            make_dir(kmer_folder);
            auto kmer = RtSeq(gp_.k_value + 1, config.components_for_kmer.substr(0, gp_.k_value + 1).c_str());
            string file_name = fs::append_path(kmer_folder, pos_name + ".dot");
            WriteKmerComponent(gp_, kmer, file_name, colorer, labeler_);
        }
    
        if (config.write_components_along_genome) {
            make_dir(pics_folder + "along_genome/");
            visualization_utils::WriteComponentsAlongPath
                    (gp_.g, path1.sequence(), pics_folder + "along_genome/", colorer, labeler_);
        }
    
        if (config.write_components_along_contigs) {
            make_dir(pics_folder + "along_contigs/");
            BasicSequenceMapper<Graph, Index> mapper(gp_.g, gp_.index, gp_.kmer_mapper);
            WriteGraphComponentsAlongContigs(gp_.g, mapper, pics_folder + "along_contigs/", colorer, labeler_);
        }

        if (!config.components_for_genome_pos.empty()) {
            string pos_loc_folder = fs::append_path(pics_folder, "pos_loc/");
            make_dir(pos_loc_folder);
            vector<string> positions;
            boost::split(positions, config.components_for_genome_pos,
                         boost::is_any_of(" ,"), boost::token_compress_on);
            for (auto it = positions.begin(); it != positions.end(); ++it) {
                boost::optional<RtSeq> close_kp1mer = FindCloseKP1mer(gp_,
                                                                                 std::stoi(*it), gp_.k_value);
                if (close_kp1mer) {
                    string locality_folder = fs::append_path(pos_loc_folder, *it + "/");
                    make_dir(locality_folder);
                    WriteKmerComponent(gp_, *close_kp1mer, fs::append_path(locality_folder, pos_name + ".dot"), colorer, labeler_);
                } else {
                    WARN(
                        "Failed to find genome kp1mer close to the one at position "
                        << *it << " in the graph. Which is " << RtSeq (gp_.k_value + 1, gp_.genome.GetSequence(), std::stoi(*it)));
                }
            }
        }
    }

    conj_graph_pack& gp_;
    const visualization::graph_labeler::GraphLabeler<Graph>& labeler_;
    string folder_;
};

inline
std::string ConstructComponentName(std::string file_name, size_t cnt) {
    stringstream ss;
    ss << cnt;
    string res = file_name;
    res.insert(res.length(), ss.str());
    return res;
}

template<class Graph>
double AvgCoverage(const Graph& g,
                   const std::vector<typename Graph::EdgeId>& edges) {
    double total_cov = 0.;
    size_t total_length = 0;
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        total_cov += g.coverage(*it) * (double) g.length(*it);
        total_length += g.length(*it);
    }
    return total_cov / (double) total_length;
}

template<class Graph>
size_t Nx(Graph &g, double percent) {
    size_t sum_edge_length = 0;
    vector<size_t> lengths;
    for (auto iterator = g.ConstEdgeBegin(); !iterator.IsEnd(); ++iterator) {
        lengths.push_back(g.length(*iterator));
        sum_edge_length += g.length(*iterator);
    }
    sort(lengths.begin(), lengths.end());
    double len_perc = (1.0 - percent * 0.01) * (double) (sum_edge_length);
    for (size_t i = 0; i < lengths.size(); i++) {
        if (lengths[i] >= len_perc)
            return lengths[i];
        else
            len_perc -= (double) lengths[i];
    }
    return 0;
}

}
}
