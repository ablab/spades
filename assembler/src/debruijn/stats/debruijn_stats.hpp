//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "statistics.hpp"
#include "debruijn_graph.hpp"

#include "graph_pack.hpp"
#include "sequence_mapper.hpp"
#include "graphio.hpp"

#include "omni/visualization/visualization.hpp"
#include "omni/edges_position_handler.hpp"
#include "omni/graph_component.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/delegating_reader_wrapper.hpp"
#include "io/io_helper.hpp"
#include "io/wrapper_collection.hpp"
#include "io/osequencestream.hpp"
#include "dataset_readers.hpp"
#include "copy_file.hpp"

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
    NewExtendedSequenceMapper<Graph, Index> srt(g, index, kmer_mapper);
    return srt.MapSequence(genome);
}

template<class graph_pack>
MappingPath<typename graph_pack::graph_t::EdgeId>
FindGenomeMappingPath(const Sequence& genome, const graph_pack& gp) {
    return FindGenomeMappingPath(genome, gp.g, gp.index, gp.kmer_mapper);
}

template <class graph_pack>
shared_ptr<omnigraph::visualization::GraphColorer<Graph>> DefaultColorer(const graph_pack& gp) {
    return omnigraph::visualization::DefaultColorer(gp.g, 
        FindGenomeMappingPath(gp.genome, gp.g, gp.index, gp.kmer_mapper).path(),
        FindGenomeMappingPath(!gp.genome, gp.g, gp.index, gp.kmer_mapper).path());
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
    GenomeMappingStat(const Graph &graph, const Index &index,	Sequence genome, size_t k) :
            graph_(graph), index_(index), genome_(genome), k_(k) {}

    virtual ~GenomeMappingStat() {}

    virtual void Count() {
        INFO("Mapping genome");
        size_t break_number = 0;
        size_t covered_kp1mers = 0;
        size_t fail = 0;
        if (genome_.size() <= k_)
            return;

        runtime_k::RtSeq cur = genome_.start<runtime_k::RtSeq>(k_ + 1);
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
                   std::shared_ptr<omnigraph::visualization::GraphColorer<Graph>> genome_colorer,
                   const omnigraph::GraphLabeler<Graph>& labeler) {
    INFO("Writing error localities for graph to folder " << folder_name);
    GraphComponent<Graph> all(g, g.begin(), g.end());
    set<EdgeId> edges = genome_colorer->ColoredWith(all.edges().begin(),
                                                    all.edges().end(), "black");
    set<typename Graph::VertexId> to_draw;
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        to_draw.insert(g.EdgeEnd(*it));
        to_draw.insert(g.EdgeStart(*it));
    }
    shared_ptr<GraphSplitter<Graph>> splitter = StandardSplitter(g, to_draw);
    WriteComponents(g, folder_name, splitter, genome_colorer, labeler);
    INFO("Error localities written written to folder " << folder_name);
}

template<class graph_pack>
void CountStats(const graph_pack& gp) {
    typedef typename graph_pack::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;
    INFO("Counting stats");
    StatList stats;
    Path<EdgeId> path1 = FindGenomeMappingPath(gp.genome, gp.g, gp.index,
                                      gp.kmer_mapper).path();
    Path<EdgeId> path2 = FindGenomeMappingPath(!gp.genome, gp.g, gp.index,
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
                                     const GraphLabeler<Graph>& labeler,
                                     const string& folder,
                                     const Path<typename Graph::EdgeId>& path1,
                                     const Path<typename Graph::EdgeId>& path2) {
    INFO("Writing graph components along genome");

    make_dir(folder);
    omnigraph::visualization::WriteComponentsAlongPath(g, path1, folder, omnigraph::visualization::DefaultColorer(g, path1, path2), labeler);

    INFO("Writing graph components along genome finished");
}

//todo refactoring needed: use graph pack instead!!!
template<class Graph, class Mapper>
void WriteGraphComponentsAlongContigs(const Graph& g,
                                      Mapper &mapper,
                                      const std::string& folder,
                                      std::shared_ptr<omnigraph::visualization::GraphColorer<Graph>> colorer,
                                      const GraphLabeler<Graph>& labeler) {
    INFO("Writing graph components along contigs");
    auto contigs_to_thread = io::EasyStream(cfg::get().pos.contigs_to_analyze, false);
    contigs_to_thread->reset();
    io::SingleRead read;
    while (!contigs_to_thread->eof()) {
        (*contigs_to_thread) >> read;
        make_dir(folder + read.name());
        omnigraph::visualization::WriteComponentsAlongPath(g, mapper.MapSequence(read.sequence()).path(), folder + read.name() + "/",
                                                           colorer, labeler);
    }
    INFO("Writing graph components along contigs finished");
}

template<class Graph>
void WriteKmerComponent(conj_graph_pack &gp, runtime_k::RtSeq const& kp1mer, const std::string& file,
                        std::shared_ptr<omnigraph::visualization::GraphColorer<Graph>> colorer,
                        const omnigraph::GraphLabeler<Graph>& labeler) {
    if(!gp.index.contains(kp1mer)) {
        WARN("no such kmer in the graph");
        return;
    }
    VERIFY(gp.index.contains(kp1mer));
    auto pos = gp.index.get(kp1mer);
    VertexId v = pos.second * 2 < gp.g.length(pos.first) ? gp.g.EdgeStart(pos.first) : gp.g.EdgeEnd(pos.first);
    GraphComponent<Graph> component = omnigraph::VertexNeighborhood<Graph>(gp.g, v);
    omnigraph::visualization::WriteComponent<Graph>(component, file, colorer, labeler);
}

inline
optional<runtime_k::RtSeq> FindCloseKP1mer(const conj_graph_pack &gp,
                                           size_t genome_pos, size_t k) {
    VERIFY(gp.genome.size() > 0);
    VERIFY(genome_pos < gp.genome.size());
    static const size_t magic_const = 200;
    for (size_t diff = 0; diff < magic_const; diff++) {
        for (int dir = -1; dir <= 1; dir += 2) {
            size_t pos = (gp.genome.size() - k + genome_pos + dir * diff) % (gp.genome.size() - k);
            runtime_k::RtSeq kp1mer = gp.kmer_mapper.Substitute(
                runtime_k::RtSeq (k + 1, gp.genome, pos));
            if (gp.index.contains(kp1mer))
                return optional<runtime_k::RtSeq>(kp1mer);
        }
    }
    return boost::none;
}

inline
void ProduceDetailedInfo(conj_graph_pack &gp,
                         const omnigraph::GraphLabeler<Graph>& labeler, const string& run_folder,
                         const string &pos_name,
                         info_printer_pos pos) {
    string base_folder = path::append_path(run_folder, "pictures/");
    make_dir(base_folder);
    string folder = path::append_path(base_folder, pos_name + "/");

    auto it = cfg::get().info_printers.find(pos);
    VERIFY(it != cfg::get().info_printers.end());

    const debruijn_config::info_printer & config = it->second;


    if (config.print_stats) {
        INFO("Printing statistics for " << details::info_printer_pos_name(pos));
        CountStats(gp);
    }

    auto path1 = FindGenomeMappingPath(gp.genome, gp.g, gp.index,
                                      gp.kmer_mapper).path();

    auto colorer = omnigraph::visualization::DefaultColorer(gp.g);

    if (config.write_error_loc ||
        config.write_full_graph ||
        config.write_full_nc_graph ||
        config.write_components ||
        !config.components_for_kmer.empty() ||
        config.write_components_along_genome ||
        config.write_components_along_contigs || config.save_full_graph ||
        !config.components_for_genome_pos.empty()) {
        colorer = DefaultColorer(gp);
        make_dir(folder);
    }

    if (config.write_error_loc) {
        make_dir(folder + "error_loc/");
        WriteErrorLoc(gp.g, folder + "error_loc/", colorer, labeler);
    }

    if (config.write_full_graph) {
        WriteComponent(GraphComponent<Graph>(gp.g, gp.g.begin(), gp.g.end()), folder + "full_graph.dot", colorer, labeler);
    }

    if (config.write_full_nc_graph) {
        WriteSimpleComponent(GraphComponent<Graph>(gp.g, gp.g.begin(), gp.g.end()), folder + "nc_full_graph.dot", colorer, labeler);
    }

    if (config.write_components) {
        make_dir(folder + "components/");
        omnigraph::visualization::WriteComponents(gp.g, folder + "components/", omnigraph::ReliableSplitter<Graph>(gp.g), colorer, labeler);
    }

    if (!config.components_for_kmer.empty()) {
        string kmer_folder = path::append_path(base_folder, "kmer_loc/");
        make_dir(kmer_folder);
        auto kmer = runtime_k::RtSeq(gp.k_value + 1, config.components_for_kmer.substr(0, gp.k_value + 1).c_str());
        string file_name = path::append_path(kmer_folder, pos_name + ".dot");
        WriteKmerComponent(gp, kmer, file_name, colorer,labeler);
    }

    if (config.write_components_along_genome) {
        make_dir(folder + "along_genome/");
        omnigraph::visualization::WriteComponentsAlongPath(gp.g, path1, folder + "along_genome/", colorer, labeler);
    }

    if (config.write_components_along_contigs) {
        make_dir(folder + "along_contigs/");
        NewExtendedSequenceMapper<Graph, Index> mapper(gp.g, gp.index, gp.kmer_mapper);
        WriteGraphComponentsAlongContigs(gp.g, mapper, folder + "along_contigs/", colorer, labeler);
    }

    if (config.save_full_graph) {
        make_dir(folder + "full_graph_save/");
        graphio::PrintGraphPack(folder + "full_graph_save/graph", gp);
    }

    if (!config.components_for_genome_pos.empty()) {
        string pos_loc_folder = path::append_path(base_folder, "pos_loc/");
        make_dir(pos_loc_folder);
        vector<string> positions;
        boost::split(positions, config.components_for_genome_pos,
                     boost::is_any_of(" ,"), boost::token_compress_on);
        for (auto it = positions.begin(); it != positions.end(); ++it) {
            optional < runtime_k::RtSeq > close_kp1mer = FindCloseKP1mer(gp,
                                                                         boost::lexical_cast<int>(*it), gp.k_value);
            if (close_kp1mer) {
                string locality_folder = path::append_path(pos_loc_folder, *it + "/");
                make_dir(locality_folder);
                WriteKmerComponent(gp, *close_kp1mer, path::append_path(locality_folder, pos_name + ".dot"), colorer, labeler);
            } else {
                WARN(
                    "Failed to find genome kp1mer close to the one at position "
                    << *it << " in the graph. Which is " << runtime_k::RtSeq (gp.k_value + 1, gp.genome, boost::lexical_cast<int>(*it)));
            }
        }
    }
}

struct detail_info_printer {
    detail_info_printer(conj_graph_pack &gp,
                        const omnigraph::GraphLabeler<Graph>& labeler, const string& folder)
            :  folder_(folder),
               func_(bind(&ProduceDetailedInfo, boost::ref(gp),
                     boost::ref(labeler), _3, _2, _1)),
               graph_(gp.g), cnt(0) {
    }

    void operator()(info_printer_pos pos,
                    string const& folder_suffix = "") {
        cnt++;
        string pos_name = details::info_printer_pos_name(pos);
        VertexEdgeStat<conj_graph_pack::graph_t> stats(graph_);
        TRACE("Number of vertices : " << stats.vertices() << ", number of edges : " << stats.edges() << ", sum length of edges : " << stats.edge_length());
        func_(pos,
              ToString(cnt, 2) + "_" + pos_name + folder_suffix,
              folder_
            //                (path::append_path(folder_, (pos_name + folder_suffix)) + "/")
              );
    }

  private:
    string folder_;
    boost::function<void(info_printer_pos, string const&, string const&)> func_;
    const conj_graph_pack::graph_t &graph_;
    size_t cnt;
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
