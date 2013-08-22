//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "statistics.hpp"
#include "debruijn_graph.hpp"
#include "graph_construction.hpp"
#include "graphio.hpp"
#include "graph_read_correction.hpp"
#include "mismatch_masker.hpp"

#include "omni/visualization/visualization.hpp"
#include "omni/edges_position_handler.hpp"
#include "omni/graph_component.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/delegating_reader_wrapper.hpp"
#include "io/easy_reader.hpp"
#include "io/splitting_wrapper.hpp"
#include "io/wrapper_collection.hpp"
#include "io/osequencestream.hpp"

#include "de/distance_estimation.hpp"
#include "de/pair_info_filters.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include "copy_file.hpp"
#include <cmath>

namespace debruijn_graph {

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
			graph_(graph), index_(index), genome_(genome), k_(k) {
	}

	virtual ~GenomeMappingStat() {
	}

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

template<class Graph, class Index>
class StatCounter: public AbstractStatCounter {
private:
	StatList stats_;
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	StatCounter(const Graph& graph, const Index& index,
	const Sequence& genome, size_t k) {
		SimpleSequenceMapper<Graph, Index> sequence_mapper(graph, index, k + 1);
		Path<EdgeId> path1 = sequence_mapper.MapSequence(Sequence(genome));
		Path<EdgeId> path2 = sequence_mapper.MapSequence(!Sequence(genome));
		stats_.AddStat(new VertexEdgeStat<Graph>(graph));
		stats_.AddStat(new BlackEdgesStat<Graph>(graph, path1, path2));
		stats_.AddStat(new NStat<Graph>(graph, path1, 50));
		stats_.AddStat(new SelfComplementStat<Graph>(graph));
		stats_.AddStat(
				new GenomeMappingStat<Graph, Index>(graph, index,
						Sequence(genome), k));
		stats_.AddStat(new IsolatedEdgesStat<Graph>(graph, path1, path2));
	}

	virtual ~StatCounter() {
		stats_.DeleteStats();
	}

	virtual void Count() {
		stats_.Count();
	}

private:
	DECL_LOGGER("StatCounter")
};

template<class Graph, class Index>
void CountStats(const Graph& g, const Index& index,
const Sequence& genome, size_t k) {
	INFO("Counting stats");
	StatCounter<Graph, Index> stat(g, index, genome, k);
	stat.Count();
	INFO("Stats counted");
}

void CountPairedInfoStats(const Graph& g,
    const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
    const PairedInfoIndexT<Graph>& paired_index,
    const PairedInfoIndexT<Graph>& etalon_index,
    const string& output_folder) {
  PairedInfoIndexT<Graph> filtered_index(g);
	PairInfoWeightFilter<Graph>(g, 40).Filter(paired_index, filtered_index);
	INFO("Counting paired info stats");
	EdgePairStat<Graph>(g, paired_index, output_folder).Count();

	//todo remove filtration if launch on etalon info is ok
	UniquePathStat<Graph>(g, filtered_index,
	                    (size_t)math::round(lib.data().mean_insert_size),
                        lib.data().read_length,
                        0.1 * lib.data().mean_insert_size).Count();
	UniqueDistanceStat<Graph>(etalon_index).Count();
	INFO("Paired info stats counted");
}

// leave only those pairs, which edges have no path in the graph between them
void FilterIndexWithExistingPaths(PairedIndexT& scaf_clustered_index,
                                  const PairedIndexT& index,
                                  const conj_graph_pack &gp,
                                  const GraphDistanceFinder<Graph>& dist_finder) {
  for (auto it = index.begin(); it != index.end(); ++it) {
    const std::set<Point>& histogram = *it;
    EdgeId e1 = it.first();
    EdgeId e2 = it.second();
    if (gp.g.OutgoingEdgeCount(gp.g.EdgeEnd(e1)) == 0 && gp.g.IncomingEdgeCount(gp.g.EdgeEnd(e1)) == 1 &&
        gp.g.IncomingEdgeCount(gp.g.EdgeStart(e2)) == 0 && gp.g.OutgoingEdgeCount(gp.g.EdgeStart(e2)) == 1)     {
      vector<size_t> dists = dist_finder.GetGraphDistancesLengths(e1, e2);
      if (dists.size() == 0)
        for (auto point_iter = histogram.begin(); point_iter != histogram.end(); ++point_iter)
          if (math::gr(point_iter->d, 0.)) {
            scaf_clustered_index.AddPairInfo(it.first(), it.second(),
                                             point_iter->d, point_iter->weight, 20.);
          }
    }
  }
}
/*
void FillAndCorrectEtalonPairedInfo(
    PairedIndexT&  corrected_etalon_index, const conj_graph_pack& gp,
    const PairedIndexT&  paired_index, size_t insert_size,
		size_t read_length, size_t delta,
		bool save_etalon_info_history = false) {
	INFO("Filling etalon paired index");
  PairedIndexT etalon_index(gp.g);
    bool successful_load = false;
    if (cfg::get().entry_point >= ws_distance_estimation) {
        string p = path::append_path(cfg::get().load_from, "../etalon");
        if (!path::is_regular_file(p + ".prd")) {
            DEBUG("file " << p + ".prd" << " does not exist");
        }
        else {
            INFO("Loading etalon pair info from the previous run...");
            Graph& graph = const_cast<Graph&>(gp.g);
            IdTrackHandler<Graph>& int_ids = const_cast<IdTrackHandler<Graph>& >(gp.int_ids);
            ScannerTraits<Graph>::Scanner scanner(graph, int_ids);
            scanner.loadPaired(p, etalon_index);
            path::files_t files;
            files.push_back(p);
            copy_files_by_prefix(files, cfg::get().output_dir);
            successful_load = true;
        }
    }
    if (!successful_load)
	    FillEtalonPairedIndex(etalon_index, gp.g,
			gp.index, gp.kmer_mapper, insert_size, read_length, delta,
			gp.genome, gp.k_value);
	INFO("Etalon paired index filled");

	INFO("Correction of etalon paired info has been started");

	INFO("Filtering etalon info");
	//leave only info between edges both present in paired_index
  PairedIndexT filtered_etalon_index(gp.g);
  for (auto iter = etalon_index.begin(); iter != etalon_index.end(); ++iter) {
    const set<Point>& histogram = *iter;
    EdgeId first_edge = iter.first();
    EdgeId second_edge = iter.second();
    if (paired_index.GetEdgePairInfo(first_edge, second_edge).size() > 0) {
      for (auto point = histogram.begin(); point != histogram.end(); ++point)
        filtered_etalon_index.AddPairInfo(first_edge, second_edge, *point);
    }
    else
      DEBUG("Filtering out pair_info " << gp.g.int_id(first_edge) << " "
                                       << gp.g.int_id(second_edge));
	}

	INFO("Pushing etalon info through estimator");
	GraphDistanceFinder<Graph> dist_finder(gp.g, insert_size, read_length, delta);
	DistanceEstimator<Graph> estimator(gp.g, filtered_etalon_index, dist_finder, 0., 4.);
	estimator.Estimate(corrected_etalon_index);
	if (save_etalon_info_history) {
		INFO("Saving etalon paired info indices on different stages");
		ConjugateDataPrinter<Graph> data_printer(gp.g, gp.int_ids);
    data_printer.savePaired(cfg::get().output_dir + "etalon", etalon_index);
		data_printer.savePaired(
				cfg::get().output_dir + "etalon_filtered_by_index",
				filtered_etalon_index);
    data_printer.savePaired(cfg::get().output_dir + "etalon_corrected_by_graph",
				corrected_etalon_index);
		INFO("Everything is saved");

        if (cfg::get().paired_info_scaffolder) {
	        GraphDistanceFinder<Graph> dist_finder(gp.g, insert_size, read_length, delta);
        	INFO("Saving paired information statistics for a scaffolding");
            PairedIndexT scaf_etalon_index(gp.g);
            FilterIndexWithExistingPaths(scaf_etalon_index, etalon_index, gp, dist_finder);
            data_printer.savePaired(
                    cfg::get().output_dir + "scaf_etalon",
                    scaf_etalon_index);
            PairedIndexT scaf_filtered_etalon_index(gp.g);
            FilterIndexWithExistingPaths(scaf_filtered_etalon_index, filtered_etalon_index, gp, dist_finder);
			data_printer.savePaired(
					cfg::get().output_dir + "scaf_etalon_filtered",
					scaf_filtered_etalon_index);
        }

        INFO("Everything saved");
	}
	INFO("Correction finished");
}*/

template<class Graph>
void GetAllDistances(const PairedInfoIndexT<Graph>& paired_index,
                     PairedInfoIndexT<Graph>& result,
                     const GraphDistanceFinder<Graph>& dist_finder) {
    for (auto iter = paired_index.begin(); iter != paired_index.end(); ++iter) {
        EdgeId e1 = iter.first();
        EdgeId e2 = iter.second();
        vector<size_t> forward = dist_finder.GetGraphDistancesLengths(e1, e2);
        for (size_t i = 0; i < forward.size(); ++i)
            result.AddPairInfo(e1, e2, (double) forward[i], -10.0, 0.0, false);
    }
}

template<class Graph>
void GetAllDistances(const Graph& g,
                     const PairedInfoIndexT<Graph>& paired_index,
                     const PairedInfoIndexT<Graph>& clustered_index,
                     const IdTrackHandler<Graph>& int_ids,
                     const GraphDistanceFinder<Graph>& dist_finder,
                     PairedInfoIndexT<Graph>& result)
{
  typedef typename Graph::EdgeId EdgeId;
  typedef set<Point> Histogram;
  typedef vector<EdgeId> Path;
	for (auto iter = paired_index.begin(); iter != paired_index.end(); ++iter) {
    EdgeId first = iter.first();
    EdgeId second = iter.second();
    const vector<Path>& raw_paths = dist_finder.GetGraphDistances(first, second);
        // adding first edge to every path
    vector<Path> paths;
        for (size_t i = 0; i < raw_paths.size(); ++i) {
      Path path;
            path.push_back(first);
      for (size_t j = 0; j < raw_paths[i].size(); ++j)
                path.push_back(raw_paths[i][j]);
            path.push_back(second);

            paths.push_back(path);
        }
        vector<size_t> path_lengths;
        vector<double> path_weights;
        for (size_t i = 0; i < paths.size(); ++i) {
            size_t len_total = 0 ;
            double weight_total = 0.;
            for (size_t j = 0; j < paths[i].size(); ++j) {
                len_total += g.length(paths[i][j]);
                size_t cur_length = 0;
                for (size_t l = j + 1; l < paths[i].size(); ++l) {
                    cur_length += g.length(paths[i][l - 1]);
          const Histogram& infos = clustered_index.GetEdgePairInfo(paths[i][j], paths[i][l]);
                    for (auto iterator = infos.begin(); iterator != infos.end(); ++iterator) {
            const Point& info = *iterator;
                        if (info.d == cur_length) {
                            weight_total += info.weight;
                            break;
                        }
                    }
                }
            }
            path_lengths.push_back(len_total - g.length(second));
            path_weights.push_back(weight_total);
        }

        for (size_t i = 0; i < paths.size(); ++i) {
            cout << int_ids.ReturnIntId(first) << "(" << g.length(first) << ") "
        << int_ids.ReturnIntId(second) << "(" << g.length(second) << ") : "
        << (i + 1) << "-th path (" << path_lengths[i] << ", " << path_weights[i] << ")   :::   ";
            for (size_t j = 0; j < paths[i].size(); ++j) {
                cout << int_ids.ReturnIntId(paths[i][j]) << "(" << g.length(paths[i][j]) << ") ";
            }
            cout << endl;
		}
    }
}

template<class Graph>
void CountAndSaveAllPaths(const Graph& g, const io::SequencingLibrary<debruijn_config::DataSetData> &lib, const IdTrackHandler<Graph>& int_ids,
    const PairedInfoIndexT<Graph>& paired_index, const PairedInfoIndexT<Graph>& /*clustered_index*/) {
  PairedIndexT all_paths(g);
  GetAllDistances<Graph>(paired_index,
                         all_paths,
                         GraphDistanceFinder<Graph>(g,
                                                    size_t(lib.data().mean_insert_size),
                                                    lib.data().read_length,
                                                    size_t(lib.data().insert_size_deviation)));

  std::string dir_name = cfg::get().output_dir + "estimation_qual/";
	make_dir(dir_name);

	typename PrinterTraits<Graph>::Printer printer(g, int_ids);
  printer.savePaired(dir_name + "paths", all_paths);

  //PairedIndexT& all_paths_2(g);
  //GetAllDistances<Graph>(g,
            //paired_index, clustered_index,
            //int_ids,
            //all_paths_2,
            //GraphDistanceFinder<Graph>(g, *cfg::get().ds.IS, *cfg::get().ds.RL,
            //size_t(*cfg::get().ds.is_var)));
	//printer.savePaired(dir_name + "paths_all", all_paths_2);
}

/*
void CountClusteredPairedInfoStats(const conj_graph_pack &gp,
    const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
    const PairedInfoIndexT<Graph> &paired_index,
    const PairedInfoIndexT<Graph> &clustered_index) {
  PairedIndexT etalon_index(gp.g);

  FillAndCorrectEtalonPairedInfo(etalon_index, gp, paired_index,
                                 (size_t)math::round(lib.data().mean_insert_size),
                                 lib.data().read_length,
                                 (size_t)math::round(lib.data().insert_size_deviation), true);

	CountAndSaveAllPaths(gp.g, lib, gp.int_ids, paired_index, clustered_index);

	INFO("Counting clustered info stats");
	EdgeQuality<Graph, Index> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
  //EstimationQualityStat<Graph> estimation_stat(gp.g, gp.int_ids, edge_qual,
                                              //paired_index, clustered_index, etalon_index);
  //estimation_stat.Count();
  //estimation_stat.SaveStats(cfg::get().output_dir + "estimation_qual/");

	INFO("Counting overall cluster stat");
	ClusterStat<Graph>(clustered_index).Count();
	INFO("Overall cluster stat");

  if (cfg::get().paired_info_scaffolder) {
		ConjugateDataPrinter<Graph> data_printer(gp.g, gp.int_ids);
    INFO("Generating the statistics of pair info for scaffolding");
    PairedIndexT scaf_clustered_index(gp.g);
    FilterIndexWithExistingPaths(scaf_clustered_index,
                                 clustered_index, gp,
                                 GraphDistanceFinder<Graph>(gp.g,
                                         (size_t)math::round(lib.data().mean_insert_size),
                                         lib.data().read_length,
                                         (size_t)math::round(lib.data().insert_size_deviation)));
    data_printer.savePaired(cfg::get().output_dir + "scaf_clustered",
                            scaf_clustered_index);
  }
  //  PairedInfoIndexT<Graph> etalon_clustered_index;
	//	DistanceEstimator<Graph> estimator(g, etalon_index, insert_size,
	//			max_read_length, cfg::get().de.delta,
	//			cfg::get().de.linkage_distance, cfg::get().de.max_distance);
	//	estimator.Estimate(etalon_clustered_index);

  //  PairedInfoIndexT<Graph> filtered_clustered_index(g);
	//	PairInfoFilter<Graph> (g, 1000.).Filter(
  //      clustered_index[>etalon_clustered_index<], filtered_clustered_index);
	INFO("Counting mate-pair transformation stat");
	MatePairTransformStat<Graph>(gp.g, //filtered_
	    clustered_index).Count();
	INFO("Mate-pair transformation stat counted");
	INFO("Clustered info stats counted");
}
*/

void WriteErrorLoc(const Graph &g,
		const string& folder_name,
		shared_ptr<omnigraph::visualization::GraphColorer<Graph>> genome_colorer,
		const omnigraph::GraphLabeler<Graph>& labeler) {
	INFO("Writing error localities for graph to folder " << folder_name);
	GraphComponent<Graph> all(g, g.begin(), g.end());
	set<EdgeId> edges = genome_colorer->ColoredWith(all.edges().begin(),
			all.edges().end(), "black");
	set<Graph::VertexId> to_draw;
	for (auto it = edges.begin(); it != edges.end(); ++it) {
		to_draw.insert(g.EdgeEnd(*it));
		to_draw.insert(g.EdgeStart(*it));
	}
	shared_ptr<GraphSplitter<Graph>> splitter = StandardSplitter(g, to_draw);
	WriteComponents(g, folder_name, splitter, genome_colorer, labeler);
	INFO("Error localities written written to folder " << folder_name);
}

template<class Graph, class Index>
Path<typename Graph::EdgeId> FindGenomePath(const Sequence& genome,
		const Graph& g, const Index& index, size_t k) {
	SimpleSequenceMapper<Graph, Index> srt(g, index, k + 1);
	return srt.MapSequence(genome);
}

template<class Graph, class Index>
MappingPath<typename Graph::EdgeId> FindGenomeMappingPath(
		const Sequence& genome, const Graph& g,
		const Index& index,
		const KmerMapper<Graph>& kmer_mapper) {
	NewExtendedSequenceMapper<Graph, Index> srt(g, index, kmer_mapper, g.k() + 1);
	return srt.MapSequence(genome);
}

void WriteGraphComponentsAlongGenome(const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& folder,
		const Path<Graph::EdgeId>& path1,
		const Path<Graph::EdgeId>& path2) {

	INFO("Writing graph components along genome");

	make_dir(folder);
	omnigraph::visualization::WriteComponentsAlongPath(g, path1, folder, omnigraph::visualization::DefaultColorer(g, path1, path2), labeler);

	INFO("Writing graph components along genome finished");
}

//todo refactoring needed: use graph pack instead!!!
template<class Graph, class Mapper>
void WriteGraphComponentsAlongContigs(const Graph& g,
		Mapper &mapper,
		const string& folder,
		shared_ptr<omnigraph::visualization::GraphColorer<Graph>> colorer,
		const GraphLabeler<Graph>& labeler) {
	INFO("Writing graph components along contigs");
	io::EasyReader contigs_to_thread(cfg::get().pos.contigs_to_analyze, false/*true*/);
	contigs_to_thread.reset();
	io::SingleRead read;
	while (!contigs_to_thread.eof()) {
		contigs_to_thread >> read;
		make_dir(folder + read.name());
		omnigraph::visualization::WriteComponentsAlongPath(g, mapper.MapSequence(read.sequence()).simple_path(), folder + read.name() + "/",
				colorer, labeler);
	}
	INFO("Writing graph components along contigs finished");
}

void WriteKmerComponent(conj_graph_pack &gp, runtime_k::RtSeq const& kp1mer, const string& file,
		shared_ptr<omnigraph::visualization::GraphColorer<Graph>> colorer,
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

optional<runtime_k::RtSeq> FindCloseKP1mer(const conj_graph_pack &gp,
		size_t genome_pos, size_t k) {
	static const size_t magic_const = 200;
	for (size_t diff = 0; diff < magic_const; diff++) {
		for (int dir = -1; dir <= 1; dir += 2) {
			size_t pos = genome_pos + dir * diff;
			runtime_k::RtSeq kp1mer = gp.kmer_mapper.Substitute(
			        runtime_k::RtSeq (k + 1, gp.genome, pos));
			if (gp.index.contains(kp1mer))
				return optional<runtime_k::RtSeq>(kp1mer);
		}
	}
	return none;
}

void ProduceDetailedInfo(conj_graph_pack &gp,
		const omnigraph::GraphLabeler<Graph>& labeler, const string& run_folder,
		const string &pos_name,
		info_printer_pos pos,
		size_t k) {
	string base_folder = path::append_path(run_folder, "pictures/");
	make_dir(base_folder);
	string folder = path::append_path(base_folder, pos_name + "/");

	auto it = cfg::get().info_printers.find(pos);
	VERIFY(it != cfg::get().info_printers.end());

	const debruijn_config::info_printer & config = it->second;


	if (config.print_stats) {
		INFO("Printing statistics for " << details::info_printer_pos_name(pos));
		CountStats(gp.g, gp.index, gp.genome, k);
	}

	typedef Path<Graph::EdgeId> path_t;
	path_t path1;
	path_t path2;
	shared_ptr<omnigraph::visualization::GraphColorer<Graph>> colorer = omnigraph::visualization::DefaultColorer(gp.g);

	if (config.write_error_loc
			|| config.write_full_graph
			|| config.write_full_nc_graph
			|| config.write_components
			|| !config.components_for_kmer.empty()
			|| config.write_components_along_genome
			|| config.write_components_along_contigs || config.save_full_graph
			|| !config.components_for_genome_pos.empty()) {
		path1 = FindGenomeMappingPath(gp.genome, gp.g, gp.index,
				gp.kmer_mapper).simple_path();
		path2 = FindGenomeMappingPath(!gp.genome, gp.g, gp.index,
				gp.kmer_mapper).simple_path();
		colorer = omnigraph::visualization::DefaultColorer(gp.g, path1, path2);
//		path1 = FindGenomePath<K>(gp.genome, gp.g, gp.index);
//		path2 = FindGenomePath<K>(!gp.genome, gp.g, gp.index);
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
		auto kmer = runtime_k::RtSeq(k + 1, config.components_for_kmer.substr(0, k + 1).c_str());
		string file_name = path::append_path(kmer_folder, pos_name + ".dot");
		WriteKmerComponent(gp, kmer, file_name, colorer,labeler);
	}

	if (config.write_components_along_genome) {
		make_dir(folder + "along_genome/");
		omnigraph::visualization::WriteComponentsAlongPath(gp.g, path1, folder + "along_genome/", colorer, labeler);
	}

	if (config.write_components_along_contigs) {
		make_dir(folder + "along_contigs/");
		NewExtendedSequenceMapper<Graph, Index> mapper(gp.g, gp.index, gp.kmer_mapper, gp.g.k() + 1);
		WriteGraphComponentsAlongContigs(gp.g, mapper, folder + "along_contigs/", colorer, labeler);
	}

	if (config.save_full_graph) {
		make_dir(folder + "full_graph_save/");
		ConjugateDataPrinter<Graph> printer(gp.g, gp.int_ids);
		PrintGraphPack(folder + "full_graph_save/graph", printer, gp);
	}

	if (!config.components_for_genome_pos.empty()) {
		string pos_loc_folder = path::append_path(base_folder, "pos_loc/");
		make_dir(pos_loc_folder);
		vector<string> positions;
		boost::split(positions, config.components_for_genome_pos,
				boost::is_any_of(" ,"), boost::token_compress_on);
		for (auto it = positions.begin(); it != positions.end(); ++it) {
			optional < runtime_k::RtSeq > close_kp1mer = FindCloseKP1mer(gp,
					boost::lexical_cast<int>(*it), k);
			if (close_kp1mer) {
				string locality_folder = path::append_path(pos_loc_folder, *it + "/");
				make_dir(locality_folder);
				WriteKmerComponent(gp, *close_kp1mer, path::append_path(locality_folder, pos_name + ".dot"), colorer, labeler);
			} else {
				WARN(
						"Failed to find genome kp1mer close to the one at position "
								<< *it << " in the graph. Which is " << runtime_k::RtSeq (k + 1, gp.genome, boost::lexical_cast<int>(*it)));
			}
		}
	}
}

struct detail_info_printer {
	detail_info_printer(conj_graph_pack &gp,
			const omnigraph::GraphLabeler<Graph>& labeler, const string& folder)

	:
			folder_(folder), func_(
					bind(&ProduceDetailedInfo, boost::ref(gp),
							boost::ref(labeler), _3, _2, _1, gp.k_value)), graph_(
					gp.g), cnt(0) {
	}

	void operator()(info_printer_pos pos,
			string const& folder_suffix = "") {
		cnt++;
		string pos_name = details::info_printer_pos_name(pos);
		VertexEdgeStat<conj_graph_pack::graph_t> stats(graph_);
		TRACE("Number of vertices : " << stats.vertices() << ", number of edges : " << stats.edges() << ", sum length of edges : " << stats.edge_length());
		func_(
				pos,
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

string ConstructComponentName(string file_name, size_t cnt) {
	stringstream ss;
	ss << cnt;
	string res = file_name;
	res.insert(res.length(), ss.str());
	return res;
}

template<class graph_pack>
int PrintGraphComponents(const string& file_name, graph_pack& gp,
    size_t split_edge_length, PairedInfoIndexT<Graph> &clustered_index) {
    shared_ptr<GraphSplitter<Graph>> inner_splitter = ReliableSplitter<Graph>(gp.g, split_edge_length);
    shared_ptr<GraphComponentFilter<Graph>> checker = make_shared<ComponentSizeFilter<Graph>>(gp.g, split_edge_length, 2, 300);
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, checker);
	size_t cnt = 1;
	while (splitter.HasNext() && cnt <= 1000) {
		string component_name = ConstructComponentName(file_name, cnt).c_str();
		auto component = splitter.Next();
		PrintWithClusteredIndex(component_name, gp, component.vertices().begin(),
				component.vertices().end(), clustered_index);
		cnt++;
	}
	return (int) cnt - 1;
}

template<class Graph>
double AvgCoverage(const Graph& g,
		const vector<typename Graph::EdgeId>& edges) {
	double total_cov = 0.;
	size_t total_length = 0;
	for (auto it = edges.begin(); it != edges.end(); ++it) {
		total_cov += g.coverage(*it) * (double) g.length(*it);
		total_length += g.length(*it);
	}
	return total_cov / (double) total_length;
}

/*
void tSeparatedStats(conj_graph_pack& gp, const Sequence& contig,
		PairedInfoIndex<conj_graph_pack::graph_t> &ind, size_t k) {
	typedef omnigraph::PairInfo<EdgeId> PairInfo;

	MappingPath<Graph::EdgeId> m_path1 = FindGenomeMappingPath(contig, gp.g,
			gp.index, gp.kmer_mapper);

	map<Graph::EdgeId, vector<pair<int, int>>> inGenomeWay;
	int CurI = 0;
	int gaps = 0;
	for (size_t i = 0; i < m_path1.size(); i++) {
		bool new_edge_added = false;
		EdgeId ei = m_path1[i].first;
		MappingRange mr = m_path1[i].second;
		int start = (int)(mr.initial_range.start_pos - mr.mapped_range.start_pos);
		if (inGenomeWay.find(ei) == inGenomeWay.end()) {
			vector<pair<int, int>> tmp;
			tmp.push_back(make_pair(CurI, start));
			inGenomeWay[ei] = tmp;
			CurI++;
			new_edge_added = true;
      DEBUG("Edge " << gp.int_ids.str(ei) << " num " << CurI << " pos " << start);
		} else {
			if (m_path1[i - 1].first == ei) {
        if (abs(start - inGenomeWay[ei][(inGenomeWay[ei].size() - 1)].second) > 50) {
					inGenomeWay[ei].push_back(make_pair(CurI, start));
					CurI++;
					new_edge_added = true;
          DEBUG("Edge " << gp.int_ids.str(ei) << " num " << CurI << " pos " << start);
				}
			} else {
				inGenomeWay[ei].push_back(make_pair(CurI, start));
				CurI++;
				new_edge_added = true;
        DEBUG("Edge " << gp.int_ids.str(ei) << " num " << CurI << " pos " << start);
			}
		}
		if (new_edge_added && (i > 0)) {
			if (gp.g.EdgeStart(ei) != gp.g.EdgeEnd(m_path1[i - 1].first)) {
				gaps++;
			}
		}
  }
  INFO("Totaly " << CurI << " edges in genome path, with " << gaps << "not adjacent conequences");

	vector<int> stats(10);
	vector<int> stats_d(10);
	int PosInfo = 0;
	int AllignedPI = 0;
	int ExactDPI = 0;
	int OurD = (int) cfg::get().ds.IS() - (int) cfg::get().ds.RL();
	for (auto p_iter = ind.begin(), p_end_iter = ind.end();
			p_iter != p_end_iter; ++p_iter) {
		vector<PairInfo> pi = *p_iter;
		for (size_t j = 0; j < pi.size(); j++) {
			EdgeId left_edge = pi[j].first;
			EdgeId right_edge = pi[j].second;
      double d = pi[j].d();
      if (d < 0.001)
				continue;
			int best_d = 100;
			int best_t = 0;
			PosInfo++;
			DEBUG(
          "PairInfo " << gp.int_ids.str(left_edge) << " -- " << gp.int_ids.str(right_edge) << " d " << d);
			bool ExactOnD = false;
			for (size_t left_i = 0; left_i < inGenomeWay[left_edge].size();
					left_i++)
				for (size_t right_i = 0;
						right_i < inGenomeWay[right_edge].size(); right_i++) {
					if (best_d
							> abs(
									inGenomeWay[right_edge][right_i].second
											- inGenomeWay[left_edge][left_i].second
                      - d)) {
						best_d = (int)math::round(abs(
								inGenomeWay[right_edge][right_i].second
										- inGenomeWay[left_edge][left_i].second
                    - d));
						best_t = inGenomeWay[right_edge][right_i].first
								- inGenomeWay[left_edge][left_i].first;
						DEBUG("best d " << best_d);
						if ((inGenomeWay[right_edge][right_i].second
								- inGenomeWay[left_edge][left_i].second
								- (int) gp.g.length(left_edge) <= OurD)
								&& (inGenomeWay[right_edge][right_i].second
										- inGenomeWay[left_edge][left_i].second
										+ (int) gp.g.length(right_edge) >= OurD))
							ExactOnD = true;
						else
							ExactOnD = false;
					}
				}
			if (best_t > 5)
				best_t = 5;
			if (best_d < 100) {
				AllignedPI++;
				stats[best_t]++;
				if (ExactOnD) {
					stats_d[best_t]++;
					ExactDPI++;
				}
			}

		}
	}INFO(
			"Total positive pair info " << PosInfo << " alligned to genome " << AllignedPI << " with exact distance " << ExactDPI);
	INFO(
			"t-separated stats Alligneg: 1 - " << stats[1] << " 2 - " << stats[2] << " 3 - " << stats[3] << " 4 - " << stats[4] << " >4 - " << stats[5]);
	INFO(
			"t-separated stats Exact: 1 - " << stats_d[1] << " 2 - " << stats_d[2] << " 3 - " << stats_d[3] << " 4 - " << stats_d[4] << " >4 - " << stats[5]);
}*/

template<class Graph, class Mapper>
class PosFiller {
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	const Mapper& mapper_;
	EdgesPositionHandler<Graph>& edge_pos_;

public:
	PosFiller(const Graph& g, const Mapper& mapper,
			EdgesPositionHandler<Graph>& edge_pos) :
			g_(g), mapper_(mapper), edge_pos_(edge_pos) {

	}

	void Process(const Sequence& s, string name) const {
		//todo stupid conversion!

		return Process(io::SingleRead(name, s.str()));
	}

	void Process(const io::SingleRead& read) const {
		MappingPath<EdgeId> path = mapper_.MapRead(read);
		const string& name = read.name();
		int cur_pos = 0;
		TRACE(
				"Contig " << name << " mapped on " << path.size()
						<< " fragments.");
		for (size_t i = 0; i < path.size(); i++) {
			EdgeId ei = path[i].first;
			MappingRange mr = path[i].second;
			int len = (int) (mr.mapped_range.end_pos - mr.mapped_range.start_pos);
			if (i > 0)
				if (path[i - 1].first != ei)
					if (g_.EdgeStart(ei) != g_.EdgeEnd(path[i - 1].first)) {
						TRACE(
								"Contig " << name
										<< " mapped on not adjacent edge. Position in contig is "
										<< path[i - 1].second.initial_range.start_pos
												+ 1
										<< "--"
										<< path[i - 1].second.initial_range.end_pos
										<< " and "
										<< mr.initial_range.start_pos + 1
										<< "--" << mr.initial_range.end_pos);
					}
            edge_pos_.AddEdgePosition(ei, (int) mr.initial_range.start_pos + 1,
                                      (int) mr.initial_range.end_pos, name,
                                      (int) mr.mapped_range.start_pos + 1,
                                      (int) mr.mapped_range.end_pos);
			cur_pos += len;
		}
	}

private:
	DECL_LOGGER("PosFiller")
	;
};

template<class Graph, class Mapper>
void FillPos(const Graph& g, const Mapper& mapper,
		EdgesPositionHandler<Graph>& edge_pos,
		io::IReader<io::SingleRead>& stream) {
	PosFiller<Graph, Mapper> filler(g, mapper, edge_pos);
	io::SingleRead read;
	while (!stream.eof()) {
		stream >> read;
		filler.Process(read);
	}
}

template<class gp_t>
void FillPos(gp_t& gp, io::IReader<io::SingleRead>& stream) {
	FillPos(gp.g, *MapperInstance(gp), gp.edge_pos, stream);
}

template<class Graph, class Mapper>
void FillPos(const Graph& g, const Mapper& mapper, EdgesPositionHandler<Graph>& edge_pos, const Sequence& s, const string& name) {
	PosFiller<Graph, Mapper>(g, mapper, edge_pos).Process(s, name);
}

template<class gp_t>
void FillPos(gp_t& gp, const Sequence& s, const string& name) {
	FillPos(gp.g, *MapperInstance(gp), gp.edge_pos, s, name);
}

//deprecated, todo remove usages!!!
template<class gp_t>
void FillPos(gp_t& gp, const string& contig_file, string prefix) {
//	typedef typename gp_t::Graph::EdgeId EdgeId;
	INFO("Threading large contigs");
	io::Reader irs(contig_file);
	while(!irs.eof()) {
		io::SingleRead read;
		irs >> read;
		DEBUG("Contig " << read.name() << ", length: " << read.size());
		if (!read.IsValid()) {
			WARN("Attention: contig " << read.name() << " contains Ns");
			continue;
		}
		Sequence contig = read.sequence();
		if (contig.size() < 1500000) {
			//		continue;
		}
		FillPos(gp, contig, prefix + read.name());
	}
}

template<class gp_t>
void FillPosWithRC(gp_t& gp, const string& contig_file, string prefix) {
//  typedef typename gp_t::Graph::EdgeId EdgeId;
	INFO("Threading large contigs");
	io::EasySplittingReader irs(contig_file, true);
	while(!irs.eof()) {
		io::SingleRead read;
		irs >> read;
		DEBUG("Contig " << read.name() << ", length: " << read.size());
		if (!read.IsValid()) {
			WARN("Attention: contig " << read.name()
					<< " is not valid (possibly contains N's)");
			continue;
		}
		Sequence contig = read.sequence();
		FillPos(gp, contig, prefix + read.name());
	}
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

