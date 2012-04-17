//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "omni/visualization_utils.hpp"
#include "statistics.hpp"
#include "new_debruijn.hpp"
#include "graphio.hpp"
#include "io/easy_reader.hpp"
#include "omni/edges_position_handler.hpp"
#include "omni/distance_estimation.hpp"
#include "omni/graph_component.hpp"
#include "io/delegating_reader_wrapper.hpp"
#include "omni/pair_info_filters.hpp"
#include "io/easy_reader.hpp"
#include "read/osequencestream.hpp"
#include "io/easy_reader.hpp"
#include "k.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

namespace debruijn_graph {

template<class Graph, size_t k>
class GenomeMappingStat: public AbstractStatCounter {
private:
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	const EdgeIndex<k + 1, Graph>& index_;
	Sequence genome_;
public:
	GenomeMappingStat(const Graph &graph, const EdgeIndex<k + 1, Graph> &index,

	Sequence genome) :
			graph_(graph), index_(index), genome_(genome) {
	}

	virtual ~GenomeMappingStat() {
	}

	virtual void Count() {
		INFO("Mapping genome");
		size_t break_number = 0;
		size_t covered_kp1mers = 0;
		size_t fail = 0;
		if (genome_.size() <= k)
			return;
		Seq < k + 1 > cur = genome_.start<k + 1>() >> 0;
		bool breaked = true;
		pair<EdgeId, size_t> cur_position;
		for (size_t cur_nucl = k; cur_nucl < genome_.size(); cur_nucl++) {
			cur = cur << genome_[cur_nucl];
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
		INFO(
				"Covered k+1-mers:" << covered_kp1mers << " of "
						<< (genome_.size() - k) << " which is "
						<< (100.0 * covered_kp1mers / (genome_.size() - k))
						<< "%");
		INFO(
				"Covered k+1-mers form " << break_number + 1
						<< " contigious parts");
		INFO("Continuity failtures " << fail);
	}
};

template<class Graph, size_t k>
class StatCounter: public AbstractStatCounter {
private:
	StatList stats_;
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	StatCounter(const Graph& graph, const EdgeIndex<k + 1, Graph>& index,
	const Sequence& genome) {
		SimpleSequenceMapper<k + 1, Graph> sequence_mapper(graph, index);
		Path<EdgeId> path1 = sequence_mapper.MapSequence(Sequence(genome));
		Path<EdgeId> path2 = sequence_mapper.MapSequence(!Sequence(genome));
		stats_.AddStat(new VertexEdgeStat<Graph>(graph));
		stats_.AddStat(new BlackEdgesStat<Graph>(graph, path1, path2));
		stats_.AddStat(new NStat<Graph>(graph, path1, 50));
		stats_.AddStat(new SelfComplementStat<Graph>(graph));
		stats_.AddStat(
				new GenomeMappingStat<Graph, k>(graph, index,
						Sequence(genome)));
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

template<size_t k, class Graph>
void CountStats(const Graph& g, const EdgeIndex<k + 1, Graph>& index,
const Sequence& genome) {
	INFO("Counting stats");
	StatCounter<Graph, k> stat(g, index, genome);
	stat.Count();
	INFO("Stats counted");
}

void CountPairedInfoStats(const Graph &g,
		const PairedInfoIndex<Graph> &paired_index,
		const PairedInfoIndex<Graph> &etalon_paired_index,
		const string &output_folder) {
	PairedInfoIndex<Graph> filtered_index(g);
	PairInfoWeightFilter<Graph>(g, 40).Filter(paired_index, filtered_index);
	INFO("Counting paired info stats");
	EdgePairStat<Graph>(g, paired_index, output_folder).Count();

	//todo remove filtration if launch on etalon info is ok
	UniquePathStat<Graph>(g, filtered_index, *cfg::get().ds.IS,	*cfg::get().ds.RL, 0.1 * (*cfg::get().ds.IS)).Count();
	UniqueDistanceStat<Graph>(etalon_paired_index).Count();
	INFO("Paired info stats counted");
}

void FillAndCorrectEtalonPairedInfo(
		paired_info_index &final_etalon_paired_index, const conj_graph_pack &gp,
		const paired_info_index &paired_index, size_t insert_size,
		size_t read_length, size_t delta,
		bool save_etalon_info_history = false) {
	INFO("Filling etalon paired index");
	paired_info_index etalon_paired_index(gp.g);
	FillEtalonPairedIndex<debruijn_graph::K>(etalon_paired_index, gp.g,
			gp.index, gp.kmer_mapper, insert_size, read_length, delta,
			gp.genome);
	INFO("Etalon paired index filled");

	INFO("Correction of etalon paired info has been started");

	INFO("Collecting data to filter etalon info");
	std::set<std::pair<Graph::EdgeId, Graph::EdgeId>> setEdgePairs;
	for (auto iter = paired_index.begin(); iter != paired_index.end(); ++iter)
		setEdgePairs.insert(
				std::make_pair((*iter)[0].first, (*iter)[0].second));

	INFO("Filtering etalon info");
	//leave only info between edges both present in paired_index
	paired_info_index filtered_etalon_index(gp.g);
	for (auto iter = etalon_paired_index.begin();
			iter != etalon_paired_index.end(); ++iter) {
		std::vector<omnigraph::PairInfo<EdgeId> > pair_info = *iter;
		if (setEdgePairs.count(
				std::make_pair(pair_info[0].first, pair_info[0].second)) > 0)
			for (auto point = pair_info.begin(); point != pair_info.end();
					point++)
				filtered_etalon_index.AddPairInfo(*point);
	}

	INFO("Pushing etalon info through estimator");
	GraphDistanceFinder<Graph> dist_finder(gp.g, insert_size, read_length,
			delta);
	DistanceEstimator<Graph> estimator(gp.g, filtered_etalon_index, dist_finder,
			0, 4);
	estimator.Estimate(final_etalon_paired_index);
	if (save_etalon_info_history) {
		INFO("Saving etalon paired info indices on different stages");
		ConjugateDataPrinter<Graph> data_printer(gp.g, gp.int_ids);
		data_printer.savePaired(cfg::get().output_dir + "etalon_paired",
				etalon_paired_index);
		data_printer.savePaired(
				cfg::get().output_dir + "etalon_paired_filtered",
				filtered_etalon_index);
		data_printer.savePaired(
				cfg::get().output_dir + "etalon_paired_corrected",
				final_etalon_paired_index);
		INFO("Everything saved");
	}
	INFO("Correction finished");
}

template<class Graph>
void GetAllDistances(const PairedInfoIndex<Graph>& paired_index,
		PairedInfoIndex<Graph>& result,
		const GraphDistanceFinder<Graph>& dist_finder) {
	for (auto iter = paired_index.begin(); iter != paired_index.end(); ++iter) {
		vector<PairInfo<EdgeId> > data = *iter;
		EdgeId first = data[0].first;
		EdgeId second = data[0].second;
		vector<size_t> forward = dist_finder.GetGraphDistances(first, second);
		//if (debug(first, second)) cout<<"i'm here"<<endl;
		for (size_t i = 0; i < forward.size(); i++)
			result.AddPairInfo(
					PairInfo < EdgeId
							> (data[0].first, data[0].second, forward[i], -10, 0.0),
					false);
	}
}

template<class Graph>
void CountAndSaveAllPaths(const Graph& g, const IdTrackHandler<Graph>& int_ids,
		const PairedInfoIndex<Graph>& paired_index) {
	paired_info_index all_paths(g);
	GetAllDistances<Graph>(
			paired_index,
			all_paths,
			GraphDistanceFinder<Graph>(g, *cfg::get().ds.IS, *cfg::get().ds.RL,
					size_t(*cfg::get().ds.is_var)));

	string dir_name = cfg::get().output_dir + "estimation_qual/";
	make_dir(dir_name);

	typename PrinterTraits<Graph>::Printer printer(g, int_ids);
	printer.savePaired(dir_name + "paths", all_paths);
}


void CountClusteredPairedInfoStats(const conj_graph_pack &gp,
		const PairedInfoIndex<Graph> &paired_index,
		const PairedInfoIndex<Graph> &clustered_index) {

	paired_info_index etalon_paired_index(gp.g);
	FillAndCorrectEtalonPairedInfo(etalon_paired_index, gp, paired_index,
			*cfg::get().ds.IS, *cfg::get().ds.RL, *cfg::get().ds.is_var, true);
	INFO(
			"Counting correlation between etalon and clustered paired infos statistics");
	std::map<std::pair<EdgeId, EdgeId>, vector<int>> my_index;
	std::map<int, int> gr;
	for (auto iter = etalon_paired_index.begin();
			iter != etalon_paired_index.end(); ++iter) {
		vector<PairInfo<EdgeId> > data = *iter;
		auto ppp = make_pair(data[0].first, data[0].second);
		for (size_t i = 0; i < data.size(); ++i)
			my_index[ppp].push_back(10 * data[i].d);
	}

	for (auto iter = clustered_index.begin(); iter != clustered_index.end();
			++iter) {
		vector<PairInfo<EdgeId> > data = *iter;
		auto ppp = make_pair(data[0].first, data[0].second);
		for (size_t i = 0; i < data.size(); ++i) {
			int min = 10000000;
			vector<int> etalon_pair = my_index[ppp];

			for (size_t j = 0; j < etalon_pair.size(); ++j) {
				if (min > abs(10 * data[i].d - etalon_pair[j]))
					min = abs(10 * data[i].d - etalon_pair[j]);
			}
			gr[min]++;
		}

	}

//    for (auto iter = gr.begin(); iter != gr.end(); ++iter){
//        cout << "Pavel " << (*iter).first / 100. << " "<< (*iter).second << endl;
//    }
	INFO("Counting clustered info stats");
	EdgeQuality<Graph> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
	EstimationQualityStat<Graph> estimation_stat(gp.g, gp.int_ids, edge_qual,
			paired_index, clustered_index, etalon_paired_index);
	estimation_stat.Count();
	estimation_stat.SaveStats();

	CountAndSaveAllPaths(gp.g, gp.int_ids, paired_index);

	INFO("Counting overall cluster stat")
	ClusterStat<Graph>(clustered_index).Count();
	INFO("Overall cluster stat")

	//	PairedInfoIndex<Graph> etalon_clustered_index;
	//	DistanceEstimator<Graph> estimator(g, etalon_paired_index, insert_size,
	//			max_read_length, cfg::get().de.delta,
	//			cfg::get().de.linkage_distance, cfg::get().de.max_distance);
	//	estimator.Estimate(etalon_clustered_index);

	//	PairedInfoIndex<Graph> filtered_clustered_index(g);
	//	PairInfoFilter<Graph> (g, 1000.).Filter(
	//			clustered_index/*etalon_clustered_index*/, filtered_clustered_index);
	INFO("Counting mate-pair transformation stat");
	MatePairTransformStat<Graph>(gp.g, /*filtered_*/clustered_index).Count();
	INFO("Mate-pair transformation stat counted");
	INFO("Clustered info stats counted");
}

void WriteToDotFile(const Graph &g,
		const omnigraph::GraphLabeler<Graph>& labeler, const string& file_name,
		string graph_name, Path<EdgeId> path1/* = Path<EdgeId> ()*/,
		Path<EdgeId> path2/* = Path<EdgeId> ()*/) {
	INFO("Writing graph '" << graph_name << "' to file " << file_name);
	omnigraph::WritePaired(g, labeler, file_name, graph_name, path1, path2);
	INFO("Graph '" << graph_name << "' written to file " << file_name);
}

void DetailedWriteToDot(const Graph &g,
		const omnigraph::GraphLabeler<Graph>& labeler, const string& file_name,
		string graph_name, Path<EdgeId> path1/* = Path<EdgeId> ()*/,
		Path<EdgeId> path2/* = Path<EdgeId> ()*/) {
	INFO("Writing graph '" << graph_name << "' to file " << file_name);
	omnigraph::WriteToFile(g, labeler, file_name, graph_name, path1, path2);
	INFO("Graph '" << graph_name << "' written to file " << file_name);
}

template<size_t k, class Graph>
Path<typename Graph::EdgeId> FindGenomePath(const Sequence& genome,
		const Graph& g, const EdgeIndex<k + 1, Graph>& index) {
	SimpleSequenceMapper<k + 1, Graph> srt(g, index);
	return srt.MapSequence(genome);
}

template<size_t k, class Graph>
MappingPath<typename Graph::EdgeId> FindGenomeMappingPath(
		const Sequence& genome, const Graph& g,
		const EdgeIndex<k + 1, Graph>& index,
		const KmerMapper<k + 1, Graph>& kmer_mapper) {
	ExtendedSequenceMapper<k + 1, Graph> srt(g, index, kmer_mapper);
	return srt.MapSequence(genome);
}

template<size_t k, class gp_t>
map<typename gp_t::graph_t::EdgeId, string> GraphColoring(const gp_t& gp) {
	return PathColorer<typename gp_t::graph_t>(
			gp.g,
			FindGenomeMappingPath<k>(gp.genome, gp.g, gp.index, gp.kmer_mapper).simple_path(),
			FindGenomeMappingPath<k>(!gp.genome, gp.g, gp.index, gp.kmer_mapper).simple_path()).ColorPath();
}

template<size_t k>
void ProduceInfo(const Graph& g, const EdgeIndex<k + 1, Graph>& index,
const omnigraph::GraphLabeler<Graph>& labeler, const Sequence& genome,
const string& file_name, const string& graph_name) {
	CountStats<k>(g, index, genome);
	Path<typename Graph::EdgeId> path1 = FindGenomePath<k>(genome, g, index);
	Path<typename Graph::EdgeId> path2 = FindGenomePath<k>(!genome, g, index);
	WriteToDotFile(g, labeler, file_name, graph_name, path1, path2);
}

template<size_t k>
void ProduceNonconjugateInfo(NCGraph& g, const EdgeIndex<k + 1, NCGraph>& index
		, const string& genome,
		const string& work_tmp_dir, const string& graph_name,
		const IdTrackHandler<NCGraph> &IdTrackLabelerResolved) {
	CountStats<k>(g, index, genome);
	//	omnigraph::WriteSimple( file_name, graph_name, g, IdTrackLabelerResolved);
	//	omnigraph::WriteSimple( work_tmp_dir, graph_name, g, IdTrackLabelerResolved);

}

void WriteGraphComponentsAlongGenome(const Graph& g,
		const IdTrackHandler<Graph>& int_ids,
		const EdgeIndex<K + 1, Graph>& index,
		const KmerMapper<K + 1, Graph>& kmer_mapper,
		const GraphLabeler<Graph>& labeler, const Sequence& genome,
		const string& folder, const string &file_name,
		size_t split_edge_length) {

	INFO("Writing graph components along genome");

	typedef MappingPath<EdgeId> map_path_t;

	map_path_t path1 = FindGenomeMappingPath<K>(genome, g, index, kmer_mapper);
	map_path_t path2 = FindGenomeMappingPath<K>(!genome, g, index, kmer_mapper);

	make_dir(folder);
	WriteComponentsAlongGenome(g, labeler, folder + file_name,
			split_edge_length, path1, path2);

	INFO("Writing graph components along genome finished");
}

//todo refactoring needed: use graph pack instead!!!
void WriteGraphComponentsAlongContigs(const Graph& g,
		const EdgeIndex<K + 1, Graph>& index,
		const KmerMapper<K + 1, Graph>& kmer_mapper,
		const GraphLabeler<Graph>& labeler, 
        const Sequence& genome,
		const string& folder,
		size_t split_edge_length) {

	INFO("Writing graph components along contigs");

	//typedef MappingPath<EdgeId> map_path_t;

	//typedef graph_pack<ConjugateDeBruijnGraph, K> gp_t;
	io::EasyReader contigs_to_thread(cfg::get().pos.contigs_to_analyze, true);
	contigs_to_thread.reset();

	NewExtendedSequenceMapper<K + 1, Graph> mapper(g, index, kmer_mapper);

	MappingPath<EdgeId> path1 = FindGenomeMappingPath<K>(genome, g, index, kmer_mapper);
	MappingPath<EdgeId> path2 = FindGenomeMappingPath<K>(!genome, g, index, kmer_mapper);


    io::SingleRead read;
    while (!contigs_to_thread.eof()) {
        contigs_to_thread >> read;
        make_dir(folder + read.name());
        WriteComponentsAlongPath(g, labeler, folder + read.name() + "/" + "g.dot"
                , /*split_edge_length*/400, mapper.MapSequence(read.sequence())
                , Path<Graph::EdgeId>(), Path<Graph::EdgeId>(), true);
                //, path1.simple_path(), path2.simple_path(), true);
    }
	INFO("Writing graph components along contigs finished");
}

void WriteKmerComponent(conj_graph_pack &gp,
		const omnigraph::GraphLabeler<Graph>& labeler, const string& folder,
		const Path<Graph::EdgeId>& path1, const Path<Graph::EdgeId>& path2,
		Seq<K + 1> const& kp1mer) {
	VERIFY(gp.index.contains(kp1mer));
	EdgeNeighborhoodFinder<Graph> splitter(gp.g, gp.index.get(kp1mer).first, 50,
			*cfg::get().ds.IS);
	ComponentSizeFilter<Graph> filter(gp.g, *cfg::get().ds.IS, 2);
	PathColorer<Graph> colorer(gp.g, path1, path2);
	WriteComponents<Graph>(gp.g, splitter, filter, folder + "kmer.dot",
			*DefaultColorer(gp.g, path1, path2), labeler);
}

optional<Seq<K + 1>> FindCloseKP1mer(const conj_graph_pack &gp,
		size_t genome_pos) {
	static const size_t magic_const = 200;
	for (size_t diff = 0; diff < magic_const; diff++) {
		for (int dir = -1; dir <= 1; dir += 2) {
			size_t pos = genome_pos + dir * diff;
			Seq < K + 1 > kp1mer = gp.kmer_mapper.Substitute(
					Seq < K + 1 > (gp.genome, pos));
			if (gp.index.contains(kp1mer))
				return optional<Seq<K + 1>>(kp1mer);
		}
	}
	return none;
}

void ProduceDetailedInfo(conj_graph_pack &gp,
		const omnigraph::GraphLabeler<Graph>& labeler, const string& folder,
		const string& file_name, const string& graph_name,
		info_printer_pos pos) {
	auto it = cfg::get().info_printers.find(pos);
	VERIFY(it != cfg::get().info_printers.end());

	const debruijn_config::info_printer & config = it->second;

	if (config.print_stats) {
		INFO("Printing statistics for " << details::info_printer_pos_name(pos));
		CountStats<K>(gp.g, gp.index, gp.genome);
	}

	typedef Path<Graph::EdgeId> path_t;
	path_t path1;
	path_t path2;

	if (config.detailed_dot_write || config.write_components
			|| !config.components_for_kmer.empty()
			|| config.write_components_along_genome
			|| config.write_components_along_contigs || config.save_full_graph
			|| !config.components_for_genome_pos.empty()) {
		path1 = FindGenomeMappingPath<K>(gp.genome, gp.g, gp.index,
				gp.kmer_mapper).simple_path();
		path2 = FindGenomeMappingPath<K>(!gp.genome, gp.g, gp.index,
				gp.kmer_mapper).simple_path();
//		path1 = FindGenomePath<K>(gp.genome, gp.g, gp.index);
//		path2 = FindGenomePath<K>(!gp.genome, gp.g, gp.index);
		make_dir(folder);
	}

	if (config.detailed_dot_write) {
		make_dir(folder + "error_loc/");
		DetailedWriteToDot(gp.g, labeler, folder + "error_loc/" + file_name,
				graph_name, path1, path2);
	}

	if (config.write_components) {
		make_dir(folder + "components/");
		size_t threshold = 500; //cfg::get().ds.IS ? *cfg::get().ds.IS : 250;
		WriteComponents(gp.g, threshold, folder + "components/" + file_name,
				*DefaultColorer(gp.g, path1, path2), labeler);
	}

	if (!config.components_for_kmer.empty()) {
		make_dir(folder + "kmer_loc/");
		WriteKmerComponent(gp, labeler, folder + "kmer_loc/", path1, path2,
				Seq < K + 1 > (config.components_for_kmer.c_str()));
	}

	if (config.write_components_along_genome) {
		make_dir(folder + "along_genome/");
		size_t threshold = 500; //cfg::get().ds.IS ? *cfg::get().ds.IS : 250;
		WriteGraphComponentsAlongGenome(gp.g, gp.int_ids, gp.index,
				gp.kmer_mapper, labeler, gp.genome, folder, "along_genome/",
				threshold);
	}

	if (config.write_components_along_contigs) {
		make_dir(folder + "along_contigs/");
		WriteGraphComponentsAlongContigs(gp.g, gp.index,
				gp.kmer_mapper, labeler, gp.genome, folder + "along_contigs/",
				*cfg::get().ds.IS);
	}

	if (config.save_full_graph) {
		make_dir(folder + "full_graph_save/");
		ConjugateDataPrinter<Graph> printer(gp.g, gp.int_ids);
		PrintGraphPack(folder + "full_graph_save/graph", printer, gp);
	}

	if (!config.components_for_genome_pos.empty()) {
		string pos_loc_folder = folder + "pos_loc/";
		make_dir(pos_loc_folder);
		vector<string> positions;
		boost::split(positions, config.components_for_genome_pos,
				boost::is_any_of(" ,"), boost::token_compress_on);
		for (auto it = positions.begin(); it != positions.end(); ++it) {
			optional < Seq < K + 1 >> close_kp1mer = FindCloseKP1mer(gp,
					boost::lexical_cast<int>(*it));
			if (close_kp1mer) {
				string locality_folder = pos_loc_folder + *it + "/";
				make_dir(locality_folder);
				WriteKmerComponent(gp, labeler, locality_folder, path1, path2,
						*close_kp1mer);
			} else {
				WARN(
						"Failed to find genome kp1mer close to the one at position "
								<< *it << " in the graph. Which is " << Seq
								< K + 1
								> (gp.genome, boost::lexical_cast<int>(*it)));
			}
		}
	}
}

struct detail_info_printer {
	detail_info_printer(conj_graph_pack &gp,
			const omnigraph::GraphLabeler<Graph>& labeler, const string& folder,
			const string& file_name)

	:
			folder_(folder), func_(
					bind(&ProduceDetailedInfo, ref(gp), ref(labeler), _3,
							file_name, _2, _1)) {
	}

	void operator()(info_printer_pos pos,
			string const& folder_suffix = "") const {
		string pos_name = details::info_printer_pos_name(pos);
		func_(
				pos,
				pos_name,
				(fs::path(folder_) / (pos_name + folder_suffix)).string()
						+ "/");
	}

private:
	string folder_;
	boost::function<void(info_printer_pos, string const&, string const&)> func_;
};

template<size_t k>
void WriteGraphComponents(const Graph& g, const EdgeIndex<k + 1, Graph>& index,
const GraphLabeler<Graph>& labeler, const Sequence& genome,
const string& folder, const string &file_name,
size_t split_edge_length) {
	make_dir(folder);
	WriteComponents(
			g,
			split_edge_length,
			folder + file_name,
			*DefaultColorer(FindGenomePath<k>(genome, g, index),
					FindGenomePath<k>(!genome, g, index)), labeler);

}

string ConstructComponentName(string file_name, size_t cnt) {
	stringstream ss;
	ss << cnt;
	string res = file_name;
	res.insert(res.length(), ss.str());
	return res;
}

template<class graph_pack>
int PrintGraphComponents(const string& file_name, graph_pack& gp,
		size_t split_edge_length, PairedInfoIndex<Graph> &clustered_index) {
	LongEdgesInclusiveSplitter<Graph> inner_splitter(gp.g, split_edge_length);
	ComponentSizeFilter<Graph> checker(gp.g, split_edge_length, 2);
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, checker);
	size_t cnt = 1;
	while (!splitter.Finished() && cnt <= 1000) {
		string component_name = ConstructComponentName(file_name, cnt).c_str();
		auto component = splitter.NextComponent();
		PrintWithClusteredIndex(component_name, gp, component.begin(),
				component.end(), clustered_index);
		cnt++;
	}
	return (cnt - 1);
}

void OutputContigs(NonconjugateDeBruijnGraph& g,
		const string& contigs_output_filename) {
	INFO("Outputting contigs to " << contigs_output_filename);
	osequencestream_cov oss(contigs_output_filename);
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		oss << g.coverage(*it);
		oss << g.EdgeNucls(*it);
	}
	DEBUG("Contigs written");
}

void OutputContigs(ConjugateDeBruijnGraph& g,
		const string& contigs_output_filename) {
	INFO("Outputting contigs to " << contigs_output_filename);
	osequencestream_cov oss(contigs_output_filename);
	set<ConjugateDeBruijnGraph::EdgeId> edges;
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		if (edges.count(*it) == 0) {
			oss << g.coverage(*it);
			oss << g.EdgeNucls(*it);
			edges.insert(g.conjugate(*it));
		}
		//		oss << g.EdgeNucls(*it);
	}
	DEBUG("Contigs written");
}

void OutputSingleFileContigs(NonconjugateDeBruijnGraph& g,
		const string& contigs_output_dir) {
	INFO("Outputting contigs to " << contigs_output_dir);
	int n = 0;
	make_dir(contigs_output_dir);
	char n_str[20];
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		sprintf(n_str, "%d.fa", n);

		osequencestream oss(contigs_output_dir + n_str);

		//		osequencestream oss(contigs_output_dir + "tst.fasta");
		oss << g.EdgeNucls(*it);
		n++;
	}
	DEBUG("SingleFileContigs written");
}

void OutputSingleFileContigs(ConjugateDeBruijnGraph& g,
		const string& contigs_output_dir) {
	INFO("Outputting contigs to " << contigs_output_dir);
	int n = 0;
	make_dir(contigs_output_dir);
	char n_str[20];
	set<ConjugateDeBruijnGraph::EdgeId> edges;
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		if (edges.count(*it) == 0) {
			sprintf(n_str, "%d.fa", n);
			edges.insert(g.conjugate(*it));
			osequencestream oss(contigs_output_dir + n_str);
			oss << g.EdgeNucls(*it);
			n++;
		}
	}
	DEBUG("SingleFileContigs(Conjugate) written");
}

void tSeparatedStats(conj_graph_pack& gp, const Sequence& contig,
		PairedInfoIndex<conj_graph_pack::graph_t> &ind) {
	typedef omnigraph::PairInfo<EdgeId> PairInfo;

	MappingPath<Graph::EdgeId> m_path1 = FindGenomeMappingPath<K>(contig, gp.g,
			gp.index, gp.kmer_mapper);

	map<Graph::EdgeId, vector<pair<int, int>>> inGenomeWay;
	int CurI = 0;
	int gaps = 0;
	for (size_t i = 0; i < m_path1.size(); i++) {
		bool new_edge_added = false;
		EdgeId ei = m_path1[i].first;
		MappingRange mr = m_path1[i].second;
		int start = mr.initial_range.start_pos - mr.mapped_range.start_pos;
		if (inGenomeWay.find(ei) == inGenomeWay.end()) {
			vector<pair<int, int>> tmp;
			tmp.push_back(make_pair(CurI, start));
			inGenomeWay[ei] = tmp;
			CurI++;
			new_edge_added = true;
			DEBUG(
					"Edge " << gp.int_ids.str(ei) << " num " << CurI << " pos "
							<< start);
		} else {
			if (m_path1[i - 1].first == ei) {
				if (abs(
						start
								- inGenomeWay[ei][(inGenomeWay[ei].size() - 1)].second)
						> 50) {
					inGenomeWay[ei].push_back(make_pair(CurI, start));
					CurI++;
					new_edge_added = true;
					DEBUG(
							"Edge " << gp.int_ids.str(ei) << " num " << CurI
									<< " pos " << start);
				}
			} else {
				inGenomeWay[ei].push_back(make_pair(CurI, start));
				CurI++;
				new_edge_added = true;
				DEBUG(
						"Edge " << gp.int_ids.str(ei) << " num " << CurI
								<< " pos " << start);
			}
		}
		if (new_edge_added && (i > 0)) {
			if (gp.g.EdgeStart(ei) != gp.g.EdgeEnd(m_path1[i - 1].first)) {
				gaps++;
			}
		}
	}
	INFO(
			"Totaly " << CurI << " edges in genome path, with " << gaps
					<< "not adjacent conequences");
	vector<int> stats(10);
	vector<int> stats_d(10);
	int PosInfo = 0;
	int AllignedPI = 0;
	int ExactDPI = 0;
	int OurD = *cfg::get().ds.IS - *cfg::get().ds.RL;
	for (auto p_iter = ind.begin(), p_end_iter = ind.end();
			p_iter != p_end_iter; ++p_iter) {
		vector<PairInfo> pi = *p_iter;
		for (size_t j = 0; j < pi.size(); j++) {
			EdgeId left_edge = pi[j].first;
			EdgeId right_edge = pi[j].second;
			int dist = pi[j].d;
			if (dist < 0.001)
				continue;
			int best_d = 100;
			int best_t = 0;
			PosInfo++;
			DEBUG(
					"PairInfo " << gp.int_ids.str(left_edge) << " -- "
							<< gp.int_ids.str(right_edge) << " dist " << dist);
			bool ExactOnD = false;
			for (size_t left_i = 0; left_i < inGenomeWay[left_edge].size();
					left_i++)
				for (size_t right_i = 0;
						right_i < inGenomeWay[right_edge].size(); right_i++) {
					if (best_d
							> abs(
									inGenomeWay[right_edge][right_i].second
											- inGenomeWay[left_edge][left_i].second
											- dist)) {
						best_d = abs(
								inGenomeWay[right_edge][right_i].second
										- inGenomeWay[left_edge][left_i].second
										- dist);
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
	}
	INFO(
			"Total positive pair info " << PosInfo << " alligned to genome "
					<< AllignedPI << " with exact distance " << ExactDPI);
	INFO(
			"t-separated stats Alligneg: 1 - " << stats[1] << " 2 - "
					<< stats[2] << " 3 - " << stats[3] << " 4 - " << stats[4]
					<< " >4 - " << stats[5]);
	INFO(
			"t-separated stats Exact: 1 - " << stats_d[1] << " 2 - "
					<< stats_d[2] << " 3 - " << stats_d[3] << " 4 - "
					<< stats_d[4] << " >4 - " << stats[5]);
}

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
//		Process(read.sequence(), read.name());
		MappingPath<EdgeId> path = mapper_.MapRead(read);
		const string& name = read.name();
		int cur_pos = 0;
		DEBUG(
				"Contig " << name << " mapped on " << path.size()
						<< " fragments.");
		for (size_t i = 0; i < path.size(); i++) {
			EdgeId ei = path[i].first;
			MappingRange mr = path[i].second;
			int len = mr.mapped_range.end_pos - mr.mapped_range.start_pos;
			if (i > 0)
				if (path[i - 1].first != ei)
					if (g_.EdgeStart(ei) != g_.EdgeEnd(path[i - 1].first)) {
						DEBUG(
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
			edge_pos_.AddEdgePosition(ei, mr.initial_range.start_pos + 1,
					mr.initial_range.end_pos, name);
			cur_pos += len;
		}
	}

private:
	DECL_LOGGER("PosFiller");
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
void FillPos(gp_t& gp,
		io::IReader<io::SingleRead>& stream) {
	typedef typename gp_t::graph_t Graph;
	typedef NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> Mapper;
	Mapper mapper(gp.g, gp.index, gp.kmer_mapper);
	FillPos<Graph, Mapper>(gp.g, mapper, gp.edge_pos, stream);
}

template<class gp_t>
void FillPos(gp_t& gp, const Sequence& s, const string& name) {
	typedef typename gp_t::graph_t Graph;
	typedef NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> Mapper;
	Mapper mapper(gp.g, gp.index, gp.kmer_mapper);
	PosFiller<Graph, Mapper>(gp.g, mapper, gp.edge_pos).Process(s, name);
}

//todo refactor!!!
class IdSettingReaderWrapper: public io::DelegatingReaderWrapper<io::SingleRead> {
	typedef io::DelegatingReaderWrapper<io::SingleRead> base;
	size_t next_id_;
public:
	IdSettingReaderWrapper(io::IReader<io::SingleRead>& reader, size_t start_id = 0) :
			base(reader), next_id_(start_id) {

	}

	/* virtual */
	IdSettingReaderWrapper& operator>>(io::SingleRead& read) {
		this->reader() >> read;
		read.ChangeName(ToString(next_id_++));
		return *this;
	}
};

class PrefixAddingReaderWrapper: public io::DelegatingReaderWrapper<io::SingleRead> {
	typedef io::DelegatingReaderWrapper<io::SingleRead> base;
	string prefix_;
public:
	PrefixAddingReaderWrapper(io::IReader<io::SingleRead>& reader,
			const string& prefix) :
			base(reader), prefix_(prefix) {

	}

	/* virtual */
	PrefixAddingReaderWrapper& operator>>(io::SingleRead& read) {
		this->reader() >> read;
		read.ChangeName(prefix_ + read.name());
		return *this;
	}
};

//deprecated, todo remove usages!!!
template<class gp_t>
void FillPos(gp_t& gp, const string& contig_file, int start_contig_id) {
//	typedef typename gp_t::Graph::EdgeId EdgeId;
	INFO("Threading large contigs");
	io::Reader irs(contig_file);
	for (int c = start_contig_id; !irs.eof(); c++) {
		io::SingleRead read;
		irs >> read;
		DEBUG("Contig #" << c << ", length: " << read.size());
		if (!read.IsValid()) {
			WARN("Attention: contig #" << c << " contains Ns");
			continue;
		}
		Sequence contig = read.sequence();
		if (contig.size() < 1500000) {
			//		continue;
		}
		FillPos(gp, contig, ToString(c));
	}
}

////template<size_t k>
////deprecated, todo remove usages
//void FillPos(conj_graph_pack& gp, const Sequence& genome) {
//	FillPos(gp, genome, 0);
//}

template<size_t k>
void OutputWrongContigs(conj_graph_pack& gp, size_t bound,
		const string &file_name) {
	OutputWrongContigs<k>(gp.g, gp.index, gp.genome, bound, file_name);
}

template<size_t k>
void OutputWrongContigs(Graph& g, EdgeIndex<k + 1, Graph>& index,
const Sequence& genome, size_t bound, const string &file_name) {
	SimpleSequenceMapper<k + 1, Graph> sequence_mapper(g, index);
	Path<EdgeId> path1 = sequence_mapper.MapSequence(Sequence(genome));
	Path<EdgeId> path2 = sequence_mapper.MapSequence(!Sequence(genome));
	set<EdgeId> path_set;
	path_set.insert(path1.begin(), path1.end());
	path_set.insert(path2.begin(), path2.end());
	osequencestream os((cfg::get().output_dir + "/" + file_name).c_str());
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		if (path_set.count(*it) == 0 && g.length(*it) > 1000) {
			const Sequence &nucls = g.EdgeNucls(*it);
			os << nucls;
		}
	}
}

/*//		Graph& g, const EdgeIndex<k + 1, Graph>& index,
 //		const Sequence& genome, EdgesPositionHandler<Graph>& edgesPos, KmerMapper<k + 1, Graph>& kmer_mapper)
 {
 Path<typename Graph::EdgeId> path1 = FindGenomePath<K> (genome, gp.g, gp.index);
 int CurPos = 0;
 for (auto it = path1.sequence().begin(); it != path1.sequence().end(); ++it) {
 EdgeId ei = *it;
 gp.edge_pos.AddEdgePosition(ei, CurPos + 1, CurPos + g.length(ei));
 CurPos += g.length(ei);
 }

 CurPos = 0;
 Path<typename Graph::EdgeId> path2 = FindGenomePath<k> (!genome, g, index);
 for (auto it = path2.sequence().begin(); it != path2.sequence().end(); ++it) {
 CurPos -= g.length(*it);
 }

 for (auto it = path2.sequence().begin(); it != path2.sequence().end(); ++it) {
 EdgeId ei = *it;
 edgesPos.AddEdgePosition(ei, CurPos, CurPos + g.length(ei) - 1);
 CurPos += g.length(ei);
 }

 }
 */

template<class Graph>
size_t Nx(Graph &g, double percent){
	size_t sum_edge_length = 0;
	vector<size_t> lengths;
	for (auto iterator = g.SmartEdgeBegin(); !iterator.IsEnd(); ++iterator) {
		lengths.push_back(g.length(*iterator));
		sum_edge_length += g.length(*iterator);
	}
	sort(lengths.begin(), lengths.end());
	double len_perc = (1 - percent * 0.01) * (sum_edge_length);
	for(size_t i = 0; i < lengths.size(); i++){
		if (lengths[i] >= len_perc)
			return lengths[i];
		else
			len_perc -= lengths[i];
	}
	return 0;
}




}
