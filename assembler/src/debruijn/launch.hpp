/*
 * launch.hpp
 *
 *  Created on: May 6, 2011
 *      Author: sergey
 */

#ifndef LAUNCH_HPP_
#define LAUNCH_HPP_

#include "common/io/reader.hpp"
#include "common/io/rc_reader_wrapper.hpp"
#include "common/io/cutting_reader_wrapper.hpp"
#include "common/io/converting_reader_wrapper.hpp"
#include "visualization_utils.hpp"

//#include "debruijn_graph.hpp"
#include "paired_info.hpp"
#include "debruijn_graph_constructor.hpp"
#include "tip_clipper.hpp"
#include "bulge_remover.hpp"
#include "erroneous_connection_remover.hpp"
#include "coverage_handler.hpp"
#include "repeat_resolver.hpp"
#include "omni_tools.hpp"
#include "seq_map.hpp"
#include "ID_track_handler.hpp"
#include "edges_position_handler.hpp"
#include "read/osequencestream.hpp"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "new_debruijn.hpp"
#include "config.hpp"
#include "graphio.hpp"
#include "rectangleRepeatResolver.hpp"
#include "distance_estimation.hpp"
#include "one_many_contigs_enlarger.hpp"
#include <cstdlib>
//#include "dijkstra.hpp"

namespace debruijn_graph {

typedef ConjugateDeBruijnGraph Graph;
typedef Graph::EdgeId EdgeId;
typedef Graph::VertexId VertexId;
typedef NonconjugateDeBruijnGraph NCGraph;
using namespace omnigraph;

template<size_t k, class Graph>
void CountStats(Graph& g, const EdgeIndex<k + 1, Graph>& index,
const Sequence& genome) {
	INFO("Counting stats");
	StatCounter<Graph, k> stat(g, index, genome);
	stat.Count();
	INFO("Stats counted");
}

void CountPairedInfoStats(Graph &g, size_t insert_size, size_t max_read_length,
		PairedInfoIndex<Graph> &paired_index, PairedInfoIndex<Graph> &etalon_paired_index, const string &output_folder) {
	INFO("Counting paired info stats");
	EdgePairStat<Graph>(g, paired_index, output_folder).Count();

	//todo remove filtration if launch on etalon info is ok
	UniquePathStat<Graph>(g, etalon_paired_index, insert_size, max_read_length, 0.1,
			40.0).Count();
	UniqueDistanceStat<Graph>(etalon_paired_index).Count();
	INFO("Paired info stats counted");
}

void CountClusteredPairedInfoStats(Graph &g, size_t insert_size, size_t max_read_length,
		PairedInfoIndex<Graph> &paired_index, PairedInfoIndex<Graph> &etalon_paired_index, const string &output_folder) {
	INFO("Counting paired info stats");
	EstimationQualityStat<Graph>(paired_index, etalon_paired_index).Count();
	ClusterStat<Graph>(paired_index).Count();
	INFO("Paired info stats counted");
}

void WriteToDotFile(Graph &g, const string& file_name, string graph_name,
		Path<EdgeId> path1/* = Path<EdgeId> ()*/,
		Path<EdgeId> path2/* = Path<EdgeId> ()*/) {
	INFO("Writing graph '" << graph_name << "' to file " << file_name);
	omnigraph::WritePaired(file_name, graph_name, g, path1, path2);
	INFO("Graph '" << graph_name << "' written to file " << file_name);
}

void DetailedWriteToDot(Graph &g, const string& file_name, string graph_name,
		Path<EdgeId> path1/* = Path<EdgeId> ()*/,
		Path<EdgeId> path2/* = Path<EdgeId> ()*/) {
	INFO("Writing graph '" << graph_name << "' to file " << file_name);
	omnigraph::WriteToFile(file_name, graph_name, g, path1, path2);
	INFO("Graph '" << graph_name << "' written to file " << file_name);
}

template<size_t k>
Path<typename Graph::EdgeId> FindGenomePath(const Sequence& genome,
		const Graph& g, const EdgeIndex<k + 1, Graph>& index) {
	SimpleSequenceMapper<k, Graph> srt(g, index);
	return srt.MapSequence(genome);
}

template<size_t k>
void ProduceInfo(Graph& g, const EdgeIndex<k + 1, Graph>& index,
const Sequence& genome, const string& file_name,
const string& graph_name) {
	CountStats<k>(g, index, genome);
	Path<typename Graph::EdgeId> path1 = FindGenomePath<k>(genome, g, index);
	Path<typename Graph::EdgeId> path2 = FindGenomePath<k>(!genome, g, index);
	WriteToDotFile(g, file_name, graph_name, path1, path2);
}

template<size_t k>
void FillEdgesPos(Graph& g, const EdgeIndex<k + 1, Graph>& index,
const Sequence& genome, EdgesPositionHandler<Graph>& edgesPos) {
	Path<typename Graph::EdgeId> path1 = FindGenomePath<k>(genome, g, index);
	int CurPos = 0;
	for (auto it = path1.sequence().begin(); it != path1.sequence().end();
			++it) {
		EdgeId ei = *it;
		edgesPos.AddEdgePosition(ei, CurPos + 1, CurPos + g.length(ei));
		CurPos += g.length(ei);
	}
	CurPos = 1000000000;
	Path<typename Graph::EdgeId> path2 = FindGenomePath<k>(!genome, g, index);
	for (auto it = path2.sequence().begin(); it != path2.sequence().end();
			++it) {
		EdgeId ei = *it;
		edgesPos.AddEdgePosition(ei, CurPos + 1, CurPos + g.length(ei));
		CurPos += g.length(ei);
	}
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

template<size_t k>
void ProduceDetailedInfo(Graph& g, const EdgeIndex<k + 1, Graph>& index,
const Sequence& genome, const string& folder, const string& file_name,
const string& graph_name) {
	CountStats<k>(g, index, genome);
	Path<typename Graph::EdgeId> path1 = FindGenomePath<k>(genome, g, index);
	Path<typename Graph::EdgeId> path2 = FindGenomePath<k>(!genome, g, index);

	mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	DetailedWriteToDot(g, folder + file_name, graph_name, path1, path2);
}

template<size_t k>
void WriteGraphComponents(Graph& g, const EdgeIndex<k + 1, Graph>& index,
const Sequence& genome, const string& folder, const string &file_name,
const string &graph_name, size_t split_edge_length) {
	Path<typename Graph::EdgeId> path1 = FindGenomePath<k>(genome, g, index);
	Path<typename Graph::EdgeId> path2 = FindGenomePath<k>(!genome, g, index);
	mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	WriteComponents(folder + file_name, graph_name, g, split_edge_length, path1,
			path2);

}

string ConstructComponentName(string file_name, size_t cnt) {
	stringstream ss;
	ss << cnt;
	string res = file_name;
	res.insert(res.length() - 4, ss.str());
	return res;
}

template<class Graph>
void PrintGraphComponents(const string& file_name, Graph& g,
		size_t split_edge_length, IdTrackHandler<Graph> &old_IDs,
		PairedInfoIndex<Graph> &paired_index,
		EdgesPositionHandler<Graph> &edges_positions) {
	LongEdgesSplitter<Graph> inner_splitter(g, split_edge_length);
	ComponentSizeFilter<Graph> checker(g, split_edge_length);
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, checker);
	size_t cnt = 1;
	while (!splitter.Finished() && cnt <= 100) {
		string component_name = ConstructComponentName(file_name, cnt).c_str();
		auto component = splitter.NextComponent();
		EdgeVertexFilter<Graph> *filter = new EdgeVertexFilter<Graph>(g,
				component);
		printGraph(g, old_IDs, component_name, paired_index, edges_positions,
				filter);
		delete filter;
		cnt++;
	}

}

template<class Graph>
void ClipTips(Graph &g) {
	INFO("-----------------------------------------");
	INFO("Clipping tips");
	TipComparator<Graph> comparator(g);
	size_t max_tip_length = CONFIG.read<size_t>("tc_max_tip_length");
	size_t max_coverage = CONFIG.read<size_t>("tc_max_coverage");
	double max_relative_coverage = CONFIG.read<double>(
			"tc_max_relative_coverage");
	TipClipper<Graph, TipComparator<Graph>> tc(g, comparator, max_tip_length,
			max_coverage, max_relative_coverage);
	tc.ClipTips();
	INFO("Clipping tips finished");
}

void RemoveBulges(Graph &g) {
	INFO("-----------------------------------------");
	INFO("Removing bulges");
	double max_coverage = CONFIG.read<double>("br_max_coverage");
	double max_relative_coverage = CONFIG.read<double>(
			"br_max_relative_coverage");
	double max_delta = CONFIG.read<double>("br_max_delta");
	double max_relative_delta = CONFIG.read<double>("br_max_relative_delta");
	size_t max_length_div_K = CONFIG.read<int>("br_max_length_div_K");
	SimplePathCondition<Graph> simple_path_condition(g);
	BulgeRemover<Graph, SimplePathCondition<Graph>> bulge_remover(g,
			max_length_div_K * g.k(), max_coverage, max_relative_coverage,
			max_delta, max_relative_delta, simple_path_condition);
	bulge_remover.RemoveBulges();
	INFO("Bulges removed");
}

void RemoveBulges2(NCGraph &g) {
	INFO("-----------------------------------------");
	INFO("Removing bulges");
	double max_coverage = CONFIG.read<double>("br_max_coverage");
	double max_relative_coverage = CONFIG.read<double>(
			"br_max_relative_coverage");
	double max_delta = CONFIG.read<double>("br_max_delta");
	double max_relative_delta = CONFIG.read<double>("br_max_relative_delta");
	size_t max_length_div_K = CONFIG.read<int>("br_max_length_div_K");
	TrivialCondition<NCGraph> trivial_condition;
	BulgeRemover<NCGraph, TrivialCondition<NCGraph>> bulge_remover(g,
			max_length_div_K * g.k(), max_coverage, max_relative_coverage,
			max_delta, max_relative_delta, trivial_condition);
	bulge_remover.RemoveBulges();
	INFO("Bulges removed");
}

template<class Graph>
void RemoveLowCoverageEdges(Graph &g) {
	INFO("-----------------------------------------");
	INFO("Removing low coverage edges");
	double max_coverage = CONFIG.read<double>("ec_max_coverage");
	int max_length_div_K = CONFIG.read<int>("ec_max_length_div_K");
	LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
			max_length_div_K * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges(g);
	INFO("Low coverage edges removed");
}

template<class Graph>
void RemoveLowCoverageEdgesForResolver(Graph &g) {
	INFO("-----------------------------------------");
	INFO("Removing low coverage edges");
	double max_coverage = CONFIG.read<double>("ec_max_coverage");
//	int max_length_div_K = CONFIG.read<int> ("ec_max_length_div_K");
	LowCoverageEdgeRemover<Graph> erroneous_edge_remover(10000000 * g.k(),
			max_coverage);
	erroneous_edge_remover.RemoveEdges(g);
	INFO("Low coverage edges removed");
}
template<class Graph>
void ResolveRepeats(Graph &g, IdTrackHandler<Graph> &old_IDs,
		PairedInfoIndex<Graph> &info, EdgesPositionHandler<Graph> &edges_pos,
		Graph &new_graph, IdTrackHandler<Graph> &new_IDs,
		EdgesPositionHandler<Graph> &edges_pos_new,
		const string& output_folder) {
	INFO("-----------------------------------------");
	INFO("Resolving primitive repeats");
	RepeatResolver<Graph> repeat_resolver(g, old_IDs, 0, info, edges_pos,
			new_graph, new_IDs, edges_pos_new);
	mkdir((output_folder).c_str(),
			S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	repeat_resolver.ResolveRepeats(output_folder);
	INFO("Primitive repeats resolved");
}
template<size_t k, class ReadStream, class Graph>
void MapPairedReads(Graph &g, ReadStream& stream,
		EdgeIndex<k + 1, Graph>& index) {
	INFO("-----------------------------------------");
	stream.reset();
	INFO("Threading reads");
	int quantity = 0, map_quantity = 0;

	SingleReadMapper<k, Graph, ReadStream> rm(g, index, stream);
	while (!stream.eof()) {
		vector<EdgeId> res = rm.GetContainingEdges();
		quantity++;
		map_quantity += res.size();
	}INFO(2 * quantity <<" reads_threaded,to "<< map_quantity <<" edges");
}

template<size_t k, class ReadStream>
void FillPairedIndex(Graph &g, PairedInfoIndex<Graph>& paired_info_index,
		ReadStream& stream, EdgeIndex<k + 1, Graph>& index) {
	INFO("-----------------------------------------");
	stream.reset();
	INFO("Counting paired info");
	PairedIndexFiller<Graph, k, ReadStream> pif(g, index, stream);
	pif.FillIndex(paired_info_index);
	INFO("Paired info counted");
}

template<size_t k>
void FillEtalonPairedIndex(Graph &g, PairedInfoIndex<Graph>& paired_info_index,
		EdgeIndex<k + 1, Graph>& index, size_t insert_size, size_t read_length,
		const Sequence& genome) {
	INFO("-----------------------------------------");
	INFO("Counting etalon paired info");
	EtalonPairedInfoCounter<k, Graph> etalon_paired_info_counter(g, index,
			insert_size, read_length, insert_size * 0.1);
	etalon_paired_info_counter.FillEtalonPairedInfo(genome, paired_info_index);
	INFO("Paired info counted");
}

template<size_t k, class ReadStream>
void FillCoverage(Graph& g/*CoverageHandler<Graph> coverage_handler*/,
		ReadStream& stream, EdgeIndex<k + 1, Graph>& index) {
	typedef SimpleSequenceMapper<k, Graph> ReadThreader;
	INFO("-----------------------------------------");
	stream.reset();
	INFO("Counting coverage");
	ReadThreader read_threader(g, index);
	//todo temporary solution!
	g.FillCoverage<ReadStream, ReadThreader>(stream, read_threader);
	//	coverage_handler.FillCoverage<k, ReadStream> (stream, index);
	INFO("Coverage counted");
}

template<size_t k, class ReadStream>
void ConstructGraph(Graph& g, EdgeIndex<k + 1, Graph>& index,
ReadStream& stream) {
	typedef SeqMap<k + 1, typename Graph::EdgeId> DeBruijn;
	INFO("-----------------------------------------");
	INFO("Constructing DeBruijn graph");
	DeBruijn& debruijn = index.inner_index();
	INFO("Filling DeBruijn graph");
	debruijn.Fill(stream);
	INFO("DeBruijn graph constructed");

	INFO("Condensing graph");
	DeBruijnGraphConstructor<k, Graph> g_c(debruijn);
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");
}

template<size_t k, class ReadStream>
void ConstructGraphWithCoverage(Graph& g, EdgeIndex<k + 1, Graph>& index,
/*CoverageHandler<Graph>& coverage_handler,*/ReadStream& stream) {
	ConstructGraph<k, ReadStream>(g, index, stream);
	FillCoverage<k, ReadStream>(g/*coverage_handler*/, stream, index);
}

template<size_t k, class PairedReadStream>
/*CoverageHandler<Graph>& coverage_handler,*/
void ConstructGraphWithPairedInfo(Graph& g, EdgeIndex<k + 1, Graph>& index,
PairedInfoIndex<Graph>& paired_index, PairedReadStream& stream) {
	typedef io::ConvertingReaderWrapper UnitedStream;
	UnitedStream united_stream(&stream);
	ConstructGraphWithCoverage<k, UnitedStream>(g, index/*, coverage_handler*/,
			united_stream);
	FillPairedIndex<k, PairedReadStream>(g, paired_index, stream, index);
}

template<size_t k, class PairedReadStream>
void ConstructGraphWithEtalonPairedInfo(Graph& g, EdgeIndex<k + 1, Graph>& index
		, PairedInfoIndex<Graph>& paired_index,
		PairedReadStream& stream, size_t insert_size, size_t read_length,
		const Sequence& genome) {
	typedef io::ConvertingReaderWrapper UnitedStream;
	UnitedStream united_stream(&stream);
	ConstructGraphWithCoverage<k, UnitedStream>(g, index/*, coverage_handler*/,
			united_stream);
	FillEtalonPairedIndex<k>(g, paired_index, index, insert_size, read_length,
			genome);
}

template<class Graph>
void printGraph(Graph & g, IdTrackHandler<Graph> &old_IDs,
		const string &file_name, PairedInfoIndex<Graph> &paired_index,
		EdgesPositionHandler<Graph> &edges_positions,
		EdgeVertexFilter<Graph> *filter) {
	DataPrinter<Graph> dataPrinter(g, old_IDs, filter);
	dataPrinter.saveGraph(file_name);
	dataPrinter.saveEdgeSequences(file_name);
	dataPrinter.saveCoverage(file_name);
	dataPrinter.savePaired(file_name, paired_index);
	dataPrinter.savePositions(file_name, edges_positions);
}

template<class Graph>
void printGraph(Graph & g, IdTrackHandler<Graph> &old_IDs,
		const string &file_name, PairedInfoIndex<Graph> &paired_index,
		EdgesPositionHandler<Graph> &edges_positions) {
	DataPrinter<Graph> dataPrinter(g, old_IDs);
	dataPrinter.saveGraph(file_name);
	dataPrinter.saveEdgeSequences(file_name);
	dataPrinter.saveCoverage(file_name);
	dataPrinter.savePaired(file_name, paired_index);
	dataPrinter.savePositions(file_name, edges_positions);
}

template<class Graph>
void printGraph(Graph & g, IdTrackHandler<Graph> &old_IDs,
		const string &file_name, PairedInfoIndex<Graph> &paired_index) {
	DataPrinter<Graph> dataPrinter(g, old_IDs);
	dataPrinter.saveGraph(file_name);
	dataPrinter.saveEdgeSequences(file_name);
	dataPrinter.saveCoverage(file_name);
	dataPrinter.savePaired(file_name, paired_index);

}

template<class Graph>
void scanNCGraph(Graph & g, IdTrackHandler<Graph> &new_IDs,
		const string &file_name, PairedInfoIndex<Graph>& paired_index,
		EdgesPositionHandler<Graph> &edges_positions) {
	DataScanner<Graph> dataScanner(g, new_IDs);
	dataScanner.loadNonConjugateGraph(file_name, true);
	dataScanner.loadCoverage(file_name);
	dataScanner.loadPaired(file_name, paired_index);
	dataScanner.loadPositions(file_name, edges_positions);
}

template<class Graph>
void scanNCGraph(Graph & g, IdTrackHandler<Graph> &new_IDs,
		const string &file_name, PairedInfoIndex<Graph>& paired_index) {
	DataScanner<Graph> dataScanner(g, new_IDs);
	dataScanner.loadNonConjugateGraph(file_name, true);
	dataScanner.loadCoverage(file_name);
	dataScanner.loadPaired(file_name, paired_index);
}

template<class Graph>
void scanConjugateGraph(Graph & g, IdTrackHandler<Graph> &new_IDs,
		const string &file_name, PairedInfoIndex<Graph>& paired_index) {
	DataScanner<Graph> dataScanner(g, new_IDs);
	dataScanner.loadConjugateGraph(file_name, true);
	dataScanner.loadCoverage(file_name);
	//	dataScanner.loadPaired(file_name, paired_index);
}

template<size_t k>
void SimplifyGraph(Graph& g, EdgeIndex<k + 1, Graph>& index,
size_t iteration_count, const Sequence& genome,
const string& output_folder) {
	INFO("-----------------------------------------");
	INFO("Graph simplification started");

	ProduceDetailedInfo<k>(g, index, genome,
			output_folder + "before_simplification/", "graph.dot",
			"non_simplified_graph");
	for (size_t i = 0; i < iteration_count; i++) {
		INFO("-----------------------------------------");
		INFO("Iteration " << i);

		ClipTips(g);
		ProduceDetailedInfo<k>(g, index, genome,
				output_folder + "tips_clipped_" + ToString(i) + "/",
				"graph.dot", "no_tip_graph");

		RemoveBulges(g);
		ProduceDetailedInfo<k>(g, index, genome,
				output_folder + "bulges_removed_" + ToString(i) + "/",
				"graph.dot", "no_bulge_graph");

		RemoveLowCoverageEdges(g);
		ProduceDetailedInfo<k>(g, index, genome,
				output_folder + "erroneous_edges_removed_" + ToString(i) + "/",
				"graph.dot", "no_erroneous_edges_graph");
	}INFO("Graph simplification finished");
}

template<class Graph>
void OutputContigs(Graph& g, const string& contigs_output_filename) {
	INFO("-----------------------------------------");
	INFO("Outputting contigs to " << contigs_output_filename);
	osequencestream oss(contigs_output_filename);
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		oss << g.EdgeNucls(*it);
	}INFO("Contigs written");
}

template<class Graph>
void OutputSingleFileContigs(Graph& g, const string& contigs_output_dir) {
	INFO("-----------------------------------------");
	INFO("Outputting contigs to " << contigs_output_dir);
	int n = 0;
	mkdir(contigs_output_dir.c_str(),
			S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	char n_str[20];
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		sprintf(n_str, "%d.fa", n);

		osequencestream oss(contigs_output_dir + n_str);

//		osequencestream oss(contigs_output_dir + "tst.fasta");
		oss << g.EdgeNucls(*it);
		n++;
	}INFO("Contigs written");
}

template<size_t k, class ReadStream>
void DeBruijnGraphTool(ReadStream& stream, const Sequence& genome,
		bool paired_mode, bool rectangle_mode, bool etalon_info_mode,
		bool from_saved, size_t insert_size, size_t max_read_length,
		const string& output_folder, const string& work_tmp_dir) {
	INFO("Edge graph construction tool started");
	INFO("Paired mode: " << (paired_mode ? "Yes" : "No"));
	INFO("Etalon paired info mode: " << (etalon_info_mode ? "Yes" : "No"))INFO(
			"From file: " << (from_saved ? "Yes" : "No"))
	mkdir(work_tmp_dir.c_str(),
			S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	Graph g(k);
	EdgeIndex<k + 1, Graph> index(g);
	IdTrackHandler<Graph> IntIds(g);
	EdgesPositionHandler<Graph> EdgePos(g);
	EdgesPosGraphLabeler<Graph> EdgePosLab(g, EdgePos);
	// if it's not paired_mode, then it'll be just unused variable -- takes O(1) to initialize from graph
//	PairedInfoIndex<Graph> paired_index(g, 5);
	PairedInfoIndex<Graph> paired_index(g, 0);
	PairedInfoIndex<Graph> etalon_paired_index(g, 0);

	if (!from_saved) {

		if (paired_mode) {
			if (etalon_info_mode) {
				ConstructGraphWithEtalonPairedInfo<k, ReadStream>(g, index,
						paired_index, stream, insert_size, max_read_length,
						genome);
			} else {
				ConstructGraphWithPairedInfo<k, ReadStream>(g, index,
						paired_index, stream);
			}
			FillEtalonPairedIndex<k>(g, etalon_paired_index, index, insert_size,
					max_read_length, genome);
		} else {
			typedef io::ConvertingReaderWrapper UnitedStream;
			UnitedStream united_stream(&stream);
			ConstructGraphWithCoverage<k, UnitedStream>(g, index,
					united_stream);
		}

		ProduceInfo<k>(g, index, genome, output_folder + "edge_graph.dot",
				"edge_graph");

		FillEdgesPos<k>(g, index, genome, EdgePos);
		omnigraph::WriteSimple(output_folder + "before_simplification_pos.dot",
				"no_repeat_graph", g, EdgePosLab);

		printGraph(g, IntIds, output_folder + "first_graph", paired_index,
				EdgePos);

		SimplifyGraph<k>(g, index, 3, genome, output_folder);
		MapPairedReads<k, ReadStream, Graph>(g, stream, index);

		ProduceInfo<k>(g, index, genome, output_folder + "simplified_graph.dot",
				"simplified_graph");

		WriteGraphComponents<k>(g, index, genome,
				output_folder + "graph_components" + "/", "graph.dot",
				"graph_component", insert_size);
		PrintGraphComponents(output_folder + "graph_components/graph", g,
				insert_size, IntIds, paired_index, EdgePos);

		if (paired_mode) {
			CountPairedInfoStats(g, insert_size, max_read_length, paired_index, etalon_paired_index,
					output_folder);
		}
	}
	//	if (paired_mode) {
	//		paired_index.OutputData(output_folder + "edges_dist.txt");
	//
	//		SimpleOfflineClusterer<Graph> clusterer(paired_index);
	//		PairedInfoIndex<Graph> clustered_paired_index(g);
	//		clusterer.cluster(clustered_paired_index);
	//	}

	omnigraph::WriteSimple(output_folder + "repeats_resolved_before_poslab.dot",
			"no_repeat_graph", g, EdgePosLab);
	omnigraph::WriteSimple(work_tmp_dir + "repeats_resolved_before_poslab.dot",
			"no_repeat_graph", g, EdgePosLab);
	PairedInfoIndex<Graph> clustered_index(g);

	printGraph(g, IntIds, output_folder + "repeats_resolved_before",
			paired_index, EdgePos);

	if (paired_mode) {
		DistanceEstimator<Graph> estimator(g, paired_index, insert_size,
				max_read_length, 10, 10, 75);
		estimator.Estimate(clustered_index);
		CountClusteredPairedInfoStats(g, insert_size, max_read_length, clustered_index, etalon_paired_index,
						output_folder);
	}

	if (paired_mode) {
		if (!from_saved) {
			INFO("before ResolveRepeats");
			RealIdGraphLabeler<Graph> IdTrackLabelerBefore(g, IntIds);
			omnigraph::WriteSimple(
					output_folder + "repeats_resolved_before.dot",
					"no_repeat_graph", g, IdTrackLabelerBefore);
			printGraph(g, IntIds, work_tmp_dir + "graph", clustered_index,
					EdgePos);
			printGraph(g, IntIds, output_folder + "graph", clustered_index,
					EdgePos);
		}

		NCGraph new_graph(k);
		IdTrackHandler<NCGraph> NewIntIds(new_graph, IntIds.MaxVertexId(),
				IntIds.MaxEdgeId());
		PairedInfoIndex<NCGraph> new_index(new_graph);
		EdgesPositionHandler<NCGraph> EdgePosBefore(new_graph);

		Graph conj_copy_graph(k);
		IdTrackHandler<Graph> conj_IntIds(conj_copy_graph, IntIds.MaxVertexId(),
				IntIds.MaxEdgeId());
		PairedInfoIndex<Graph> conj_copy_index(conj_copy_graph);
		/*
		 scanConjugateGraph(conj_copy_graph, conj_IntIds,
		 work_tmp_dir + "graph", conj_copy_index);
		 printGraph(conj_copy_graph, conj_IntIds, work_tmp_dir + "graph_copy",
		 conj_copy_index);
		 */
		scanNCGraph(new_graph, NewIntIds, work_tmp_dir + "graph", new_index,
				EdgePosBefore);

		RealIdGraphLabeler<NCGraph> IdTrackLabelerAfter(new_graph, NewIntIds);

		omnigraph::WriteSimple(
				work_tmp_dir + "repeats_resolved_nonconjugate_copy.dot",
				"no_repeat_graph", new_graph, IdTrackLabelerAfter);
		INFO("repeat resolved graph written");

		NonconjugateDeBruijnGraph resolved_graph(k);
		IdTrackHandler<NCGraph> Resolved_IntIds(resolved_graph);
		EdgesPositionHandler<NCGraph> EdgePosAfter(resolved_graph);

		DEBUG("New index size: "<< new_index.size());
		if (rectangle_mode) {
			void RectangleResolve(
					PairedInfoIndex<NonconjugateDeBruijnGraph>& index,
					NonconjugateDeBruijnGraph& graph,
					const string& work_tmp_dir, const string& output_folder);
			RectangleResolve(new_index, new_graph, work_tmp_dir, output_folder);
		}

		ResolveRepeats(new_graph, NewIntIds, new_index, EdgePosBefore,
				resolved_graph, Resolved_IntIds, EdgePosAfter,
				output_folder + "resolve/");

		RealIdGraphLabeler<NCGraph> IdTrackLabelerResolved(resolved_graph,
				Resolved_IntIds);

		omnigraph::WriteSimple(work_tmp_dir + "repeats_resolved_after.dot",
				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);
		omnigraph::WriteSimple(output_folder + "repeats_resolved_after.dot",
				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);
		EdgesPosGraphLabeler<NCGraph> EdgePosLAfterLab(resolved_graph,
				EdgePosAfter);

		omnigraph::WriteSimple(work_tmp_dir + "repeats_resolved_after_pos.dot",
				"no_repeat_graph", resolved_graph, EdgePosLAfterLab);
		omnigraph::WriteSimple(output_folder + "repeats_resolved_after_pos.dot",
				"no_repeat_graph", resolved_graph, EdgePosLAfterLab);

		ClipTips(resolved_graph);
		RemoveLowCoverageEdgesForResolver(resolved_graph);

		omnigraph::WriteSimple(
				work_tmp_dir + "repeats_resolved_after_und_cleared_pos.dot",
				"no_repeat_graph", resolved_graph, EdgePosLAfterLab);
		omnigraph::WriteSimple(
				output_folder + "repeats_resolved_after_und_cleared_pos.dot",
				"no_repeat_graph", resolved_graph, EdgePosLAfterLab);

		omnigraph::WriteSimple(
				work_tmp_dir + "repeats_resolved_und_cleared.dot",
				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);
		omnigraph::WriteSimple(
				output_folder + "repeats_resolved_und_cleared.dot",
				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);
//		one_many_contigs_enlarger<NCGraph> N50enlarger(resolved_graph);
//		N50enlarger.one_many_resolve();

		omnigraph::WriteSimple(
				output_folder
						+ "repeats_resolved_und_cleared_und_simplified.dot",
				"no_repeat_graph", resolved_graph, IdTrackLabelerResolved);
		INFO("repeat resolved grpah written");
		EdgeIndex<k + 1, NCGraph> aux_index(resolved_graph);

		//		SimplifyGraph<k>(resolved_graph, aux_index, 3, genome, output_folder);

		//CountStats<k, NCGraph> (resolved_graph, aux_index, genome);
		//		ProduceNonconjugateInfo<k> (resolved_graph, aux_index, genome, output_folder + "repeats_resolved.dot",
		//				"no_repeat_graph");
		//		ProduceNonconjugateInfo<k> (resolved_graph, aux_index, genome, work_tmp_dir + "repeats_resolved.dot",
		//				"no_repeat_graph");sss

		OutputContigs(resolved_graph, output_folder + "contigs.fasta");

		OutputSingleFileContigs(resolved_graph, output_folder + "consensus/");
		OutputContigs(new_graph,
				output_folder + "contigs_before_resolve.fasta");

	}
	if (!paired_mode)
		OutputContigs(g, output_folder + "contigs.fasta");
	INFO("Tool finished");
}

void RectangleResolve(PairedInfoIndex<NonconjugateDeBruijnGraph>& index,
		NonconjugateDeBruijnGraph& graph, const string& work_tmp_dir,
		const string& output_folder) {

	NonconjugateDeBruijnGraph resolvedGraph(graph.k());
	typedef NonconjugateDeBruijnGraph::EdgeId NCEdgeId;
	PairInfoIndexData<NCEdgeId> piid;
	for (auto iter = index.begin(); iter != index.end(); ++iter) {

		vector<PairInfo<NCEdgeId> > pi = *iter;
		for (size_t i = 0; i < pi.size(); ++i) {
			if (pi[i].d >= 0)
				piid.AddPairInfo(pi[i], 1);
		}
	}
	RectangleRepeatResolver<NonconjugateDeBruijnGraph> rectangleResolver(graph,
			piid, resolvedGraph, (size_t) 30);
	rectangleResolver.Process();

	ClipTips(resolvedGraph);
	RemoveLowCoverageEdges(resolvedGraph);
	IdTrackHandler<NCGraph> Resolved_IntIds(resolvedGraph);
	RealIdGraphLabeler<NCGraph> IdTrackLabelerResolved(resolvedGraph,
			Resolved_IntIds);

	ClipTips(resolvedGraph);
	RemoveLowCoverageEdges(resolvedGraph);
	EmptyGraphLabeler<NonconjugateDeBruijnGraph> emptyLabeler;

	omnigraph::WriteSimple(work_tmp_dir + "rectgraph.dot", "rectgraph",
			resolvedGraph, IdTrackLabelerResolved);
	INFO("rect graph written: " + work_tmp_dir + "rectgraph.dot");

	omnigraph::WriteSimple(work_tmp_dir + "before-rectgraph.dot",
			"before-rectgraph", graph, emptyLabeler);
	INFO("rect graph written: " + work_tmp_dir + "before-rectgraph.dot");

	OutputContigs(resolvedGraph, output_folder + "rectcontig.fasta");
	OutputContigs(graph, output_folder + "before-rectcontig.fasta");
}

}

#endif /* LAUNCH_HPP_ */
