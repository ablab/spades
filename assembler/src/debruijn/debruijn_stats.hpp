#pragma once

#include "omni/visualization_utils.hpp"
#include "statistics.hpp"
#include "new_debruijn.hpp"
#include "graphio.hpp"
#include "omni/edges_position_handler.hpp"
#include "omni/distance_estimation.hpp"
#include "omni/graph_component.hpp"
//#include <boost/filesystem.hpp>
#include "read/osequencestream.hpp"
#include "k.hpp"
#include <sys/types.h>
#include <sys/stat.h>

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
		Seq<k + 1> cur = genome_.start<k + 1>() >> 0;
		bool breaked = true;
		pair<EdgeId, size_t> cur_position;
		for (size_t cur_nucl = k; cur_nucl < genome_.size(); cur_nucl++) {
			cur = cur << genome_[cur_nucl];
			if (index_.containsInIndex(cur)) {
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
		}INFO("Genome mapped");
		INFO("Genome mapping results:");
		INFO(
				"Covered k+1-mers:" << covered_kp1mers << " of " << (genome_.size() - k) << " which is " << (100.0 * covered_kp1mers / (genome_.size() - k)) << "%");
		INFO(
				"Covered k+1-mers form " << break_number + 1 << " contigious parts");
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
	PairInfoFilter<Graph>(g, 40).Filter(paired_index, filtered_index);
	INFO("Counting paired info stats");
	EdgePairStat<Graph>(g, paired_index, output_folder).Count();

	//todo remove filtration if launch on etalon info is ok
	UniquePathStat<Graph>(g, filtered_index, cfg::get().ds.IS, cfg::get().ds.RL,
			0.1 * cfg::get().ds.IS).Count();
	UniqueDistanceStat<Graph>(etalon_paired_index).Count();
	INFO("Paired info stats counted");
}

void CountClusteredPairedInfoStats(const conj_graph_pack &gp,
		const PairedInfoIndex<Graph> &paired_index,
		const PairedInfoIndex<Graph> &clustered_index,
		const PairedInfoIndex<Graph> &etalon_paired_index,
		const string &output_folder) {
	INFO("Counting clustered info stats");

	EdgeQuality<Graph> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
	EstimationQualityStat<Graph> estimation_stat(gp.g, gp.int_ids, edge_qual, paired_index,
			clustered_index, etalon_paired_index);
	estimation_stat.Count();
	estimation_stat.SaveStats();

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

template<size_t k>
Path<typename Graph::EdgeId> FindGenomePath(const Sequence& genome,
		const Graph& g, const EdgeIndex<k + 1, Graph>& index) {
	SimpleSequenceMapper<k + 1, Graph> srt(g, index);
	return srt.MapSequence(genome);
}

template<size_t k>
MappingPath<EdgeId> NewFindGenomePath(const Sequence& genome, const Graph& g,
		const IdTrackHandler<Graph>& int_ids,
		const EdgeIndex<k + 1, Graph>& index,
		const KmerMapper<k + 1, Graph>& kmer_mapper) {
	ExtendedSequenceMapper<k + 1, Graph> srt(g, int_ids, index, kmer_mapper);
	return srt.MapSequence(genome);
}

template<size_t k>
MappingPath<typename Graph::EdgeId> FindGenomeMappingPath(
		const Sequence& genome, const Graph& g,
		const IdTrackHandler<Graph>& int_ids,
		const EdgeIndex<k + 1, Graph>& index,
		KmerMapper<k + 1, Graph>& kmer_mapper) {
	ExtendedSequenceMapper<k + 1, Graph> srt(g, int_ids, index, kmer_mapper);
	return srt.MapSequence(genome);
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
		const string &graph_name, size_t split_edge_length) {

    INFO("Writing graph components along genome");

    typedef MappingPath<EdgeId> map_path_t;

    map_path_t path1 = NewFindGenomePath<K>( genome, g, int_ids, index, kmer_mapper);
    map_path_t path2 = NewFindGenomePath<K>(!genome, g, int_ids, index, kmer_mapper);

	make_dir(folder);
	WriteComponentsAlongGenome(g, labeler, folder + file_name, graph_name, split_edge_length, path1, path2);

	INFO("Writing graph components along genome finished");
}

void WriteKmerComponent(
    conj_graph_pack &gp,
	const omnigraph::GraphLabeler<Graph>& labeler,
	const string& folder,
	const string& graph_name,
	Path<Graph::EdgeId> path1,
	Path<Graph::EdgeId> path2,
	Seq<K+1> const& kp1mer)
{
	KMerNeighborhoodFinder<Graph, K> splitter(gp.g, kp1mer, gp.index, 50, cfg::get().ds.IS);
	WriteComponents<Graph>(gp.g, labeler, folder + "kmer.dot", graph_name, cfg::get().ds.IS, splitter, path1, path2);
}

void ProduceDetailedInfo(
	conj_graph_pack &gp,
	const omnigraph::GraphLabeler<Graph>& labeler,
	const string& folder,
	const string& file_name,
	const string& graph_name,
	info_printer_pos pos)
{
    auto it = cfg::get().info_printers.find(pos);
    VERIFY(it != cfg::get().info_printers.end());

    debruijn_config::info_printer const& config = it->second;

    make_dir(folder);

    if (config.print_stats)
    {
        INFO("Printing statistics for " << details::info_printer_pos_name(pos));
        CountStats<K>(gp.g, gp.index, gp.genome);
    }

	typedef Path<Graph::EdgeId> path_t;
	path_t path1;
	path_t path2;

	if (config.detailed_dot_write           ||
	    config.write_components             ||
	    !config.components_for_kmer.empty() ||
	    config.write_components_along_genome)
	{
	    path1 = FindGenomePath<K>( gp.genome, gp.g, gp.index);
        path2 = FindGenomePath<K>(!gp.genome, gp.g, gp.index);
        make_dir(folder);
	}

	if (config.detailed_dot_write)
	    DetailedWriteToDot(gp.g, labeler, folder + file_name, graph_name, path1, path2);

	if (config.write_components)
	    WriteComponents(gp.g, labeler, folder + file_name, graph_name, cfg::get().ds.IS, path1, path2);

	if (!config.components_for_kmer.empty())
	    WriteKmerComponent(gp, labeler, folder, graph_name, path1, path2, Seq<K + 1>(config.components_for_kmer.c_str()));

	if (config.write_components_along_genome)
	    WriteGraphComponentsAlongGenome(
            gp.g,
            gp.int_ids,
            gp.index,
            gp.kmer_mapper,
            labeler,
            gp.genome,
            folder,
            file_name,
            "components_along_genome",
            cfg::get().ds.IS);

	if (config.save_full_graph) {
		ConjugateDataPrinter<Graph> printer(gp.g, gp.int_ids);
		PrintGraphPack(folder + "graph", printer, gp);
	}
}

struct detail_info_printer
{
    detail_info_printer(
        conj_graph_pack &gp,
        const omnigraph::GraphLabeler<Graph>& labeler,
        const string& folder,
        const string& file_name)

        : folder_   (folder)
        , func_     (bind(&ProduceDetailedInfo, ref(gp), ref(labeler), _3, file_name, _2, _1))
    {
    }

    void operator()(info_printer_pos pos, string const& folder_suffix = "") const
    {
        string pos_name = details::info_printer_pos_name(pos);
        func_(pos, pos_name, (fs::path(folder_) / (pos_name + folder_suffix)).string() + "/");
    }

private:
    string  folder_;
    boost::function<void(info_printer_pos, string const&, string const&)>  func_;
};

template<size_t k>
void WriteGraphComponents(const Graph& g, const EdgeIndex<k + 1, Graph>& index,
const GraphLabeler<Graph>& labeler, const Sequence& genome,
const string& folder, const string &file_name,
const string &graph_name, size_t split_edge_length) {
	Path<typename Graph::EdgeId> path1 = FindGenomePath<k>(genome, g, index);
	Path<typename Graph::EdgeId> path2 = FindGenomePath<k>(!genome, g, index);
	make_dir(folder);
	WriteComponents(g, labeler, folder + file_name, graph_name,
			split_edge_length, path1, path2);

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
		size_t split_edge_length,
		PairedInfoIndex<Graph> &clustered_index) {
	LongEdgesInclusiveSplitter<Graph> inner_splitter(gp.g, split_edge_length);
	ComponentSizeFilter<Graph> checker(gp.g, split_edge_length, 2);
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, checker);
	size_t cnt = 1;
	while (!splitter.Finished() && cnt <= 1000) {
		string component_name = ConstructComponentName(file_name, cnt).c_str();
		auto component = splitter.NextComponent();
		PrintWithClusteredIndex(component_name, gp, component.begin(), component.end(), clustered_index);
		cnt++;
	}
	return (cnt - 1);
}

void OutputContigs(NonconjugateDeBruijnGraph& g,
		const string& contigs_output_filename) {
	INFO("-----------------------------------------");
	INFO("Outputting contigs to " << contigs_output_filename);
	osequencestream_cov oss(contigs_output_filename);
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		oss << g.coverage(*it);
		oss << g.EdgeNucls(*it);
	}INFO("Contigs written");
}

void OutputContigs(ConjugateDeBruijnGraph& g,
		const string& contigs_output_filename) {
	INFO("-----------------------------------------");
	INFO("Outputting contigs to " << contigs_output_filename);
	osequencestream_cov oss(contigs_output_filename);
	set<ConjugateDeBruijnGraph::EdgeId> edges;
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		if (edges.find(*it) == edges.end()) {
			oss << g.coverage(*it);
			oss << g.EdgeNucls(*it);
			edges.insert(g.conjugate(*it));
		}
		//		oss << g.EdgeNucls(*it);
	}INFO("Contigs written");
}

void OutputSingleFileContigs(NonconjugateDeBruijnGraph& g,
		const string& contigs_output_dir) {
	INFO("-----------------------------------------");
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
	}INFO("SingleFileContigs written");
}

void OutputSingleFileContigs(ConjugateDeBruijnGraph& g,
		const string& contigs_output_dir) {
	INFO("-----------------------------------------");
	INFO("Outputting contigs to " << contigs_output_dir);
	int n = 0;
	make_dir(contigs_output_dir);
	char n_str[20];
	set<ConjugateDeBruijnGraph::EdgeId> edges;
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		if (edges.find(*it) == edges.end()) {
			sprintf(n_str, "%d.fa", n);
			edges.insert(g.conjugate(*it));
			osequencestream oss(contigs_output_dir + n_str);
			oss << g.EdgeNucls(*it);
			n++;
		}
	}INFO("SingleFileContigs(Conjugate) written");
}

void tSeparatedStats(conj_graph_pack& gp, const Sequence& contig,
		PairedInfoIndex<conj_graph_pack::graph_t> &ind) {
	typedef omnigraph::PairInfo<EdgeId> PairInfo;

	MappingPath<Graph::EdgeId> m_path1 = FindGenomeMappingPath<K>(contig, gp.g,
			gp.int_ids, gp.index, gp.kmer_mapper);

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
			DEBUG("Edge "<<gp.int_ids.str(ei)<< " num "<< CurI<<" pos "<<start);
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
							"Edge "<<gp.int_ids.str(ei)<< " num "<< CurI<<" pos "<<start);
				}
			} else {
				inGenomeWay[ei].push_back(make_pair(CurI, start));
				CurI++;
				new_edge_added = true;
				DEBUG(
						"Edge "<<gp.int_ids.str(ei)<< " num "<< CurI<<" pos "<<start);
			}
		}
		if (new_edge_added && (i > 0)) {
			if (gp.g.EdgeStart(ei) != gp.g.EdgeEnd(m_path1[i - 1].first)) {
				gaps++;
			}
		}
	}INFO(
			"Totaly "<<CurI<<" edges in genome path, with "<<gaps<< "not adjacent conequences");
	vector<int> stats(10);
	vector<int> stats_d(10);
	int PosInfo = 0;
	int AllignedPI = 0;
	int ExactDPI = 0;
	int OurD = cfg::get().ds.IS - cfg::get().ds.RL;
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
					"PairInfo "<< gp.int_ids.str(left_edge)<<" -- "<< gp.int_ids.str(right_edge)<<" dist "<<dist);
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
						DEBUG("best d "<<best_d);
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
			"Total positive pair info "<<PosInfo<< " alligned to genome "<< AllignedPI <<" with exact distance "<< ExactDPI);
	INFO(
			"t-separated stats Alligneg: 1 - "<<stats[1]<<" 2 - "<<stats[2]<<" 3 - "<<stats[3]<<" 4 - "<<stats[4]<<" >4 - "<<stats[5]);
	INFO(
			"t-separated stats Exact: 1 - "<<stats_d[1]<<" 2 - "<<stats_d[2]<<" 3 - "<<stats_d[3]<<" 4 - "<<stats_d[4]<<" >4 - "<<stats[5]);
}

//template<size_t k>
void FillEdgesPos(conj_graph_pack& gp, const Sequence& contig, int contigId)

//		Graph& g, const EdgeIndex<k + 1, Graph>& index,
//		const Sequence& genome, EdgesPositionHandler<Graph>& edgesPos, KmerMapper<k + 1, Graph>& kmer_mapper, int contigId)

		{
	MappingPath<Graph::EdgeId> m_path1 = FindGenomeMappingPath<K>(contig, gp.g,
			gp.int_ids, gp.index, gp.kmer_mapper);
	int CurPos = 0;
	DEBUG("Contig "<<contigId<< " maped on "<<m_path1.size()<<" fragments.");
	for (size_t i = 0; i < m_path1.size(); i++) {
		EdgeId ei = m_path1[i].first;
		MappingRange mr = m_path1[i].second;
		int len = mr.mapped_range.end_pos - mr.mapped_range.start_pos;
		if (i > 0)
			if (m_path1[i - 1].first != ei)
				if (gp.g.EdgeStart(ei) != gp.g.EdgeEnd(m_path1[i - 1].first)) {
					DEBUG(
							"Contig "<<contigId<<" maped on not adjacent edge. Position in contig is "<<m_path1[i-1].second.initial_range.start_pos+1<<"--" <<m_path1[i-1].second.initial_range.end_pos<< " and "<<mr.initial_range.start_pos+1<<"--" <<mr.initial_range.end_pos);
				}
		gp.edge_pos.AddEdgePosition(ei, mr.initial_range.start_pos + 1,
				mr.initial_range.end_pos, contigId);
		CurPos += len;
	}
	//CurPos = 1000000000;
}

//template<size_t k>
void FillEdgesPos(conj_graph_pack& gp, const string& contig_file,
		int start_contig_id)

//(Graph& g, const EdgeIndex<k + 1, Graph>& index,
//		const string& contig_file, EdgesPositionHandler<Graph>& edgesPos, KmerMapper<k + 1, Graph>& kmer_mapper, int start_contig_id)

		{
	INFO("Threading large contigs");
	io::Reader<io::SingleRead> irs(contig_file);
	for (int c = start_contig_id; !irs.eof(); c++) {
		io::SingleRead read;
		irs >> read;
		DEBUG("Contig #" << c << ", length: " << read.size());
		read.ClearQuality();
		if (!read.IsValid()) {
			WARN("Attention: contig #" << c << " contains Ns");
			continue;
		}
		Sequence contig = read.sequence();
		if (contig.size() < 1500000) {
			//		continue;
		}
		FillEdgesPos(gp, contig, c);
	}
}

//template<size_t k>
void FillEdgesPos(conj_graph_pack& gp, const Sequence& genome) {
	FillEdgesPos(gp, genome, 0);
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

}
