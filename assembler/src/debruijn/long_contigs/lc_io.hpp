/*
 * lc_io.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: andrey
 */

#ifndef LC_IO_HPP_
#define LC_IO_HPP_

#include "path_utils.hpp"

namespace long_contigs {

using namespace debruijn_graph;
using debruijn::K;

template<size_t k>
void LoadFromFile(std::string fileName, Graph* g,  IdTrackHandler<Graph>* conj_IntIds,	Sequence& sequence) {

	string input_dir = cfg::get().input_dir;
	string dataset = cfg::get().dataset_name;
	string genome_filename = input_dir
			+ cfg::get().reference_genome;
	checkFileExistenceFATAL(genome_filename);
	int dataset_len = cfg::get().ds.LEN;

	typedef io::Reader<io::SingleRead> ReadStream;
	typedef io::Reader<io::PairedRead> PairedReadStream;
	typedef io::RCReaderWrapper<io::PairedRead> RCStream;

	// read data ('genome')
	std::string genome;
	{
		ReadStream genome_stream(genome_filename);
		io::SingleRead full_genome;
		genome_stream >> full_genome;
		genome = full_genome.GetSequenceString().substr(0, dataset_len); // cropped
	}
	sequence = Sequence(genome);

	INFO("Reading graph");
	omnigraph::scanConjugateGraph(g, conj_IntIds, fileName);
	INFO("Graph read")
}

template<size_t k>
void AddEtalonInfo(const Graph& g, EdgeIndex<k+1, Graph>& index, const Sequence& genome, PairedInfoIndices& pairedInfos) {
	for (auto el = lc_cfg::get().etalon_libs.begin(); el != lc_cfg::get().etalon_libs.end(); ++el) {
		INFO("Generating info with read size " << el->read_size << ", insert size " << el->insert_size);

		pairedInfos.push_back(PairedInfoIndexLibrary(el->read_size, el->insert_size, lc_cfg::get().es.etalon_distance_dev, new PairedInfoIndex<Graph>(g, 0)));
		FillEtalonPairedIndex<k> (g, *pairedInfos.back().pairedInfoIndex, index, el->insert_size, el->read_size, genome);
	}
}

template<size_t k>
void AddRealInfo(Graph& g, EdgeIndex<k+1, Graph>& index, IdTrackHandler<Graph>& conj_IntIds, PairedInfoIndices& pairedInfos) {
	for (auto rl = lc_cfg::get().real_libs.begin(); rl != lc_cfg::get().real_libs.end(); ++rl) {
		size_t insertSize = rl->insert_size;
		size_t readSize = rl->read_size;
		size_t var = rl->var;
		string dataset = cfg::get().dataset_name;
		pairedInfos.push_back(PairedInfoIndexLibrary(readSize, insertSize, var, new PairedInfoIndex<Graph>(g, 0)));

		INFO("Reading additional info with read size " << readSize << ", insert size " << insertSize);

		if (rl->ds.precounted) {
			//Reading saved paired info
			DataScanner<Graph> dataScanner(g, conj_IntIds);
			dataScanner.loadPaired(rl->ds.precounted_path, *pairedInfos.back().pairedInfoIndex);
		}
		else {
			//Reading paired info from fastq files
			string reads_filename1 = rl->ds.first;
			string reads_filename2 = rl->ds.second;
			checkFileExistenceFATAL(reads_filename1);
			checkFileExistenceFATAL(reads_filename2);

			typedef io::Reader<io::PairedRead> PairedReadStream;
			typedef io::RCReaderWrapper<io::PairedRead> RCStream;

			PairedReadStream pairStream(std::pair<std::string,
									  std::string>(reads_filename1,
												   reads_filename2),
												   insertSize);

			RCStream rcStream(&pairStream);

			if (lc_cfg::get().use_new_metrics) {
				KmerMapper<k+1, Graph> mapper(g);
				FillPairedIndexWithReadCountMetric<k, RCStream>(g, index, mapper,*pairedInfos.back().pairedInfoIndex, rcStream);
			} else {
				FillPairedIndex<k, RCStream>(g, index, *pairedInfos.back().pairedInfoIndex, rcStream);
			}
		}
		INFO("Done");
	}
}

void SavePairedInfo(Graph& g, PairedInfoIndices& pairedInfos, IdTrackHandler<Graph>& old_IDs, const std::string& fileNamePrefix) {
	INFO("Saving paired info");
	DataPrinter<Graph> dataPrinter(g, old_IDs);
	for (auto lib = pairedInfos.begin(); lib != pairedInfos.end(); ++lib) {
		std::string fileName = fileNamePrefix + "IS" + ToString(lib->insertSize) + "_RS" + ToString(lib->readSize);
		INFO("Saving to " << fileName);

		if (lc_cfg::get().cluster_paired_info) {
			PairedInfoIndex<Graph> clustered_index(g);
			DistanceEstimator<Graph> estimator(g, *(lib->pairedInfoIndex), lib->insertSize, lib->readSize, cfg::get().de.delta,
					cfg::get().de.linkage_distance,
					cfg::get().de.max_distance);
			estimator.Estimate(clustered_index);

			dataPrinter.savePaired(fileName, clustered_index);
		} else {
			dataPrinter.savePaired(fileName, *(lib->pairedInfoIndex));
		}


	}
	INFO("Saved");
}

void SaveGraph(Graph& g, IdTrackHandler<Graph>& old_IDs, const std::string& fileName) {
	DataPrinter<Graph> dataPrinter(g, old_IDs);
	INFO("Saving graph file to " << fileName);
	dataPrinter.saveGraph(fileName);
}

void DeleteAdditionalInfo(PairedInfoIndices& pairedInfos) {
	while (pairedInfos.size() > 1) {
		delete pairedInfos.back().pairedInfoIndex;
		pairedInfos.pop_back();
	}
}

void PrintPairedInfo(Graph& g, PairedInfoIndex<Graph>& pairedInfo, IdTrackHandler<Graph>& intIds) {
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		for (auto iter2 = g.SmartEdgeBegin(); !iter2.IsEnd(); ++iter2) {
			omnigraph::PairedInfoIndex<Graph>::PairInfos pairs = pairedInfo.GetEdgePairInfo(*iter, *iter2);
			for (auto pair = pairs.begin(); pair != pairs.end(); ++pair) {
				INFO(intIds.ReturnIntId(pair->first) << " - " << intIds.ReturnIntId(pair->second) << " : " << pair->weight);
			}
		}
	}
}

//Convert path to contig sequence
Sequence PathToSequence(Graph& g, BidirectionalPath& path) {
	SequenceBuilder result;

	if (!path.empty()) {
		result.append(g.EdgeNucls(path[0]).Subseq(0, K));
	}
	for (auto edge = path.begin(); edge != path.end(); ++edge) {
		result.append(g.EdgeNucls(*edge).Subseq(K));
	}

	return result.BuildSequence();
}

//Output
void OutputPathsAsContigs(Graph& g, std::vector<BidirectionalPath> paths, const string& filename) {
	INFO("Writing contigs to " << filename);
	osequencestream oss(filename);
	for (auto path = paths.begin(); path != paths.end(); ++path ) {
		oss << PathToSequence(g, *path);
	}
	INFO("Contigs written");
}


//Output only one half of edges
void OutputContigsNoComplement(Graph& g, const std::string& filename) {
	std::set<EdgeId> filtered;
	FilterComlementEdges(g, filtered);

	INFO("Outputting contigs to " << filename);
	osequencestream oss(filename);
	for (auto it = filtered.begin(); it != filtered.end(); ++it) {
		oss << g.EdgeNucls(*it);
	}
	INFO("Contigs written");
}



void OutputPathsAsContigsNoComplement(Graph& g, std::vector<BidirectionalPath> paths, const string& filename) {
	INFO("Writing contigs to " << filename);
	osequencestream oss(filename);

	std::vector<BidirectionalPath> temp(paths.size());
	std::copy(paths.begin(), paths.end(), temp.begin());

	for (auto path = paths.begin(); path < paths.end(); path += 2) {
		oss << PathToSequence(g, *path);
	}
	INFO("Contigs written");
}


} // namespace long_contigs

#endif /* LC_IO_HPP_ */
