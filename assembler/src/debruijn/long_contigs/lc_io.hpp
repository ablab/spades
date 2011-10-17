/*
 * lc_io.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: andrey
 */

#ifndef LC_IO_HPP_
#define LC_IO_HPP_

#include "io/reader.hpp"
#include "io/parser.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/cutting_reader_wrapper.hpp"
#include "io/filtering_reader_wrapper.hpp"
#include "path_utils.hpp"

namespace long_contigs {

using namespace debruijn_graph;

enum EdgeType { TOTALY_ISOLATED, SELF_LOOP, HAS_NEIGHBOUR, REGULAR };

class PairedInfoSimpleSymmertrizer {
private:
	Graph &g;

	int FindInfoByDistance(omnigraph::PairedInfoIndex<Graph>::PairInfos pairs, int d) {
		for (int i = 0; i < (int) pairs.size(); ++i) {
			if (d >= rounded_d(pairs[i]) - 1 && d <= rounded_d(pairs[i]) + 1) {
				return i;
			}
		}
		INFO("Can not find info with distance " << d);
		return -1;
	}

public:
	PairedInfoSimpleSymmertrizer(Graph& g_): g(g_) {
	}

	void MakeSymmetricInfo(PairedInfoIndexLibrary& lib) {
		INFO("Making symmetric info");
		for (auto edge1 = g.SmartEdgeBegin(); !edge1.IsEnd(); ++edge1) {
			for (auto edge2 = g.SmartEdgeBegin(); !edge2.IsEnd(); ++edge2) {

				std::vector<PairInfo<EdgeId>> pairs = lib.pairedInfoIndex->GetEdgePairInfo(*edge1, *edge2);
				std::vector<PairInfo<EdgeId>> conjPairs = lib.pairedInfoIndex->GetEdgePairInfo(g.conjugate(*edge2), g.conjugate(*edge1));

				if (pairs.size() != conjPairs.size()) {
					INFO("Paired info count is not the same in conjugate");
				}

				for (auto iter1 = pairs.begin(); iter1 != pairs.end(); ++iter1) {
					int d = rounded_d(*iter1);
					int conjD = d - g.length(*edge1) + g.length(*edge2);
					int i = FindInfoByDistance(conjPairs, conjD);

					if (i == -1) {
						INFO("Edges: " << g.length(*edge1) << ", " << g.length(*edge2)  << ". Distances: " << d << ", " << conjD);
						continue;
					}

					if (conjPairs[i].d != conjD) {
						INFO("Changing distances");
						conjPairs[i].d = conjD;

						std::vector<PairInfo<EdgeId>> conjSymPairs = lib.pairedInfoIndex->GetEdgePairInfo(g.conjugate(*edge1), g.conjugate(*edge2));
						int j = FindInfoByDistance(conjSymPairs, -conjD);
						if (j == -1) {
							INFO("SYMETIC CONJ! Edges: " << g.length(*edge1) << ", " << g.length(*edge2)  << ". Distances: " << d << ", " << - conjD);
						}
						conjSymPairs[j].d = - conjD;
					}

					if (conjPairs[i].weight != iter1->weight) {
						INFO("Changing weights");
						double newWeight = (conjPairs[i].weight + iter1->weight) / 2;
						iter1->weight = newWeight;
						conjPairs[i].weight = newWeight;


						std::vector<PairInfo<EdgeId>> symPairs = lib.pairedInfoIndex->GetEdgePairInfo(*edge2, *edge1);
						int j = FindInfoByDistance(symPairs, -d);
						if (j == -1) {
							INFO("SYMETIC! Edges: " << g.length(*edge1) << ", " << g.length(*edge2)  << ". Distances: " << d << ", " << -d);
						}
						symPairs[j].weight = newWeight;

						std::vector<PairInfo<EdgeId>> conjSymPairs = lib.pairedInfoIndex->GetEdgePairInfo(g.conjugate(*edge1), g.conjugate(*edge2));
						j = FindInfoByDistance(conjSymPairs, -conjD);
						if (j == -1) {
							INFO("SYMETIC CONJ! Edges: " << g.length(*edge1) << ", " << g.length(*edge2)  << ". Distances: " << d << ", " << - conjD);
						}
						conjSymPairs[j].weight = newWeight;
					}
				}
			}
		}
		INFO("Done")
	}

};


EdgeType ClassifyEdge(Graph& g, EdgeId e) {
	EdgeType type;
	if (g.IncomingEdgeCount(g.EdgeStart(e)) == 0 && g.IncomingEdgeCount(g.EdgeEnd(e)) == 1 &&
			g.OutgoingEdgeCount(g.EdgeStart(e)) == 1 && g.OutgoingEdgeCount(g.EdgeEnd(e)) == 0) {
		type = TOTALY_ISOLATED;
	} else if ((g.IncomingEdgeCount(g.EdgeStart(e)) == 0 && g.OutgoingEdgeCount(g.EdgeEnd(e)) == 0) &&
			 (g.IncomingEdgeCount(g.EdgeEnd(e)) > 1 || g.OutgoingEdgeCount(g.EdgeStart(e)) > 1)) {
		type = HAS_NEIGHBOUR;
	} else if (g.IncomingEdgeCount(g.EdgeStart(e)) == 1 && g.OutgoingEdgeCount(g.EdgeEnd(e)) == 1
			&& g.GetUniqueIncomingEdge(g.EdgeStart(e)) == e && g.GetUniqueOutgoingEdge(g.EdgeEnd(e)) == e) {

		type = SELF_LOOP;
	} else {
		type = REGULAR;
	}
	return type;
}

std::string ToStr(EdgeType t) {
	std::string res;
	switch(t) {
	case TOTALY_ISOLATED: {
		res = "Isolated";
		break;
	}
	case HAS_NEIGHBOUR: {
		res = "With neighbours";
		break;
	}
	case SELF_LOOP: {
		res = "Self looped";
		break;
	}
	case REGULAR: {
		res = "Regular";
		break;
	}
	}
	return res;
}

void PrintEdgesStats(std::map<EdgeType, size_t>& amount, std::map<EdgeType, size_t>& lengths, const std::string& title) {
	INFO(title);
	for (auto iter = amount.begin(); iter != amount.end(); ++iter) {
		INFO(ToStr(iter->first) << " edges: " << iter->second << " with total length " << lengths[iter->first]);
	}
}



void CheckPairedInfo(Graph& g, PairedInfoIndex<Graph>& index, size_t lenThreshold) {
	size_t total = 0, totalLength = 0;
	std::map<EdgeType, size_t> lengths;
	std::map<EdgeType, size_t> amount;

	std::map<EdgeType, size_t> noPairedInfoOtherLen;
	std::map<EdgeType, size_t> noPairedInfoOther;

	std::map<EdgeType, size_t> noPairedInfoOtherLongLen;
	std::map<EdgeType, size_t> noPairedInfoOtherLong;

	std::map<EdgeType, size_t> noPairedInfoLen;
	std::map<EdgeType, size_t> noPairedInfo;

	std::map<EdgeType, size_t> noPairedInfoLongLen;
	std::map<EdgeType, size_t> noPairedInfoLong;


	INFO("Checking paired info");
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		EdgeId e = *iter;
		size_t len = g.length(e);
		double dist = len / 2.0 + 1.0;
		auto pi = index.GetEdgeInfo(e);

		++total;
		totalLength += len;

		EdgeType et = ClassifyEdge(g, e);
		amount[et] += 1;
		lengths[et] += len;

		if (pi.size() == 0 || (pi.size() == 1 && index.GetEdgePairInfo(e,e).size() == 1 &&	math::le(std::abs(pi[0].d), dist))) {

			noPairedInfoOther[et] += 1;
			noPairedInfoOtherLen[et] += len;

			if (len > lenThreshold) {
				noPairedInfoOtherLong[et] += 1;
				noPairedInfoOtherLongLen[et] += len;
			}
		}

		if (pi.size() == 0 || (pi.size() == 1 && index.GetEdgePairInfo(e,e).size() == 1 && math::le(pi[0].weight, 0.1))) {
			noPairedInfo[et] += 1;
			noPairedInfoLen[et] += len;

			if (len > lenThreshold) {
				noPairedInfoLong[et] += 1;
				noPairedInfoLongLen[et] += len;
			}
		}


	}
	INFO("Total " << total << " with length " << totalLength);
	PrintEdgesStats(amount, lengths, "All edges");
	PrintEdgesStats(noPairedInfoOther, noPairedInfoOtherLen, "No paired info to other edges");
	PrintEdgesStats(noPairedInfoOtherLong, noPairedInfoOtherLongLen, "No paired info to other edges for edges longer than " + ToString(lenThreshold));
	PrintEdgesStats(noPairedInfo, noPairedInfoLen, "No paired info at all");
	PrintEdgesStats(noPairedInfoLong, noPairedInfoLongLen, "No paired info at all for edges longer than " + ToString(lenThreshold));
}

void LoadFromFile(std::string fileName, Graph* g,  IdTrackHandler<Graph>* conj_IntIds,	Sequence& sequence, KmerMapper<K+1, Graph> * mapper) {
	string input_dir = cfg::get().input_dir;
	string dataset = cfg::get().dataset_name;
	string genome_filename = input_dir
			+ cfg::get().ds.reference_genome;
	checkFileExistenceFATAL(genome_filename);
	int dataset_len = cfg::get().ds.LEN;

	typedef io::Reader<io::SingleRead> ReadStream;

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
		debruijn_graph::scanConjugateGraph(g, conj_IntIds, fileName, (PairedInfoIndex<Graph>*) 0,
			(EdgesPositionHandler<Graph> *) NULL, (PairedInfoIndex<Graph>*) 0, (PairedInfoIndex<Graph>*) 0, mapper);
	INFO("Graph read");
}

template<size_t k>
void AddEtalonInfo(const Graph& g, EdgeIndex<k+1, Graph>& index, KmerMapper<k+1, Graph>& mapper, const Sequence& genome, PairedInfoIndices& pairedInfos) {
	for (auto el = lc_cfg::get().etalon_libs.begin(); el != lc_cfg::get().etalon_libs.end(); ++el) {
		INFO("Generating info with read size " << el->read_size << ", insert size " << el->insert_size);

		pairedInfos.push_back(PairedInfoIndexLibrary(g, el->read_size, el->insert_size, 0, lc_cfg::get().es.etalon_distance_dev, new PairedInfoIndex<Graph>(g, 0)));
		FillEtalonPairedIndex<k> (*pairedInfos.back().pairedInfoIndex, g, index, mapper, el->insert_size, el->read_size, genome);
	}
}

template<size_t k>
void AddRealInfo(Graph& g, EdgeIndex<k+1, Graph>& index, IdTrackHandler<Graph>& conj_IntIds, PairedInfoIndices& pairedInfos, KmerMapper<k+1, Graph>& mapper,
		bool useNewMetrics) {
	PairedInfoSimpleSymmertrizer sym(g);

	for (auto rl = lc_cfg::get().real_libs.begin(); rl != lc_cfg::get().real_libs.end(); ++rl) {
		size_t insertSize = rl->insert_size;
		size_t readSize = rl->read_size;
		size_t delta = rl->is_delta;
		size_t var = rl->var;
		string dataset = cfg::get().dataset_name;
		pairedInfos.push_back(PairedInfoIndexLibrary(g, readSize, insertSize, delta, var, new PairedInfoIndex<Graph>(g, 0)));

		INFO("Reading additional info with read size " << readSize << ", insert size " << insertSize);

		if (rl->ds.precounted) {
			//Reading saved paired info
			DataScanner<Graph> dataScanner(g, conj_IntIds);
			dataScanner.loadPaired(rl->ds.precounted_path, *pairedInfos.back().pairedInfoIndex);
			CheckPairedInfo(g, *pairedInfos.back().pairedInfoIndex, insertSize - 2 * readSize + K);

			pairedInfos.back().raw = new PairedInfoIndex<Graph>(g, 0);
			if (!lc_cfg::get().paired_info_only) {
				//dataScanner.loadPaired(rl->ds.raw, *pairedInfos.back().raw);
		
				pairedInfos.back().has_advanced = rl->ds.has_advanced;
				if (rl->ds.has_advanced) {
					pairedInfos.back().advanced = new PairedInfoIndexLibrary(g, readSize, insertSize, delta, var, new PairedInfoIndex<Graph>(g, 0));
					DataScanner<Graph> advDataScanner(g, conj_IntIds);
					advDataScanner.loadPaired(rl->ds.advanced, *pairedInfos.back().advanced->pairedInfoIndex);
				} else {
					pairedInfos.back().advanced = 0;
				}
			}
		}
		else {
			//Reading paired info from fastq files
			string reads_filename1 = rl->ds.first;
			string reads_filename2 = rl->ds.second;
			checkFileExistenceFATAL(reads_filename1);
			checkFileExistenceFATAL(reads_filename2);

			typedef io::Reader<io::SingleRead> ReadStream;
			typedef io::Reader<io::PairedRead> PairedReadStream;
			typedef io::RCReaderWrapper<io::PairedRead> RCStream;
			typedef io::FilteringReaderWrapper<io::PairedRead> FilteringStream;

			PairedReadStream pairStream(std::pair<std::string,
									  std::string>(reads_filename1,
												   reads_filename2),
												   insertSize);

			FilteringStream filter_stream(pairStream);

			RCStream rcStream(filter_stream);

			if (useNewMetrics) {
				FillPairedIndexWithReadCountMetric<k, RCStream>(g, index, mapper, *pairedInfos.back().pairedInfoIndex, rcStream);
			} else {
				FillPairedIndexWithProductMetric<k, RCStream>(g, index, mapper, *pairedInfos.back().pairedInfoIndex, rcStream);
			}
		}
		INFO("Done");

		if (lc_cfg::get().syminfo) {
			sym.MakeSymmetricInfo(pairedInfos.back());
		}
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
			DistanceEstimator<Graph> estimator(g, *(lib->pairedInfoIndex), old_IDs, lib->insertSize, lib->readSize, cfg::get().de.delta,
					cfg::get().de.linkage_distance,
					cfg::get().de.max_distance);
			estimator.Estimate(clustered_index);

			dataPrinter.savePaired(fileName + "_clustered", clustered_index);

			PairedInfoIndex<Graph> clustered_index_(g);
            AdvancedDistanceEstimator<Graph> estimator_(g, *(lib->pairedInfoIndex), old_IDs,
            		lib->insertSize, lib->readSize, cfg::get().de.delta,
                    cfg::get().de.linkage_distance, cfg::get().de.max_distance, cfg::get().ade.threshold, cfg::get().ade.range_coeff, cfg::get().ade.delta_coeff, cfg::get().ade.cutoff, cfg::get().ade.minpeakpoints, cfg::get().ade.inv_density, cfg::get().ade.percentage, cfg::get().ade.derivative_threshold);
            estimator_.Estimate(clustered_index_);

            dataPrinter.savePaired(fileName + "_acl", clustered_index_);

		}
		dataPrinter.savePaired(fileName, *(lib->pairedInfoIndex));
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



void OutputPathsAsContigsNoComplement(Graph& g, std::vector<BidirectionalPath>& paths,
		const string& filename,
		std::set<int> notToPrint) {

	INFO("Writing contigs to " << filename);
	osequencestream oss(filename);

	for (int i = 0; i < (int) paths.size(); i += 2) {
		if (notToPrint.count(i) || notToPrint.count(i + 1) || paths[i].size() == 0) {
			continue;
		}

		oss.ptr = (void*) &paths[i];
		oss << PathToSequence(g, paths[i]);
	}

	INFO("Contigs written");
}


} // namespace long_contigs

#endif /* LC_IO_HPP_ */
