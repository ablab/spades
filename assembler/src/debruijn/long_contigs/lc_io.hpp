/*
 * lc_io.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: andrey
 */

#ifndef LC_IO_HPP_
#define LC_IO_HPP_

#include "../graph_pack.hpp"
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
		DEBUG("Can not find info with distance " << d);
		return -1;
	}

public:
	PairedInfoSimpleSymmertrizer(Graph& g_): g(g_) {
	}

	void MakeSymmetricInfo(PairedInfoIndexLibrary& lib) {
		DEBUG("Making symmetric info");
		for (auto edge1 = g.SmartEdgeBegin(); !edge1.IsEnd(); ++edge1) {
			for (auto edge2 = g.SmartEdgeBegin(); !edge2.IsEnd(); ++edge2) {

				std::vector<PairInfo<EdgeId>> pairs = lib.pairedInfoIndex->GetEdgePairInfo(*edge1, *edge2);
				std::vector<PairInfo<EdgeId>> conjPairs = lib.pairedInfoIndex->GetEdgePairInfo(g.conjugate(*edge2), g.conjugate(*edge1));

				if (pairs.size() != conjPairs.size()) {
					DEBUG("Paired info count is not the same in conjugate");
				}

				for (auto iter1 = pairs.begin(); iter1 != pairs.end(); ++iter1) {
					int d = rounded_d(*iter1);
					int conjD = d - g.length(*edge1) + g.length(*edge2);
					int i = FindInfoByDistance(conjPairs, conjD);

					if (i == -1) {
						DEBUG("Edges: " << g.length(*edge1) << ", " << g.length(*edge2)  << ". Distances: " << d << ", " << conjD);
						continue;
					}

					if (conjPairs[i].d != conjD) {
						DEBUG("Changing distances");
						conjPairs[i].d = conjD;

						std::vector<PairInfo<EdgeId>> conjSymPairs = lib.pairedInfoIndex->GetEdgePairInfo(g.conjugate(*edge1), g.conjugate(*edge2));
						int j = FindInfoByDistance(conjSymPairs, -conjD);
						if (j == -1) {
							DEBUG("SYMETIC CONJ! Edges: " << g.length(*edge1) << ", " << g.length(*edge2)  << ". Distances: " << d << ", " << - conjD);
						}
						conjSymPairs[j].d = - conjD;
					}

					if (conjPairs[i].weight != iter1->weight) {
						DEBUG("Changing weights");
						double newWeight = (conjPairs[i].weight + iter1->weight) / 2;
						iter1->weight = newWeight;
						conjPairs[i].weight = newWeight;


						std::vector<PairInfo<EdgeId>> symPairs = lib.pairedInfoIndex->GetEdgePairInfo(*edge2, *edge1);
						int j = FindInfoByDistance(symPairs, -d);
						if (j == -1) {
							DEBUG("SYMETIC! Edges: " << g.length(*edge1) << ", " << g.length(*edge2)  << ". Distances: " << d << ", " << -d);
						}
						symPairs[j].weight = newWeight;

						std::vector<PairInfo<EdgeId>> conjSymPairs = lib.pairedInfoIndex->GetEdgePairInfo(g.conjugate(*edge1), g.conjugate(*edge2));
						j = FindInfoByDistance(conjSymPairs, -conjD);
						if (j == -1) {
							DEBUG("SYMETIC CONJ! Edges: " << g.length(*edge1) << ", " << g.length(*edge2)  << ". Distances: " << d << ", " << - conjD);
						}
						conjSymPairs[j].weight = newWeight;
					}
				}
			}
		}
		DEBUG("Done")
	}

};


EdgeType ClassifyEdge(const Graph& g, EdgeId e) {
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



void CheckPairedInfo(const Graph& g, const PairedInfoIndex<Graph>& index, size_t lenThreshold) {
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

Sequence load_genome() {
	string genome_filename = lc_cfg::get().ds.reference_genome;
	std::string genome;
	if (genome_filename.length() > 0) {
		checkFileExistenceFATAL(genome_filename);
		io::Reader genome_stream(genome_filename);
		io::SingleRead full_genome;
		genome_stream >> full_genome;
		genome = full_genome.GetSequenceString().substr(0, lc_cfg::get().ds.LEN);
	}
	if (io::SingleRead::IsValid(genome)) {
		return Sequence(genome);
	}
	else {
		INFO("Reference genome (" + genome_filename + ") have non-ACGT characters. Skipping it.");
		return Sequence();
	}
}

void LoadFromFile(std::string fileName, Graph& g, IdTrackHandler<Graph>& intIds, KmerMapper<K+1, Graph>& mapper) {
	string dataset = lc_cfg::get().dataset_name;
	INFO("Reading graph");
	debruijn_graph::ScanBasicGraph(fileName, g, intIds);
	debruijn_graph::LoadKmerMapper(fileName, mapper);
	INFO("Graph read");
}

template<size_t k>
void AddEtalonInfo(const Graph& g, EdgeIndex<k+1, Graph>& index, KmerMapper<k+1, Graph>& mapper, const Sequence& genome, PairedInfoIndices& pairedInfos) {
	for (auto el = lc_cfg::get().ds.etalon_paired_lib.begin(); el != lc_cfg::get().ds.etalon_paired_lib.end(); ++el) {
		INFO("Generating info with read size " << el->read_size << ", insert size " << el->insert_size);

		pairedInfos.push_back(PairedInfoIndexLibrary(g, el->read_size, el->insert_size, el->is_delta, el->de_delta, params.ps.es.etalon_distance_dev, new PairedInfoIndex<Graph>(g, 0)));
		FillEtalonPairedIndex<k> (*pairedInfos.back().pairedInfoIndex, g, index, mapper, el->insert_size, el->read_size, el->de_delta, genome);
	}
}

template<size_t k>
void AddRealInfo(Graph& g, EdgeIndex<k+1, Graph>& index, IdTrackHandler<Graph>& conj_IntIds, PairedInfoIndices& pairedInfos, KmerMapper<k+1, Graph>& mapper,
		bool useNewMetrics) {
	//PairedInfoSimpleSymmertrizer sym(g);

	for (auto rl = lc_cfg::get().ds.paired_lib.begin(); rl != lc_cfg::get().ds.paired_lib.end(); ++rl) {
		size_t insertSize = rl->insert_size;
		size_t readSize = rl->read_size;
		size_t delta = rl->is_delta;
		size_t de_delta = rl->de_delta;
		size_t var = rl->var;
		string dataset = lc_cfg::get().dataset_name;
		pairedInfos.push_back(PairedInfoIndexLibrary(g, readSize, insertSize, delta, de_delta, var, new PairedInfoIndex<Graph>(g)));

		INFO("Reading additional info with read size " << readSize << ", insert size " << insertSize);

		if (rl->precounted && !lc_cfg::get().paired_info_only) {
			//Reading saved paired info
		    typename ScannerTraits<Graph>::Scanner scanner(g, conj_IntIds);
			ScanPairedIndex(rl->precounted_path, scanner, *pairedInfos.back().pairedInfoIndex);
			CheckPairedInfo(g, *pairedInfos.back().pairedInfoIndex, insertSize - 2 * readSize + K);

			pairedInfos.back().raw = new PairedInfoIndex<Graph>(g, 0);
			if (!lc_cfg::get().paired_info_only) {
				//dataScanner.loadPaired(rl->ds.raw, *pairedInfos.back().raw);
		
				pairedInfos.back().has_advanced = rl->has_advanced;
				if (rl->has_advanced) {
					pairedInfos.back().advanced = new PairedInfoIndexLibrary(g, readSize, insertSize, delta, de_delta, var, new PairedInfoIndex<Graph>(g, 0));
					ScanPairedIndex(rl->precounted_path, scanner, *pairedInfos.back().advanced->pairedInfoIndex);
				}
				else {
					pairedInfos.back().advanced = 0;
				}
			}
		}
		else {
			string reads_filename_1 = rl->first;
			string reads_filename_2 = rl->second;

			checkFileExistenceFATAL(reads_filename_1);
			checkFileExistenceFATAL(reads_filename_2);

			io::PairedEasyReader stream(
					std::make_pair(reads_filename_1, reads_filename_2),
					rl->insert_size);


			if (useNewMetrics) {
				FillPairedIndexWithReadCountMetric<k>(g, conj_IntIds, index, mapper, *pairedInfos.back().pairedInfoIndex, stream);
			} else {
				FillPairedIndexWithProductMetric<k>(g, conj_IntIds, index, mapper, *pairedInfos.back().pairedInfoIndex, stream);
			}
		}
		INFO("Done");

//		if (lc_cfg::get().syminfo) {
//			sym.MakeSymmetricInfo(pairedInfos.back());
//		}
	}
}

//void SavePairedInfo(const Graph& g, IdTrackHandler<Graph>& intIds, PairedInfoIndices& pairedInfos, const std::string& fileNamePrefix,
//		bool advEstimator = false) {
//
//	INFO("Saving paired info");
//
//	PrinterTraits<Graph>::Printer printer(g, g.begin(), g.end(), intIds);
//
//	for (auto lib = pairedInfos.begin(); lib != pairedInfos.end(); ++lib) {
//		std::string fileName = fileNamePrefix + "IS" + ToString(lib->insertSize) + "_RS" + ToString(lib->readSize);
//		INFO("Saving to " << fileName);
//
//		if (lc_cfg::get().cluster_paired_info) {
//			GraphDistanceFinder<Graph> dist_finder(g, lib->insertSize, lib->readSize, lib->deDelta);
//
//			if (!advEstimator) {
//				PairedInfoIndex<Graph> clustered_index(g);
//
//				DistanceEstimator<Graph> estimator(g, *(lib->pairedInfoIndex), dist_finder, cfg::get().de.linkage_distance, cfg::get().de.max_distance);
//				estimator.Estimate(clustered_index);
//
//				PrintPairedIndex(fileName + "_clustered", printer, clustered_index);
//			} else {
//
//				PairedInfoIndex<Graph> clustered_index_(g);
//				AdvancedDistanceEstimator<Graph> estimator(g, *(lib->pairedInfoIndex), dist_finder,
//						cfg::get().de.linkage_distance,
//						cfg::get().ade.threshold,
//						cfg::get().ade.range_coeff, cfg::get().ade.delta_coeff,
//						cfg::get().ade.cutoff, cfg::get().ade.minpeakpoints,
//						cfg::get().ade.inv_density, cfg::get().ade.percentage,
//						cfg::get().ade.derivative_threshold);
//
//				estimator.Estimate(clustered_index_);
//
//				PrintPairedIndex(fileName + "acl", printer, clustered_index_);
//			}
//		}
//		if (params.write_raw_paired_info) {
//			PrintPairedIndex(fileName, printer, *(lib->pairedInfoIndex));
//		}
//	}
//
//	INFO("Saved");
//}

void DeleteAdditionalInfo(PairedInfoIndices& pairedInfos) {
	while (pairedInfos.size() > 1) {
		delete pairedInfos.back().pairedInfoIndex;
		pairedInfos.pop_back();
	}
}

void PrintPairedInfo(const Graph& g, PairedInfoIndex<Graph>& pairedInfo, IdTrackHandler<Graph>& intIds) {
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
Sequence PathToSequence(const Graph& g, BidirectionalPath& path) {
	SequenceBuilder result;

	if (!path.empty()) {
		result.append(g.EdgeNucls(path[0]).Subseq(0, K));
	}
	for (auto edge = path.begin(); edge != path.end(); ++edge) {
		result.append(g.EdgeNucls(*edge).Subseq(K));
	}

	return result.BuildSequence();
}

double PathCoverage(const Graph& g, BidirectionalPath& path) {
	double cov = 0.0;
	double len = 0.0;

	for (int i = 0; i < (int) path.size(); ++i) {
		cov += g.coverage(path[i]) * g.length(path[i]);
		len += g.length(path[i]);
	}
	return cov / len;
}

//Output
void OutputPathsAsContigs(const Graph& g, std::vector<BidirectionalPath> paths, const string& filename) {
	INFO("Writing contigs to " << filename);
	osequencestream oss(filename);
	for (auto path = paths.begin(); path != paths.end(); ++path ) {
		oss << PathToSequence(g, *path);
	}
	INFO("Contigs written");
}


//Output only one half of edges
void OutputContigsNoComplement(const Graph& g, const std::string& filename) {
	std::set<EdgeId, Graph::Comparator> filtered(g.ReliableComparatorInstance());
	FilterComlementEdges(g, filtered);

	INFO("Outputting contigs to " << filename);
	osequencestream_with_data_for_scaffold oss(filename);
	for (auto it = filtered.begin(); it != filtered.end(); ++it) {

		oss.setCoverage(g.coverage(*it));
		oss << g.EdgeNucls(*it);
	}
	INFO("Contigs written");
}



void OutputPathsAsContigsNoComplement(const Graph& g, std::vector<BidirectionalPath>& paths,
		const string& filename,
		std::set<int> notToPrint) {

	INFO("Writing contigs to " << filename);
	osequencestream_with_data_for_scaffold oss(filename);

	for (int i = 0; i < (int) paths.size(); i += 2) {
		if (notToPrint.count(i) || notToPrint.count(i + 1) || paths[i].size() == 0) {
			continue;
		}

		oss.setID(paths[i].uid);
		oss.setCoverage(PathCoverage(g, paths[i]));
		oss << PathToSequence(g, paths[i]);
	}

	INFO("Contigs written");
}


} // namespace long_contigs

#endif /* LC_IO_HPP_ */
