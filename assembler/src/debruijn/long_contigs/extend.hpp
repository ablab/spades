/*
 * extend.hpp
 *
 *  Created on: Oct 2, 2011
 *      Author: andrey
 */

#ifndef EXTEND_HPP_
#define EXTEND_HPP_

#include <cmath>

#include "lc_common.hpp"
#include "loop.hpp"


namespace long_contigs {

using namespace debruijn_graph;

// ====== Weight functions ======
//Weight filter
double WeightFunction(double weight) {
	return math::gr(weight, 0.0) ? 1.0 : 0.0;
}

//Weight from a set of libraries
double EdgeLengthExtentionWeight(Graph& g, BidirectionalPath& path, PathLengths& lengths, EdgeId e, PairedInfoIndices& pairedInfo,
		size_t edgesToExclude, bool forward, bool useWeightFunction = false, size_t additionalGapLength = 0) {

	double weight = 0;
	for (auto lib = pairedInfo.begin(); lib != pairedInfo.end(); ++lib) {
//		weight += ExtentionWeight(g, path, lengths, e, *lib, edgesToExclude, forward, useWeightFunction, additionalGapLength);
	}
	return weight;
}


double CorrectWeightByAdvanced(double weight, double advWeight) {
	return math::gr(advWeight, 0.0) ? weight * lc_cfg::get().ps.es.advanced_coeff
			: weight;
}



//Calculate weight
double GetWeight(omnigraph::PairedInfoIndex<Graph>::PairInfos pairs, PairedInfoIndexLibrary& pairedInfoLibrary,
		int distance, int distanceDev, bool useWeightFunction = false,
		omnigraph::PairedInfoIndex<Graph>::PairInfos* ad_pairs = 0) {
	double weight = 0;

	for (auto iter = pairs.begin(); iter != pairs.end(); ++iter) {


		int pairedDistance = rounded_d(*iter);

		if (lc_cfg::get().ps.es.use_delta_first) {
			if (iter->variance != 0) {
				distanceDev = iter->variance;
			}
		}

		//Can be modified according to distance comparison
		if (pairedDistance >= distance - distanceDev &&
				pairedDistance <= distance + distanceDev) {
			double w = iter->weight;

			if (lc_cfg::get().ps.es.fix_weight) {
				w = pairedInfoLibrary.NormalizeWeight(*iter);
			}

			weight += w;
		}

		if (ad_pairs) {
			for (auto aiter = ad_pairs->begin(); aiter != ad_pairs->end(); ++aiter) {
				int ad_pairedDistance = rounded_d(*aiter);
				//Can be modified according to distance comparison
				if (ad_pairedDistance >= distance - distanceDev &&
						ad_pairedDistance <= distance + distanceDev) {

					weight = CorrectWeightByAdvanced(weight, aiter->weight);
				}
			}
		}
	}




	return useWeightFunction ? WeightFunction(weight) : weight;
}

//Weight fixing coefficient, sum weight of ideal info
int FixingCoefficient(Graph& g, const BidirectionalPath& path, EdgeId edge, size_t edgesToExclude, PairedInfoIndexLibrary& pairedInfoLibrary, bool forward) {
	int pathLen = 0;
	size_t start = forward ? 0 : edgesToExclude;
	size_t end = forward ? path.size() - edgesToExclude : path.size();
	for (size_t i = start; i < end; ++i) {
		pathLen += g.length(path[i]);
	}
	int exclLen = PathLength(g, path) - pathLen;
	int edgeLen = g.length(edge);

	int is = pairedInfoLibrary.insertSize;
	int rs = pairedInfoLibrary.readSize;

	int right = std::min(is, exclLen + edgeLen + rs);
	int left = std::max(exclLen - is, - rs - pathLen) + is;

	int delta = right - left + 1 - K + pairedInfoLibrary.is_delta;

	return delta; // > 0 ? delta : -1;
}

//Fixing weight value
double WeightFixing(Graph& g, const BidirectionalPath& path, EdgeId edge, size_t edgesToExclude, PairedInfoIndexLibrary& pairedInfoLibrary, double weight, bool forward) {
	int coeff = FixingCoefficient(g, path, edge, edgesToExclude, pairedInfoLibrary, forward);
	if (coeff < 0 && weight != 0) {
		INFO("Strange fixing!!! Weight: " << weight << ", c = " << coeff << ", edge: " << edge << " = " << g.length(edge));
		PrintPath(g ,path);
		return 0;
	}

	return weight/(double) coeff;
}

//Calculate weight for particular path extension from one library
double ExtentionWeight(Graph& g, BidirectionalPath& path, PathLengths& lengths, EdgeId e, PairedInfoIndexLibrary& pairedInfoLibrary,
		size_t edgesToExclude, bool forward, bool useWeightFunction = false, size_t additionalGapLength = 0) {
	double weight = 0;
	int edgeLength = forward ? 0 : g.length(e);
	size_t start = forward ? 0 : edgesToExclude;
	size_t end = forward ? path.size() - edgesToExclude : path.size();

	static int DISTANCE_DEV = cfg::get().etalon_info_mode ? lc_cfg::get().ps.es.etalon_distance_dev : pairedInfoLibrary.var;

	for(size_t i = start; i < end; ++i) {
		EdgeId edge = path[i];
		omnigraph::PairedInfoIndex<Graph>::PairInfos pairs =
				forward ? pairedInfoLibrary.pairedInfoIndex->GetEdgePairInfo(edge, e) : pairedInfoLibrary.pairedInfoIndex->GetEdgePairInfo(e, edge);
		int distance = lengths[i] + edgeLength + additionalGapLength;

		double w = 0;
		if (pairedInfoLibrary.has_advanced) {
			omnigraph::PairedInfoIndex<Graph>::PairInfos ad_pairs = forward ? pairedInfoLibrary.advanced->pairedInfoIndex->GetEdgePairInfo(edge, e) : pairedInfoLibrary.advanced->pairedInfoIndex->GetEdgePairInfo(e, edge);
			w = GetWeight(pairs, pairedInfoLibrary, distance, DISTANCE_DEV, useWeightFunction, lc_cfg::get().ps.es.use_advanced ? &ad_pairs : 0);
		} else {
			w = GetWeight(pairs, pairedInfoLibrary, distance, DISTANCE_DEV, useWeightFunction, 0);
		}
		weight += w;
	}

	return weight;
}

//Weight from a set of libraries
double ExtentionWeight(Graph& g, BidirectionalPath& path, PathLengths& lengths, EdgeId e, PairedInfoIndices& pairedInfo,
		size_t edgesToExclude, bool forward, bool useWeightFunction = false, size_t additionalGapLength = 0) {

	double weight = 0;
	for (auto lib = pairedInfo.begin(); lib != pairedInfo.end(); ++lib) {
		weight += ExtentionWeight(g, path, lengths, e, *lib, edgesToExclude, forward, useWeightFunction, additionalGapLength);
	}
	return weight;
}

// ====== Extension functions ======

//Check whether selected extension is good enough
EdgeId ExtensionGoodEnough(EdgeId edge, double weight, double threshold) {
	//Condition of passing threshold is to be done
	return weight > threshold ? edge : 0;
}

//Check whether selected extension is good enough
EdgeId ExtensionGoodEnough(EdgeId edge, double weight, double threshold, Graph& g, BidirectionalPath& path, PathStopHandler& handler, bool forward) {
	//Condition of passing threshold is to be done
	if (weight > threshold) {
		return edge;
	} else {
		handler.AddStop(&path, WEAK_EXTENSION, forward);
		return 0;
	}
}



void FindEdges(Graph& g, EdgeId edge, int depth, std::vector<EdgeId>& result, std::vector<int>& distances, bool forward) {
	std::vector<int> depths;
	result.clear();
	distances.clear();
	int i = 0;

	result.push_back(edge);
	depths.push_back(i);
	distances.push_back(0);

	while (i < depth) {
		int j = result.size() - 1;
		while (j >= 0 && depths[j] == i) {
			auto edges = forward ? g.OutgoingEdges(g.EdgeEnd(result[j])) : g.IncomingEdges(g.EdgeStart(result[j]));
			int len = g.length(result[j]);

			for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
				result.push_back(*iter);
				depths.push_back(i + 1);

				distances.push_back(distances[j] + len);
			}
			--j;
		}
		++i;
	}
	INFO("== Depth info == ");
	PrintPath(g, result);
	for (int i = 0; i < (int) result.size(); ++i) {
		INFO("D = " << distances[i] << ", DEPTH = " << depths[i]);
	}
}

//Select only best extensions using forward weights
double FilterExtentionsDeep(Graph& g, BidirectionalPath& path, std::vector<EdgeId>& edges,
		PathLengths& lengths, PairedInfoIndices& pairedInfo, size_t edgesToExclude, bool forward,
		LoopDetector& detector,
		int depth = 1) {

	std::multimap<double, EdgeId> weights;
	std::vector<EdgeId> result;
	std::vector<int> distances;
	for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
		double weight = 0;
		FindEdges(g, *iter, depth, result, distances, forward);

		for (int i = 0; i < (int) result.size(); ++i) {
			weight += ExtentionWeight(g, path, lengths, result[i], pairedInfo, edgesToExclude, forward, false, distances[i]);
		}

		weights.insert(std::make_pair(weight, *iter));
		detector.temp.weights[*iter] = weight;
	}

	DETAILED_INFO("Choosing weights deeper (" << depth << "): " << (forward ? "forward" : "backward"))
	for (auto iter = weights.begin(); iter != weights.end(); ++iter) {
		DETAILED_INFO(iter->second << " (" << g.length(iter->second) << ") = " << iter->first);
	}

	//Filling maximum edges
	edges.clear();
	auto bestEdge = weights.lower_bound((--weights.end())->first / lc_cfg::get().ps.es.priority_coeff);
	for (auto maxEdge = bestEdge; maxEdge != weights.end(); ++maxEdge) {
		edges.push_back(maxEdge->second);
	}

	return bestEdge->first;
}

//Select only best extensions
double FilterExtentions(Graph& g, BidirectionalPath& path, std::vector<EdgeId>& edges,
		PathLengths& lengths, PairedInfoIndices& pairedInfo, size_t edgesToExclude, bool forward,
		LoopDetector& detector,
		bool useWeightFunction = false) {

	std::multimap<double, EdgeId> weights;

	for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
		double weight = ExtentionWeight(g, path, lengths, *iter, pairedInfo, edgesToExclude, forward, useWeightFunction);
		weights.insert(std::make_pair(weight, *iter));
		detector.temp.AddAlternative(*iter, weight);
	}

	DETAILED_INFO("Choosing weights " << (forward ? "forward" : "backward"))
	for (auto iter = weights.begin(); iter != weights.end(); ++iter) {
		DETAILED_INFO(iter->second << " (" << g.length(iter->second) << ") = " << iter->first);
	}

	//Filling maximum edges
	edges.clear();
	auto bestEdge = weights.lower_bound((--weights.end())->first / lc_cfg::get().ps.es.priority_coeff);
	for (auto maxEdge = bestEdge; maxEdge != weights.end(); ++maxEdge) {
		edges.push_back(maxEdge->second);
	}

	return bestEdge->first;
}

//Choose best matching extension
//Threshold to be discussed
EdgeId ChooseExtension(Graph& g, BidirectionalPath& path, std::vector<EdgeId>& edges,
		PathLengths& lengths, PairedInfoIndices& pairedInfo, double * maxWeight, size_t edgesToExclude, bool forward,
		LoopDetector& detector,
		PathStopHandler& handler) {

	detector.temp.clear();

	if (edges.size() == 0) {
		handler.AddStop(&path, NO_EXTENSION, forward);
		//TODO: scafolder mode here
		return 0;
	}

	if (edges.size() == 1) {

		if (lc_cfg::get().ps.ss.check_trusted) {
			double weight =
					ExtentionWeight(g, path, lengths, edges.back(), pairedInfo, 0, forward, false);

			if (ExtensionGoodEnough(edges.back(), weight, lc_cfg::get().ps.ss.trusted_threshold) == 0) {
				return 0;
			}
		}

		detector.temp.AddAlternative(edges.back(), 1);
		return edges.back();
	}
	EdgeId toReturn = 0;
	if (lc_cfg::get().rs.research_mode && lc_cfg::get().rs.force_to_cycle) {
		for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
			if (g.length(*edge) == lc_cfg::get().rs.cycle_priority_edge) {
				toReturn = *edge;
			}
		}
	}

	static bool useWeightFunctionFirst = lc_cfg::get().ps.es.use_weight_function_first;

	if (useWeightFunctionFirst) {
		FilterExtentions(g, path, edges, lengths, pairedInfo, edgesToExclude, forward, detector, true);

		if (edges.size() == 1) {
			static double weightFunThreshold = lc_cfg::get().ps.es.weight_fun_threshold;
			*maxWeight = ExtentionWeight(g, path, lengths, edges.back(), pairedInfo, edgesToExclude, forward);

			return toReturn == 0 ? ExtensionGoodEnough(edges.back(), *maxWeight, weightFunThreshold, g, path, handler, forward) : toReturn;
		}
	}

	*maxWeight = FilterExtentions(g, path, edges, lengths, pairedInfo, edgesToExclude, forward, detector);

	static double weightThreshold = lc_cfg::get().ps.es.weight_threshold;
	if (edges.size() == 1) {
		return toReturn == 0 ? ExtensionGoodEnough(edges.back(), *maxWeight, weightThreshold, g, path, handler, forward) : toReturn;
	}
	else if (edges.size() > 1) {
		if (ExtensionGoodEnough(edges.back(), *maxWeight, weightThreshold) == 0) {
			DETAILED_INFO("No good extension");
			handler.AddStop(&path, NO_GOOD_EXTENSION, forward);
		} else {
			DETAILED_INFO("Cannot choose extension, no obvious maximum");

			static int maxDepth = lc_cfg::get().ps.es.max_depth;
			for (int depth = 1; depth <= maxDepth; ++depth) {
				DETAILED_INFO("Trying to look deeper to " << depth);
				*maxWeight = FilterExtentionsDeep(g, path, edges, lengths, pairedInfo, edgesToExclude, forward, detector, depth);

				if (edges.size() == 1) {
					return toReturn == 0 ? ExtensionGoodEnough(edges.back(), *maxWeight, weightThreshold, g, path, handler, forward) : toReturn;
				}
			}
			INFO("Still no obvious selection, will stop growing");
			handler.AddStop(&path, MANY_GOOD_EXTENSIONS, forward);
		}
	}
	return toReturn;
}

//Count edges to be excluded
size_t EdgesToExcludeForward(Graph& g, BidirectionalPath& path, int from = -1) {
	static bool maxCycles = lc_cfg::get().ps.ss.max_cycles;
	static LoopDetector detector;
	detector.clear();

	if (path.empty()) {
		return 0;
	}

	VertexId currentVertex = (from == -1) ? g.EdgeEnd(path.back()) : g.EdgeEnd(path[from]);
	size_t toExclude = 0;

	while (g.CheckUniqueIncomingEdge(currentVertex)) {
		EdgeId e = g.GetUniqueIncomingEdge(currentVertex);
		currentVertex = g.EdgeStart(e);
		++toExclude;

		detector.temp.clear();
		detector.temp.AddAlternative(e);
		detector.AddNewEdge(e, toExclude);
		if (CheckCycle(path, e, detector, maxCycles)) {
			INFO("Cycled trivial path");
			return 0;
		}
	}

	return std::min(toExclude, path.size());
}

//Count edges to be excludeD
size_t EdgesToExcludeBackward(Graph& g, BidirectionalPath& path, int from = -1) {
	static bool maxCycles = lc_cfg::get().ps.ss.max_cycles;
	static LoopDetector detector;
	detector.clear();

	if (path.empty()) {
		return 0;
	}
	VertexId currentVertex = (from == -1) ? g.EdgeStart(path.front()) : g.EdgeStart(path[from]);
	size_t toExclude = 0;

	while (g.CheckUniqueOutgoingEdge(currentVertex)) {
		EdgeId e = g.GetUniqueOutgoingEdge(currentVertex);
		currentVertex = g.EdgeEnd(e);
		++toExclude;

		detector.temp.clear();
		detector.temp.AddAlternative(e);
		detector.AddNewEdge(e, toExclude);
		if (CheckCycle(path, e, detector, maxCycles)) {
			INFO("Cycled trivial path");
			return 0;
		}
	}

	return std::min(toExclude, path.size());
}

void RecountDetectorForward(Graph& g, BidirectionalPath& path, PairedInfoIndices& pairedInfo, LoopDetector& detector) {
	BidirectionalPath emulPath;
	PathLengths emulLengths;
	detector.clear();

	DETAILED_INFO("Recounting detector forward");

	for (int i = 0; i < (int) path.size(); ++i) {
		DETAILED_INFO(i);
		size_t edgesToExclude = EdgesToExcludeForward(g, emulPath);

		detector.temp.clear();
		if (g.OutgoingEdgeCount(g.EdgeStart(path[i])) == 1) {
			detector.temp.AddAlternative(path[i]);
			detector.AddNewEdge(path[i], i);
		}
		else {
			auto edges = g.OutgoingEdges(g.EdgeStart(path[i]));

			for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
				double weight = (i == 0 || (int) edgesToExclude >= i) ?
						1.0 : ExtentionWeight(g, emulPath, emulLengths, *iter, pairedInfo, edgesToExclude, true, false);

				detector.temp.AddAlternative(*iter, weight);
			}
			detector.AddNewEdge(path[i], i, detector.temp.weights[path[i]]);
		}

		emulPath.push_back(path[i]);
		IncreaseLengths(g, emulLengths, path[i], true);
	}
}

void RecountDetectorBackward(Graph& g, BidirectionalPath& path, PairedInfoIndices& pairedInfo, LoopDetector& detector) {
	BidirectionalPath emulPath;
	PathLengths emulLengths;
	detector.clear();

	DETAILED_INFO("Recounting detector backward");

	for (int i = path.size() - 1; i >= 0; --i) {
		size_t edgesToExclude = EdgesToExcludeBackward(g, emulPath);

		detector.temp.clear();
		if (g.IncomingEdgeCount(g.EdgeEnd(path[i])) == 1) {
			detector.temp.AddAlternative(path[i]);
			detector.AddNewEdge(path[i], path.size() - 1 - i);
		}
		else {
			auto edges = g.IncomingEdges(g.EdgeEnd(path[i]));

			for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
				double weight = (i == (int) path.size() - 1 || (int) edgesToExclude >= (int) path.size() - 1 - i) ?
						1.0 : ExtentionWeight(g, emulPath, emulLengths, *iter, pairedInfo, edgesToExclude, false, false);

				detector.temp.AddAlternative(*iter, weight);
			}

			detector.AddNewEdge(path[i], path.size() - 1 - i, detector.temp.weights[path[i]]);
		}

		emulPath.push_front(path[i]);
		IncreaseLengths(g, emulLengths, path[i], false);
	}
}

EdgeId FindScafoldExtension(Graph& g, BidirectionalPath& path,
		PathLengths& lengths, PairedInfoIndices& pairedInfo, double * maxWeight, bool forward) {


	return 0;
}

} //namespace long_contigs


#endif /* EXTEND_HPP_ */
