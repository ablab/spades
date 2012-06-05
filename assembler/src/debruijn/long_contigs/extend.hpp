//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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
#include "omni/omni_utils.hpp"

namespace long_contigs {

using namespace debruijn_graph;

template<class Graph>
class JumpingHero {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	const Graph& g_;
	const PairedInfoIndex<Graph>& jump_index_;
	//always contains at most 2 edges
	size_t length_bound_;
	size_t invalidation_length_;
	bool forward_direct_;
	std::vector<EdgeId> edges_;
	size_t length_from_long_;
	double weight_threshold_;
	set<VertexId> vertices_to_reach_dest_;

	void ClearHistory() {
		edges_.clear();
		length_from_long_ = 0;
		vertices_to_reach_dest_.clear();
	}

	void FillVerticesToReachDest() {
		static const size_t max_vertex_bound = 2000;
		VERIFY(edges_.size() == 2);
		if (forward_direct_) {
			omnigraph::BackwardReliableBoundedDijkstra<Graph> dijkstra(g_,
					invalidation_length_, max_vertex_bound);
			dijkstra.run(g_.EdgeStart(edges_[1]));
			vertices_to_reach_dest_.clear();
			auto pv = dijkstra.ProcessedVertices();
			vertices_to_reach_dest_.insert(pv.begin(), pv.end()) ;
		} else {
			omnigraph::ReliableBoundedDijkstra<Graph> dijkstra(g_,
					invalidation_length_, max_vertex_bound);
			dijkstra.run(g_.EdgeEnd(edges_[1]));
			auto pv = dijkstra.ProcessedVertices();
			vertices_to_reach_dest_.clear();
			vertices_to_reach_dest_.insert(pv.begin(), pv.end()) ;
		}

		//todo make it more natural
		if (vertices_to_reach_dest_.size() > 1000) {
			WARN("Very many vertices reached");
		}
	}

	bool ConnectedForwardWithInfo(EdgeId e1, EdgeId e2) const {
		TRACE(
				"Checking if edges " << g_.str(e1) << " and " << g_.str(e2) << " are connected forward with pair info");
		vector<PairInfo<EdgeId>> infos = jump_index_.GetEdgePairInfo(e1, e2);
		double s = 0.;
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			if (math::gr(it->d, 0.)) {
				s += it->weight;
			}
		}TRACE(
				"Weight sum is " << s << " while weight_threshold is " << weight_threshold_);
		bool answer = math::gr(s, weight_threshold_);
		TRACE("The answer is " << (answer ? "YES" : "NO"));
		return answer;
	}

	bool ConnectedWithInfo(EdgeId e1, EdgeId e2, bool forward_direct) const {
		TRACE(
				"Checking if edges " << g_.str(e1) << " and " << g_.str(e2) << " are connected with pair info " << (forward_direct ? "forward" : "backward"));
		return forward_direct ?
				ConnectedForwardWithInfo(e1, e2) :
				ConnectedForwardWithInfo(e2, e1);
	}

	boost::optional<EdgeId> UniqueLongProlongation(EdgeId e) const {
		TRACE("Trying to find long prolongation for edge " << g_.str(e));
		boost::optional<EdgeId> answer;
		vector<PairInfo<EdgeId>> infos = jump_index_.GetEdgeInfo(e);
		set<EdgeId> paired_edges;
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			paired_edges.insert(it->second);
		}TRACE(
				"Trying to find long prolongation for edge " << g_.str(e) << " among edges " << g_.str(paired_edges));

		for (auto it = paired_edges.begin(); it != paired_edges.end(); ++it) {
			if (g_.length(*it) > length_bound_ && (*it) != e
					&& ConnectedWithInfo(e, *it, forward_direct_)) {
				if (!answer) {
					answer.reset(*it);
				} else {
					WARN(
							"Several long prolongations were found for edge " << g_.str(e));
					return boost::none;
				}
			}
		}
		if (answer) {
			TRACE(
					"Unique " << (forward_direct_ ? "forward" : "backward") << " prolongation for edge " << g_.str(e) << " found. It is edge ");// << g_.str(*answer));
		} else {
			TRACE(
					"Unique " << (forward_direct_ ? "forward" : "backward") << " prolongation for edge " << g_.str(e) << " wasn't found.");
		}
		return answer;
	}

	void RefreshLongEdges(EdgeId e) {
		VERIFY(g_.length(e) >= length_bound_);

		TRACE("Refreshing information for long edge " << g_.str(e));
		VERIFY(edges_.size() <= 2);
		TRACE("Clearing history");
		ClearHistory();

		TRACE("Adding first edge " << g_.str(e) << " to buffer");
		edges_.push_back(e);

		boost::optional<EdgeId> prolong = UniqueLongProlongation(e);
		if (prolong) {
			TRACE("Adding second edge " << g_.str(e) << " to buffer");
			edges_.push_back(*prolong);
			FillVerticesToReachDest();
		} else {
			TRACE("Couldn't add second edge");
		}
	}

	template<class It>
	void Init(It begin, It end) {
		TRACE("Initializing jumping hero with path: " << g_.str(begin, end));
		for (auto it = begin; it != end; ++it) {
			ProcessEdge(*it);
		}TRACE("Initialization finished");
	}

	template<class Container>
	void Init(const Container& path) {
		if (forward_direct_)
			Init(path.begin(), path.end());
		else
			Init(path.rbegin(), path.rend());
	}

public:
//	JumpingHero(const Graph& g, const PairedInfoIndex<Graph>& jump_index,
//			size_t length_bound, size_t invalidation_length, bool forward_direct,
//			double weight_threshold):
//		g_(g), length_bound_(length_bound), invalidation_length_(invalidation_length),
//		forward_direct_(forward_direct), length_from_long_(0), weight_threshold_(weight_threshold) {
//
//	}

	template<class Container>
	JumpingHero(const Graph& g, const Container& path,
			const PairedInfoIndex<Graph>& jump_index, size_t length_bound,
			size_t invalidation_length, bool forward_direct,
			double weight_threshold) :
			g_(g), jump_index_(jump_index), length_bound_(length_bound), invalidation_length_(
					invalidation_length), forward_direct_(forward_direct), length_from_long_(
					0), weight_threshold_(weight_threshold) {
		TRACE(
				"Creating jumping hero directed " << (forward_direct ? "forward" : "backward") << ". length_bound = " << length_bound_ << ". invalidation_length = " << invalidation_length << ". weight_threshold = " << weight_threshold_);
		Init(path);
	}

	void ProcessEdge(EdgeId e) {
		TRACE("Processing edge " << g_.str(e));
		if (g_.length(e) > length_bound_) {
			TRACE(
					"Edge is longer than length bound of " << length_bound_ << ". Refreshing buffer of long edges.");
			RefreshLongEdges(e);
		} else {
			TRACE("Edge is shorter than length bound of " << length_bound_);
			length_from_long_ += g_.length(e);
			if (length_from_long_ > invalidation_length_) {
				TRACE(
						"length_from_long = " << length_from_long_ << " exceeded invalidation_length = " << invalidation_length_ << ". Clearing history.");
				if (edges_.size() > 1) {
					WARN(
							"Strange situation happened while going from edge " << g_.str(edges_[0]) << " to edge " << g_.str(edges_[1]));
				}
				ClearHistory();
			}
		}TRACE("length_from_long = " << length_from_long_);
	}

	//returns number of long edges that has pair info to e
	size_t CheckEdge(EdgeId e) const {
		TRACE("Checking edge " << g_.str(e) << " for jumping pair info");
		size_t answer = 0;
		for (size_t i = 0; i < edges_.size(); ++i) {
			if (ConnectedWithInfo(edges_[i], e, forward_direct_ != bool(i))) {
				answer++;
			}
		}TRACE("Check result was " << answer);
		return answer;
	}

	bool CheckDestinationPathExistence(VertexId v) {
		if (edges_.size() < 2)
			return true;
		return vertices_to_reach_dest_.count(v) > 0;
	}

private:
	DECL_LOGGER("JumpingHero")
	;
};

// ====== Weight functions ======
//Weight filter
double WeightFunction(double weight) {
	return math::gr(weight, 0.0) ? 1.0 : 0.0;
}

//Weight from a set of libraries
double EdgeLengthExtentionWeight(const Graph& g, BidirectionalPath& path,
		PathLengths& lengths, EdgeId e, PairedInfoIndices& pairedInfo,
		size_t edgesToExclude, bool forward, bool useWeightFunction = false,
		size_t additionalGapLength = 0) {

	double weight = 0;
	for (auto lib = pairedInfo.begin(); lib != pairedInfo.end(); ++lib) {
//		weight += ExtentionWeight(g, path, lengths, e, *lib, edgesToExclude, forward, useWeightFunction, additionalGapLength);
	}
	return weight;
}

double CorrectWeightByAdvanced(double weight, double advWeight) {
	return math::gr(advWeight, 0.0) ?
			weight * params.ps.es.advanced_coeff : weight;
}

//Calculate weight
double GetWeight(omnigraph::PairedInfoIndex<Graph>::PairInfos pairs,
		PairedInfoIndexLibrary& pairedInfoLibrary, int distance,
		int distanceDev, bool useWeightFunction = false,
		omnigraph::PairedInfoIndex<Graph>::PairInfos* ad_pairs = 0) {
	double weight = 0;

	for (auto iter = pairs.begin(); iter != pairs.end(); ++iter) {

		int pairedDistance = rounded_d(*iter);

		if (params.ps.es.use_delta_first) {
			if (iter->variance != 0) {
				distanceDev = iter->variance;
			}
		}

		//Can be modified according to distance comparison
		if (pairedDistance >= distance - distanceDev
				&& pairedDistance <= distance + distanceDev) {
			double w = iter->weight;

			if (params.ps.es.fix_weight) {
				w = pairedInfoLibrary.NormalizeWeight(*iter);
			}

			weight += w;
		}

		if (ad_pairs) {
			for (auto aiter = ad_pairs->begin(); aiter != ad_pairs->end();
					++aiter) {
				int ad_pairedDistance = rounded_d(*aiter);
				//Can be modified according to distance comparison
				if (ad_pairedDistance >= distance - distanceDev
						&& ad_pairedDistance <= distance + distanceDev) {

					weight = CorrectWeightByAdvanced(weight, aiter->weight);
				}
			}
		}
	}

	return useWeightFunction ? WeightFunction(weight) : weight;
}

//Weight fixing coefficient, sum weight of ideal info
int FixingCoefficient(const Graph& g, const BidirectionalPath& path,
		EdgeId edge, size_t edgesToExclude,
		PairedInfoIndexLibrary& pairedInfoLibrary, bool forward) {
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
	int left = std::max(exclLen - is, -rs - pathLen) + is;

	int delta = right - left + 1 - K + pairedInfoLibrary.is_delta;

	return delta; // > 0 ? delta : -1;
}

//Fixing weight value
double WeightFixing(const Graph& g, const BidirectionalPath& path, EdgeId edge,
		size_t edgesToExclude, PairedInfoIndexLibrary& pairedInfoLibrary,
		double weight, bool forward) {
	int coeff = FixingCoefficient(g, path, edge, edgesToExclude,
			pairedInfoLibrary, forward);
	if (coeff < 0 && weight != 0) {
		DEBUG(
				"Strange fixing!!! Weight: " << weight << ", c = " << coeff << ", edge: " << edge << " = " << g.length(edge));
		PrintPath(g, path);
		return 0;
	}

	return weight / (double) coeff;
}

//Calculate weight for particular path extension from one library
double ExtentionWeight(const Graph& g, BidirectionalPath& path,
		PathLengths& lengths, EdgeId e,
		PairedInfoIndexLibrary& pairedInfoLibrary, size_t edgesToExclude,
		bool forward, bool useWeightFunction = false,
		size_t additionalGapLength = 0) {
	double weight = 0;
	int edgeLength = forward ? 0 : g.length(e);
	size_t start = forward ? 0 : edgesToExclude;
	size_t end = forward ? path.size() - edgesToExclude : path.size();

	static int DISTANCE_DEV =
			/*cfg::get().etalon_info_mode ?
					params.ps.es.etalon_distance_dev : */pairedInfoLibrary.var;

	for (size_t i = start; i < end; ++i) {
		EdgeId edge = path[i];
		omnigraph::PairedInfoIndex<Graph>::PairInfos pairs =
				forward ?
						pairedInfoLibrary.pairedInfoIndex->GetEdgePairInfo(edge,
								e) :
						pairedInfoLibrary.pairedInfoIndex->GetEdgePairInfo(e,
								edge);
		int distance = lengths[i] + edgeLength + additionalGapLength;

		double w = 0;
		if (params.ps.es.use_advanced && pairedInfoLibrary.has_advanced) {
			omnigraph::PairedInfoIndex<Graph>::PairInfos ad_pairs =
					forward ?
							pairedInfoLibrary.advanced->pairedInfoIndex->GetEdgePairInfo(
									edge, e) :
							pairedInfoLibrary.advanced->pairedInfoIndex->GetEdgePairInfo(
									e, edge);
			w = GetWeight(pairs, pairedInfoLibrary, distance, DISTANCE_DEV,
					useWeightFunction,
					params.ps.es.use_advanced ? &ad_pairs : 0);
		} else {
			w = GetWeight(pairs, pairedInfoLibrary, distance, DISTANCE_DEV,
					useWeightFunction, 0);
		}
		weight += w;
	}

	return weight;
}

//Weight from a set of libraries
double ExtentionWeight(const Graph& g, BidirectionalPath& path,
		PathLengths& lengths, EdgeId e, PairedInfoIndices& pairedInfo,
		size_t edgesToExclude, bool forward, bool useWeightFunction = false,
		size_t additionalGapLength = 0) {

	double weight = 0;
	for (auto lib = pairedInfo.begin(); lib != pairedInfo.end(); ++lib) {
		weight += ExtentionWeight(g, path, lengths, e, *lib, edgesToExclude,
				forward, useWeightFunction, additionalGapLength);
	}
	return weight;
}

// ====== Extension functions ======

//Check whether selected extension is good enough
EdgeId ExtensionGoodEnough(EdgeId edge, double weight, double threshold) {
	//Condition of passing threshold is to be done
	return weight > threshold ? edge : EdgeId(0);
}

//Check whether selected extension is good enough
EdgeId ExtensionGoodEnough(EdgeId edge, double weight, double threshold,
		const Graph& g, BidirectionalPath& path, PathStopHandler& handler,
		bool forward) {
	//Condition of passing threshold is to be done
	if (weight > threshold) {
		return edge;
	} else {
		handler.AddStop(&path, WEAK_EXTENSION, forward);
		return EdgeId(0);
	}
}

void FindEdges(const Graph& g, EdgeId edge, int depth,
		std::vector<EdgeId>& result, std::vector<int>& distances,
		bool forward) {
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
			auto edges =
					forward ?
							g.OutgoingEdges(g.EdgeEnd(result[j])) :
							g.IncomingEdges(g.EdgeStart(result[j]));
			int len = g.length(result[j]);

			for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
				result.push_back(*iter);
				depths.push_back(i + 1);

				distances.push_back(distances[j] + len);
			}
			--j;
		}
		++i;
	}DEBUG("== Depth info == ");
	long_contigs::PrintAnyPath(g, result);
	for (int i = 0; i < (int) result.size(); ++i) {
		DEBUG("D = " << distances[i] << ", DEPTH = " << depths[i]);
	}
}

//Select only best extensions using forward weights
double FilterExtentionsDeep(const Graph& g, BidirectionalPath& path,
		std::vector<EdgeId>& edges, PathLengths& lengths,
		PairedInfoIndices& pairedInfo, size_t edgesToExclude, bool forward,
		LoopDetector& detector, int depth = 1) {

	std::multimap<double, EdgeId> weights;
	std::vector<EdgeId> result;
	std::vector<int> distances;
	for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
		double weight = 0;
		FindEdges(g, *iter, depth, result, distances, forward);

		for (int i = 0; i < (int) result.size(); ++i) {
			weight += ExtentionWeight(g, path, lengths, result[i], pairedInfo,
					edgesToExclude, forward, false, distances[i]);
		}

		weights.insert(std::make_pair(weight, *iter));
		detector.temp.weights[*iter] = weight;
	}

	DETAILED_DEBUG(
			"Choosing weights deeper (" << depth << "): " << (forward ? "forward" : "backward"))
	for (auto iter = weights.begin(); iter != weights.end(); ++iter) {
		DETAILED_DEBUG(
				iter->second << " (" << g.length(iter->second) << ") = " << iter->first);
	}

	//Filling maximum edges
	edges.clear();
	auto bestEdge = weights.lower_bound(
			(--weights.end())->first / params.ps.es.priority_coeff);
	for (auto maxEdge = bestEdge; maxEdge != weights.end(); ++maxEdge) {
		edges.push_back(maxEdge->second);
	}

	return bestEdge->first;
}

//Select only best extensions
double FilterExtentions(const Graph& g, BidirectionalPath& path,
		std::vector<EdgeId>& edges, PathLengths& lengths,
		PairedInfoIndices& pairedInfo, size_t edgesToExclude, bool forward,
		LoopDetector& detector, JumpingHero<Graph>& hero,
		bool useWeightFunction = false) {

	static const double magic_constant = 1000.0;

	DETAILED_DEBUG(
			"Filtering " << (forward ? "forward" : "backward") << " extensions based on pair info weight");

	std::multimap<double, EdgeId> weights;

	for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
		double weight = ExtentionWeight(g, path, lengths, *iter, pairedInfo,
				edgesToExclude, forward, useWeightFunction);
		weights.insert(std::make_pair(weight, *iter));
		detector.temp.AddAlternative(*iter, weight);
	}

	for (auto iter = weights.begin(); iter != weights.end(); ++iter) {
		DETAILED_DEBUG(
				"Weight for " << g.str(iter->second) << " = " << iter->first);
	}

	//Filling maximum edges
//	edges.clear();
	std::vector<EdgeId> tmp;
	DETAILED_DEBUG(
			"Choosing edges with max weight with coeff " << params.ps.es.priority_coeff);
	auto bestEdge = weights.lower_bound(
			(--weights.end())->first / params.ps.es.priority_coeff);
	for (auto maxEdge = bestEdge; maxEdge != weights.end(); ++maxEdge) {
		tmp.push_back(maxEdge->second);
	}DETAILED_DEBUG("Edges with max weight: " << g.str(tmp));
	if (bestEdge != weights.end()) {
		DETAILED_DEBUG(
				"Best edge: " << g.str(bestEdge->second) << " with weight " << bestEdge->first);
	}

	//todo discuss logic with Anton and Andrew
	if (true/*tmp.size() > 1*/) {
		DETAILED_DEBUG("Processing with JumpingHero");
		size_t max_long_paired = 0;
		bool valid = false;
		EdgeId best;
		//todo think of logic!!!
		for (auto it = edges/*tmp*/.begin(); it != edges/*tmp*/.end(); ++it) {
			size_t curr = hero.CheckEdge(*it);
			if (curr == max_long_paired) {
				valid = false;
			}
			if (curr > max_long_paired) {
				max_long_paired = curr;
				valid = true;
				best = *it;
			}
		}
		if (valid && max_long_paired > 0) {
			DETAILED_DEBUG("Jumping hero found best extension edge " << best);
			edges.clear();
			edges.push_back(best);
			return magic_constant;
		} else if (valid && max_long_paired == 0) {
			DETAILED_DEBUG("Jumping hero failed to find best extension edge based on jump info, trying path criterion");
			bool is_good_edge = false;
			EdgeId good_edge;
			for (auto it = edges/*tmp*/.begin(); it != edges/*tmp*/.end(); ++it) {
				if (hero.CheckDestinationPathExistence(g.EdgeEnd(*it))) {
					if (!is_good_edge) {
						is_good_edge = true;
						good_edge = *it;
					} else {
						is_good_edge = false;
						break;
					}
				}
			}
			if (is_good_edge) {
				DETAILED_DEBUG("Jumping found best extension edge based on path criterion");
				edges.clear();
				edges.push_back(good_edge);
				return magic_constant;
			} else {
				DETAILED_DEBUG("Jumping hero failed to find best extension edge based path criterion");
			}
		} else {
			DETAILED_DEBUG("Jumping hero failed to find best extension edge");
		}
	}
	edges.clear();
	DETAILED_DEBUG("Returning edges with max weight: " << g.str(tmp));

	edges.insert(edges.end(), tmp.begin(), tmp.end());
//	cout << "Best weight is " << bestEdge->first << endl;

	DETAILED_DEBUG("Best edge's weight: " << bestEdge->first);
	return bestEdge->first;
}

//Choose best matching extension
//Threshold to be discussed
EdgeId ChooseExtension(const Graph& g, BidirectionalPath& path,
		std::vector<EdgeId>& edges, PathLengths& lengths,
		PairedInfoIndices& pairedInfo, double * maxWeight,
		size_t edgesToExclude, bool forward, LoopDetector& detector,
		PathStopHandler& handler, JumpingHero<Graph>& hero) {
//	cout << "here" << endl;

	DETAILED_DEBUG(
			"Choosing " << (forward ? "forward" : "backward") << "  extension among edges " << g.str(edges));

	detector.temp.clear();

	if (edges.size() == 0) {
		DETAILED_DEBUG("No edges to choose from");

		handler.AddStop(&path, NO_EXTENSION, forward);
		//TODO: scafolder mode here
		return EdgeId(0);
	}
//	cout << "here" << endl;
	if (edges.size() == 1) {
		DETAILED_DEBUG(
				"Single possible extension edge " << g.str(*edges.begin()));

		if (params.ps.ss.check_trusted) {
			DETAILED_DEBUG("Checking if trusted");
			double weight = ExtentionWeight(g, path, lengths, edges.back(),
					pairedInfo, 0, forward, false);

			if (ExtensionGoodEnough(edges.back(), weight,
					params.ps.ss.trusted_threshold) == EdgeId(0)) {
				DETAILED_DEBUG("No");
				return EdgeId(0);
			} else {
				DETAILED_DEBUG("Yes");
			}
		}

		detector.temp.AddAlternative(edges.back(), 1);
		return edges.back();
	}
//	cout << "here" << endl;

	EdgeId toReturn(0);
	if (params.rs.research_mode && params.rs.force_to_cycle) {
		WARN("Strange mode launched");
		for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
			if (g.length(*edge) == params.rs.cycle_priority_edge) {
				toReturn = *edge;
			}
		}
	}
	static bool useWeightFunctionFirst = params.ps.es.use_weight_function_first;

	if (useWeightFunctionFirst) {
		WARN("Strange path chosen!!! use_weight_function_first set to true!!!");
		FilterExtentions(g, path, edges, lengths, pairedInfo, edgesToExclude,
				forward, detector, hero, true);

		if (edges.size() == 1) {
			static double weightFunThreshold = params.ps.es.weight_fun_threshold;
			*maxWeight = ExtentionWeight(g, path, lengths, edges.back(),
					pairedInfo, edgesToExclude, forward);

			WARN("Shouldn't be here!!!");
			return toReturn == EdgeId(0) ?
					ExtensionGoodEnough(edges.back(), *maxWeight,
							weightFunThreshold, g, path, handler, forward) :
					toReturn;
		}
	}

	*maxWeight = FilterExtentions(g, path, edges, lengths, pairedInfo,
			edgesToExclude, forward, detector, hero);

	static double weightThreshold = params.ps.es.weight_threshold;
	if (edges.size() == 1) {
		DETAILED_DEBUG(
				"Single extension passed filtering: " << g.str(*edges.begin()) << " with weight " << *maxWeight);

		DETAILED_DEBUG(
				"Checking if good enough with threshold " << weightThreshold);

		return toReturn == EdgeId(0) ?
				ExtensionGoodEnough(edges.back(), *maxWeight, weightThreshold,
						g, path, handler, forward) :
				toReturn;
	} else if (edges.size() > 1) {
		DETAILED_DEBUG(
				"Several extension passed filtering: " << g.str(edges) << " with best weight " << *maxWeight);

		if (ExtensionGoodEnough(edges.back(), *maxWeight, weightThreshold)
				== EdgeId(0)) {
//			cout << "here2" << endl;
			DETAILED_DEBUG(
					"Best weight extension didn't pass threshold " << weightThreshold);
			DETAILED_DEBUG("No good extension");
			handler.AddStop(&path, NO_GOOD_EXTENSION, forward);
		} else {
			DETAILED_DEBUG(
					"Best weight extension passed threshold " << weightThreshold);
//			cout << "here3" << endl;
			DETAILED_DEBUG("Cannot choose extension, no obvious maximum");

			DETAILED_DEBUG(
					"Hopefully doing nothing and will return no extension");
			static int maxDepth = params.ps.es.max_depth;
			for (int depth = 1; depth <= maxDepth; ++depth) {
				WARN("Shouldn't be here!!!");
				DETAILED_DEBUG("Trying to look deeper to " << depth);
				*maxWeight = FilterExtentionsDeep(g, path, edges, lengths,
						pairedInfo, edgesToExclude, forward, detector, depth);

				if (edges.size() == 1) {
					return toReturn == EdgeId(0) ?
							ExtensionGoodEnough(edges.back(), *maxWeight,
									weightThreshold, g, path, handler,
									forward) :
							toReturn;
				}
			}DEBUG("Still no obvious selection, will stop growing");
			handler.AddStop(&path, MANY_GOOD_EXTENSIONS, forward);
		}
	}
	return toReturn;
}

//Count edges to be excluded
size_t EdgesToExcludeForward(const Graph& g, BidirectionalPath& path, int from =
		-1) {
	static bool maxCycles = params.ps.ss.max_cycles;
	static LoopDetector detector(g);
	detector.clear();

	if (path.empty()) {
		return 0;
	}

	VertexId currentVertex =
			(from == -1) ? g.EdgeEnd(path.back()) : g.EdgeEnd(path[from]);
	size_t toExclude = 0;

	while (g.CheckUniqueIncomingEdge(currentVertex)) {
		EdgeId e = g.GetUniqueIncomingEdge(currentVertex);
		currentVertex = g.EdgeStart(e);
		++toExclude;

		detector.temp.clear();
		detector.temp.AddAlternative(e);
		detector.AddNewEdge(e, toExclude);
		if (CheckCycle(path, e, detector, maxCycles)) {
			DEBUG("Cycled trivial path");
			return 0;
		}
	}

	return std::min(toExclude, path.size());
}

//Count edges to be excludeD
size_t EdgesToExcludeBackward(const Graph& g, BidirectionalPath& path,
		int from = -1) {
	static bool maxCycles = params.ps.ss.max_cycles;
	static LoopDetector detector(g);
	detector.clear();

	if (path.empty()) {
		return 0;
	}
	VertexId currentVertex =
			(from == -1) ? g.EdgeStart(path.front()) : g.EdgeStart(path[from]);
	size_t toExclude = 0;

	while (g.CheckUniqueOutgoingEdge(currentVertex)) {
		EdgeId e = g.GetUniqueOutgoingEdge(currentVertex);
		currentVertex = g.EdgeEnd(e);
		++toExclude;

		detector.temp.clear();
		detector.temp.AddAlternative(e);
		detector.AddNewEdge(e, toExclude);
		if (CheckCycle(path, e, detector, maxCycles)) {
			DEBUG("Cycled trivial path");
			return 0;
		}
	}

	return std::min(toExclude, path.size());
}

void RecountDetectorForward(const Graph& g, BidirectionalPath& path,
		PairedInfoIndices& pairedInfo, LoopDetector& detector) {
	BidirectionalPath emulPath;
	PathLengths emulLengths;
	detector.clear();

	DETAILED_DEBUG("Recounting detector forward");

	for (int i = 0; i < (int) path.size(); ++i) {
		DETAILED_DEBUG(i);
		size_t edgesToExclude = EdgesToExcludeForward(g, emulPath);

		detector.temp.clear();
		if (g.OutgoingEdgeCount(g.EdgeStart(path[i])) == 1) {
			detector.temp.AddAlternative(path[i]);
			detector.AddNewEdge(path[i], i);
		} else {
			auto edges = g.OutgoingEdges(g.EdgeStart(path[i]));

			for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
				double weight =
						(i == 0 || (int) edgesToExclude >= i) ?
								1.0 :
								ExtentionWeight(g, emulPath, emulLengths, *iter,
										pairedInfo, edgesToExclude, true,
										false);

				detector.temp.AddAlternative(*iter, weight);
			}
			detector.AddNewEdge(path[i], i, detector.temp.weights[path[i]]);
		}

		emulPath.push_back(path[i]);
		IncreaseLengths(g, emulLengths, path[i], true);
	}
}

void RecountDetectorBackward(const Graph& g, BidirectionalPath& path,
		PairedInfoIndices& pairedInfo, LoopDetector& detector) {
	BidirectionalPath emulPath;
	PathLengths emulLengths;
	detector.clear();

	DETAILED_DEBUG("Recounting detector backward");

	for (int i = path.size() - 1; i >= 0; --i) {
		size_t edgesToExclude = EdgesToExcludeBackward(g, emulPath);

		detector.temp.clear();
		if (g.IncomingEdgeCount(g.EdgeEnd(path[i])) == 1) {
			detector.temp.AddAlternative(path[i]);
			detector.AddNewEdge(path[i], path.size() - 1 - i);
		} else {
			auto edges = g.IncomingEdges(g.EdgeEnd(path[i]));

			for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
				double weight =
						(i == (int) path.size() - 1
								|| (int) edgesToExclude
										>= (int) path.size() - 1 - i) ?
								1.0 :
								ExtentionWeight(g, emulPath, emulLengths, *iter,
										pairedInfo, edgesToExclude, false,
										false);

				detector.temp.AddAlternative(*iter, weight);
			}

			detector.AddNewEdge(path[i], path.size() - 1 - i,
					detector.temp.weights[path[i]]);
		}

		emulPath.push_front(path[i]);
		IncreaseLengths(g, emulLengths, path[i], false);
	}
}

EdgeId FindScafoldExtension(const Graph& g, BidirectionalPath& path,
		PathLengths& lengths, PairedInfoIndices& pairedInfo, double * maxWeight,
		bool forward) {

	return EdgeId(0);
}

} //namespace long_contigs

#endif /* EXTEND_HPP_ */
