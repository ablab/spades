/*
 * tip_clipper.hpp
 *
 *  Created on: Mar 25, 2011
 *      Author: sergey
 */

#ifndef TIP_CLIPPER_HPP_
#define TIP_CLIPPER_HPP_

#include <set>
//#include "edge_graph.hpp"
//#include "utils.hpp"
#include "omni_utils.hpp"
#include "xmath.h"
//
//#define DEFAULT_COVERAGE_BOUND 1000
//#define DEFAULT_RELATIVE_COVERAGE_BOUND 2.0
//#define DEFAULT_MAX_TIP_LENGTH 50

namespace omnigraph {

/**
 * This class removes tips from given graph with the following algorithm: it iterates through all edges of
 * the graph(in order defined by certain comparator) and for each edge checks if this edge is likely to be
 * a tip and if edge is judged to be one it is removed.
 */
template<class Graph, typename Comparator>
class AbstractTipClipper : private boost::noncopyable {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph &graph_;
	Comparator comparator_;
	const size_t max_tip_length_;

	boost::function<void(EdgeId)> removal_handler_;
	//boost::function<double (EdgeId)> qual_handler_;
protected:

	/**
	 * Create TipClipper with specified parameters. Those parameters could probably be replaced later with
	 * certain generic checker class.
	 */
	AbstractTipClipper(Graph &graph, Comparator comparator,
			size_t max_tip_length,
			boost::function<void(EdgeId)> removal_handler = 0) :
			graph_(graph), comparator_(comparator), max_tip_length_(
					max_tip_length), removal_handler_(removal_handler) {

	}
    //
	//TipClipper(Graph &graph, Comparator comparator, size_t max_tip_length,
			//size_t max_coverage, double max_relative_coverage, boost::function<void (EdgeId)> removal_handler = 0, boost::function<double (EdgeId)> qual_handler = 0) :
				//graph_(graph), comparator_(comparator),
				//max_tip_length_(max_tip_length), max_coverage_(max_coverage),
				//max_relative_coverage_(max_relative_coverage), removal_handler_(removal_handler), qual_handler_(qual_handler)  {
	//}

	const Graph& graph() const {
		return graph_;
	}

	/**
	 * This method checks if given vertex topologically looks like end of tip
	 * @param v vertex to be checked
	 * @return true if vertex judged to be tip and false otherwise.
	 */
	bool IsTip(VertexId v) const {
		return graph_.IncomingEdgeCount(v) + graph_.OutgoingEdgeCount(v) == 1;
	}

/**
	 * This method checks if given edge topologically looks like a tip.
	 * @param edge edge vertex to be checked
	 * @return true if edge judged to be tip and false otherwise.
	 */
	bool IsTip(EdgeId edge) const {
		return graph_.length(edge) <= max_tip_length_
				&& (IsTip(graph_.EdgeEnd(edge)) || IsTip(graph_.EdgeStart(edge)))
				&& (graph_.OutgoingEdgeCount(graph_.EdgeStart(edge))
						+ graph_.IncomingEdgeCount(graph_.EdgeEnd(edge)) > 2);
	}

	void CompressSplitVertex(VertexId splitVertex) {
		if (graph_.CanCompressVertex(splitVertex)) {
			EdgeId edge1 = graph_.GetUniqueOutgoingEdge(splitVertex);
			EdgeId edge2 = graph_.GetUniqueIncomingEdge(splitVertex);
			if (IsTip(edge1) || IsTip(edge2)) {
				graph_.CompressVertex(splitVertex);
			}
		}
	}

	bool DeleteTipVertex(VertexId vertex) {
		if (graph_.IsDeadEnd(vertex) && graph_.IsDeadStart(vertex)) {
			graph_.DeleteVertex(vertex);
			return true;
		}
		return false;
	}

	void ProcessVertex(VertexId v) {
		if (!DeleteTipVertex(v)) {
			CompressSplitVertex(v);
		}
	}

	void RemoveTip(EdgeId tip) {
		VertexId start = graph_.EdgeStart(tip);
		VertexId end = graph_.EdgeEnd(tip);
		if (removal_handler_) {
			removal_handler_(tip);
		}
		graph_.DeleteEdge(tip);
		ProcessVertex(start);
		ProcessVertex(end);
	}

	/**
	 * This method checks if given edge is a tip and thus should be removed
	 * @param tip edge to check
	 */
	virtual bool AdditionalCondition(EdgeId e) const = 0;

public:
	/**
	 * Method clips tips of the graph.
	 */
	vector<EdgeId> ClipTips(boost::function<double(EdgeId)> get_total_weight = 0) {
        vector<EdgeId> ans;
        size_t removed = 0;
        size_t removed_with_check = 0;
		TRACE("Tip clipping started");
        TipChecker<Graph> tipchecker(graph_, cfg::get().simp.tc.max_iterations, cfg::get().simp.tc.max_levenshtein);
		for (auto iterator = graph_.SmartEdgeBegin(comparator_); !iterator.IsEnd(); ) {
			EdgeId tip = *iterator;
			TRACE("Checking edge for being tip " << tip);
			if (IsTip(tip)) {
				TRACE(
						"Edge " << tip << " judged to look like tip topologically");
				if (AdditionalCondition(tip)) {
                    TRACE("Additional sequence comparing");
                    removed++;
                    //if (get_total_weight && math::gr(get_total_weight(tip), 0.)){ 
                        //if (qual_handler_ && math::gr(qual_handler_(tip), 0.)){
                            //INFO("Pair INFO FOR GOOD TIP " << graph_.int_id(tip) << " " << get_total_weight(tip) << " " << qual_handler_(tip));
                        //}else{ 
                            //INFO("Pair INFO FOR BAD TIP " << graph_.int_id(tip) << " " << get_total_weight(tip));
                        //}
                    //}

                    if (!cfg::get().simp.tc.advanced_checks || tipchecker.TipCanBeProjected(tip)){
					    TRACE("Edge " << tip << " judged to be tip");
    					removed_with_check++;
                        //if (qual_handler_ && math::eq(qual_handler_(tip), 0.))
                            RemoveTip(tip);
					    TRACE("Edge " << tip << " removed as tip");
                    }else 
                        ans.push_back(tip);
				} else {
					TRACE("Edge " << tip << " judged NOT to be tip");
				}
			} else {
				TRACE("Edge " << tip << " judged NOT to look like tip topologically");
			}
			TRACE("Try to find next edge");
			++iterator;
			TRACE("Use next edge");
		}
		TRACE("Tip clipping finished");
        INFO("REMOVED STATS " << removed_with_check << " " << removed);
        tipchecker.PrintTimeStats();
		Compressor<Graph> compressor(graph_);
		compressor.CompressAllVertices();
        return ans;
	}

// -------------------------------------Clipping tips for Resolver-------------------------------------


	vector<EdgeId> ClipTipsForResolver(boost::function<void(EdgeId)> plotter = 0) {
        vector<EdgeId> ans;
        size_t removed = 0;
        size_t removed_with_check = 0;
		TRACE("Tip clipping started");
        TipChecker<Graph> tipchecker(graph_, cfg::get().simp.tc.max_iterations, cfg::get().simp.tc.max_levenshtein);
		for (auto iterator = graph_.SmartEdgeBegin(comparator_); !iterator.IsEnd(); ) {
			EdgeId tip = *iterator;
			TRACE("Checking edge for being tip " << tip);
			if (IsTip(tip)) {
				TRACE("Edge " << tip << " judged to look like tip topologically");
                if (AdditionalCondition(tip)){
                    if (plotter) 
                        plotter(tip);
                    
                    TRACE("Additional sequence comparing");
                    removed++;
                    //if (get_total_weight && math::gr(get_total_weight(tip), 0.)){ 
                        //if (qual_handler_ && math::gr(qual_handler_(tip), 0.)){
                            //INFO("Pair INFO FOR GOOD TIP " << graph_.int_id(tip) << " " << get_total_weight(tip) << " " << qual_handler_(tip));
                        //}else{ 
                            //INFO("Pair INFO FOR BAD TIP " << graph_.int_id(tip) << " " << get_total_weight(tip));
                        //}
                    //}
                    //if (qual_handler_ && math::gr(qual_handler_(tip), 0.)) 
                        //INFO("Good edge " << graph_.int_id(tip));

                    if (!cfg::get().simp.tc.advanced_checks || tipchecker.TipCanBeProjected(tip)){
					    TRACE("Edge " << tip << " judged to be tip");
    					removed_with_check++;
                        RemoveTip(tip);
					    TRACE("Edge " << tip << " removed as tip");
                    }else 
                        ans.push_back(tip);
				} else {
					TRACE("Edge " << tip << " judged NOT to be tip");
				}
			} else {
				TRACE(
						"Edge " << tip << " judged NOT to look like tip topologically");
			}TRACE("Try to find next edge");
			++iterator;
			TRACE("Use next edge");
		}
		TRACE("Tip clipping finished");
        INFO("REMOVED STATS " << removed_with_check << " " << removed);
        tipchecker.PrintTimeStats();
		Compressor<Graph> compressor(graph_);
		compressor.CompressAllVertices();
        return ans;
	}
};

template<class Graph, typename Comparator>
class TipClipper: public AbstractTipClipper<Graph, Comparator> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef AbstractTipClipper<Graph, Comparator> base;

	const size_t max_coverage_;
	const double max_relative_coverage_;

	//	void FindTips() {
	//		for (Graph::VertexIterator it = graph_.begin(); it
	//				!= graph_.begin(); ++it) {
	//			if (isTip(*it)) {
	//				tipQueue_.offer(graph_.GetUniqueIncomingEdge(*it));
	//			}
	//		}
	//	}

	double MaxCompetitorCoverage(EdgeId tip, vector<EdgeId> competitors) const {
		double result = 0;
		for (auto it = competitors.begin(); it != competitors.end(); ++it) {
			if (*it != tip)
				result = max(result, this->graph().coverage(*it));
		}
		return result;
	}

	//todo strange semantics, discuss with Anton
	double MaxCompetitorCoverage(EdgeId tip) const {
		return max(
				MaxCompetitorCoverage(
						tip,
						this->graph().OutgoingEdges(
								this->graph().EdgeStart(tip))),
				MaxCompetitorCoverage(
						tip,
						this->graph().IncomingEdges(
								this->graph().EdgeEnd(tip))));
	}

	/*virtual*/
	bool AdditionalCondition(EdgeId tip) const {
		if (this->graph().coverage(tip) > max_coverage_)
			return false;
		double max_coverage = MaxCompetitorCoverage(tip);
		return math::le(this->graph().coverage(tip),
				max_relative_coverage_ * max_coverage);
	}

public:

	TipClipper(Graph &graph, Comparator comparator, size_t max_tip_length,
			size_t max_coverage, double max_relative_coverage,
			boost::function<void(EdgeId)> removal_handler = 0) :
			base(graph, comparator, max_tip_length, removal_handler), max_coverage_(
					max_coverage), max_relative_coverage_(max_relative_coverage) {
	}



private:
	DECL_LOGGER("TipClipper")
};

template<class Graph, typename Comparator>
class TopologyTipClipper : public AbstractTipClipper<Graph, Comparator> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef AbstractTipClipper<Graph, Comparator> base;

	size_t uniqueness_length_;
	UniquePathFinder<Graph> unique_path_finder_;

	boost::optional<EdgeId> PathStart(EdgeId tip, bool outgoing_tip) const {
		vector<EdgeId> edges;
		if (outgoing_tip) {
			edges = this->graph().IncomingEdges(this->graph().EdgeStart(tip));
		} else {
			edges = this->graph().OutgoingEdges(this->graph().EdgeEnd(tip));
		}
		if (edges.size() == 1) {
			return boost::optional<EdgeId>(*edges.begin());
		} else {
			return boost::none;
		}
	}
protected:

	/*virtual*/
	bool AdditionalCondition(EdgeId tip) const {
		vector<EdgeId> unique_path;
		if (this->graph().IsDeadEnd(this->graph().EdgeEnd(tip))
				&& this->graph().CheckUniqueIncomingEdge(this->graph().EdgeStart(tip))) {
			unique_path = unique_path_finder_.UniquePathBackward(
					this->graph().GetUniqueIncomingEdge(this->graph().EdgeStart(tip)));
		} else if (this->graph().IsDeadStart(this->graph().EdgeStart(tip))
				&& this->graph().CheckUniqueOutgoingEdge(this->graph().EdgeEnd(tip))) {
			unique_path = unique_path_finder_.UniquePathForward(
					this->graph().GetUniqueOutgoingEdge(this->graph().EdgeEnd(tip)));
		}
		return CummulativeLength(this->graph(), unique_path) >= uniqueness_length_;
	}

public:
	TopologyTipClipper(Graph &graph, Comparator comparator,
			size_t max_tip_length, size_t uniqueness_length,
			boost::function<void(EdgeId)> removal_handler = 0) :
			base(graph, comparator, max_tip_length, removal_handler), uniqueness_length_(
					uniqueness_length), unique_path_finder_(graph) {
	}

private:
	DECL_LOGGER("TopologyTipClipper")
};
}

#endif /* TIP_CLIPPER_HPP_ */
