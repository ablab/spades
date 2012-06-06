//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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
public:
	const size_t max_tip_length_;
private:	
    boost::function<void(EdgeId)> removal_handler_;
	//boost::function<double (EdgeId)> qual_handler_;
protected:


	/**
	 * Create TipClipper with specified parameters. Those parameters could probably be replaced later with
	 * certain generic checker class.
	 */
	AbstractTipClipper(Graph &graph, Comparator comparator,
			size_t max_tip_length,
			boost::function<void(EdgeId)> removal_handler = 0, boost::function<double(EdgeId)> qual_f = 0) :
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

	const Graph& graph() const{
		return graph_;
	}

    Graph& graph(){
        return graph_;   
    }

	const Comparator& comparator() const {
		return comparator_;
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
//	EdgeId edge1 = graph_.GetUniqueOutgoingEdge(splitVertex);
//	EdgeId edge2 = graph_.GetUniqueIncomingEdge(splitVertex);
//	if (IsTip(edge1) || IsTip(edge2)) {
				graph_.CompressVertex(splitVertex);
//	}
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
    virtual void ClipTips() = 0;   


private:
	DECL_LOGGER("AbstractTipClipper")
};

template<class Graph, typename Comparator>
class DefaultTipClipper: public AbstractTipClipper<Graph, Comparator> {

typedef typename Graph::EdgeId EdgeId;
typedef typename Graph::VertexId VertexId;
typedef AbstractTipClipper<Graph, Comparator> base;
        
private:
    const size_t max_coverage_;
    const double max_relative_coverage_;

	double MaxCompetitorCoverage(EdgeId tip, vector<EdgeId> competitors) const {
		double result = 0;
		for (auto it = competitors.begin(); it != competitors.end(); ++it) {
			if (*it != tip)
				result = max(result, this->graph().coverage(*it));
		}
		return result;
	}

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

	bool AdditionalCondition(EdgeId tip) const {
		if (this->graph().coverage(tip) > max_coverage_)
			return false;
		double max_coverage = MaxCompetitorCoverage(tip) + 1;
		return math::le(this->graph().coverage(tip),
				max_relative_coverage_ * max_coverage);
	}

public:

	DefaultTipClipper(Graph &graph, Comparator comparator, size_t max_tip_length, size_t max_coverage, double max_relative_coverage, 
			boost::function<void(EdgeId)> removal_handler = 0, boost::function<double(EdgeId)> qual_f = 0) :
			base(graph, comparator, max_tip_length, removal_handler, qual_f), max_coverage_(
					max_coverage), max_relative_coverage_(max_relative_coverage) {
	}

    void ClipTips(){
		TRACE("Tip clipping (maximal corruption mode) started");
        size_t removed = 0;
        for (auto iterator = this->graph().SmartEdgeBegin(this->comparator()); !iterator.IsEnd(); ) {
			EdgeId tip = *iterator;
			TRACE("Checking edge for being tip "  << this->graph().str(tip));
			if (IsTip(tip)) {
				TRACE("Edge "  << this->graph().str(tip) << " judged to look like a tip topologically");
				if (AdditionalCondition(tip)) {
                    TRACE("Edge "  << this->graph().str(tip) << " judged to be a tip");
                    RemoveTip(tip);
                    removed++;
                    TRACE("Edge "  << tip << " removed as a tip");
				} else {
					TRACE("Edge "  << this->graph().str(tip) << " judged NOT to be tip");
				}
			} else {
				TRACE("Edge "  << this->graph().str(tip) << " judged NOT to look like tip topologically");
			}
			TRACE("Try to find next edge");
			++iterator;
			TRACE("Use next edge");
		}
		TRACE("Tip clipping finished");

        DEBUG("REMOVED " << removed);
		Compressor<Graph> compressor(this->graph());
		compressor.CompressAllVertices();
    }

private:
    DECL_LOGGER("DefaultTipClipper")
};


template<class Graph, typename Comparator>
class AdvancedTipClipper: public AbstractTipClipper<Graph, Comparator> {

private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef AbstractTipClipper<Graph, Comparator> base;
    TipLock<EdgeId> tip_lock;

    const size_t max_coverage_;
    const double max_relative_coverage_;
    const size_t max_iterations_;
    const size_t max_levenshtein_;
    const size_t max_ec_length_;

	double MaxCompetitorCoverage(EdgeId tip, vector<EdgeId> competitors) const {
		double result = 0;
		for (auto it = competitors.begin(); it != competitors.end(); ++it) {
			if (*it != tip)
				result = max(result, this->graph().coverage(*it));
		}
		return result;
	}

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

	double MinCompetitorCoverage(EdgeId tip, vector<EdgeId> competitors) const {
		double result = 1000000; //inf
		for (auto it = competitors.begin(); it != competitors.end(); ++it) {
			if (*it != tip)
				result = min(result, this->graph().coverage(*it)/this->graph().length(*it));
		}
        //if (result == 1000000) WARN("INFINITY REACHED, WHILE SEEKING FOR MINIMUM");
		return result;
	}

	double MinCompetitorCoverage(EdgeId tip) const {
		return min(
				MinCompetitorCoverage(
						tip,
						this->graph().OutgoingEdges(
								this->graph().EdgeStart(tip))),
				MinCompetitorCoverage(
						tip,
						this->graph().IncomingEdges(
								this->graph().EdgeEnd(tip))));
	}

    bool CheckAllAlternativesAreTips(EdgeId tip) const {
        VertexId start = this->graph().EdgeStart(tip);
        TRACE("Check started");
        VertexId end = this->graph().EdgeEnd(tip);
        for (size_t i = 0; i<this->graph().OutgoingEdgeCount(start); ++i){
            EdgeId edge = this->graph().OutgoingEdges(start)[i];
            if (edge != tip){
                if (!IsTip(edge))
                    return false;
            }
        }
        for (size_t i = 0; i<this->graph().IncomingEdgeCount(end); ++i){
            EdgeId edge = this->graph().IncomingEdges(end)[i];
            if (edge != tip){
                if (!IsTip(edge))
                    return false;
            }
        }
        TRACE("Check finished");
        return true;
    }

    //TODO: remove constants
    /// checking whether the next edge after tip is very long, then we'd rather remove it
    bool CheckUniqueExtension(EdgeId tip) const {
        static const size_t mid_edge = 200;
        static const size_t long_edge = 1500;
        bool backward = IsTip(this->graph().EdgeStart(tip));
        if (backward){
            VertexId vertex = this->graph().EdgeEnd(tip);
            for (size_t i = 0; i<this->graph().IncomingEdgeCount(vertex); ++i) 
                if (this->graph().length(this->graph().IncomingEdges(vertex)[i]) < mid_edge) return false;
            if (this->graph().IncomingEdgeCount(vertex) == 2 && this->graph().OutgoingEdgeCount(vertex) == 1)
                return (this->graph().length(this->graph().OutgoingEdges(vertex)[0]) > long_edge);
        }else{
            VertexId vertex = this->graph().EdgeStart(tip);
            for (size_t i = 0; i<this->graph().OutgoingEdgeCount(vertex); ++i) 
                if (this->graph().length(this->graph().OutgoingEdges(vertex)[i]) < mid_edge) return false;
            if (this->graph().OutgoingEdgeCount(vertex) == 2 && this->graph().IncomingEdgeCount(vertex) == 1)
                return (this->graph().length(this->graph().IncomingEdges(vertex)[0]) > long_edge);
        }
        return false;
    }

    bool TipHasVeryLowRelativeCoverage(EdgeId tip) const {
        double max_coverage = MaxCompetitorCoverage(tip);
        return math::ls(200.*this->graph().coverage(tip), max_coverage);
    }

    bool TipHasLowRelativeCoverage(EdgeId tip, bool final_stage = false) const {
        double min_covlen = MinCompetitorCoverage(tip);  
        if (final_stage && this->graph().length(tip) < this->graph().k() / 2) 
            return true;
        return math::ls(this->graph().coverage(tip)/this->graph().length(tip), min_covlen);
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

	AdvancedTipClipper(Graph &graph, Comparator comparator, size_t max_tip_length,
			size_t max_coverage, double max_relative_coverage, size_t max_iterations, size_t max_levenshtein, size_t max_ec_length,
			boost::function<void(EdgeId)> removal_handler = 0) :
			base(graph, comparator, max_tip_length, removal_handler), max_coverage_(
					max_coverage), max_relative_coverage_(max_relative_coverage), max_iterations_(max_iterations),
                            max_levenshtein_(max_levenshtein), max_ec_length_(max_ec_length) {
	}

    // Method deletes tips from the graph carefully, its work depends on the number of simplification iteration

	void ClipTips(bool final_stage){
        size_t removed = 0;
        size_t removed_with_check = 0;
        size_t locked = 0;
		TRACE("Tip clipping started");
        
        TipChecker<Graph> tipchecker(this->graph(), tip_lock, max_iterations_, max_levenshtein_, this->max_tip_length_, max_ec_length_);
         
		for (auto iterator = this->graph().SmartEdgeBegin(this->comparator()); !iterator.IsEnd(); ++iterator) {
			EdgeId tip = *iterator;
			TRACE("Checking edge for being a tip "  << this->graph().str(tip));
			if (IsTip(tip)) {
				TRACE("Edge "  << this->graph().str(tip) << " judged to look like tip topologically");
				if (AdditionalCondition(tip)) {
                    TRACE("Additional checking");
                    removed++;
                    
                    // if tip was locked, we should not delete it
                    if (tip_lock.IsLocked(tip)){
                        TRACE("Tip " << this->graph().str(tip) << " was locked => can not remove it");
                        locked++;
                        continue;
                    }
                    
                    // removing only if stage is final
                    if (final_stage && CheckUniqueExtension(tip)){
                        TRACE("Edge " << this->graph().str(tip) << " has a unique extension");
                        RemoveTip(tip);
                        TRACE("Edge " << tip << " was removed"); 
                        continue;
                    }

                    // tricky condition -- not removing short edges with high coverage until the final stage
                    if (!TipHasLowRelativeCoverage(tip, final_stage)){
                        TRACE("Tip is covered well too much => not removing");
                        continue;
                    }
                    
                    // now we delete tip if we are not in the final stage
                    if (!final_stage){
					    TRACE("Edge "  << this->graph().str(tip) << " judged to be a tip with a very low coverage");
                        removed_with_check++;
                        RemoveTip(tip);
                        TRACE("Edge "  << tip << " removed as a tip");
                        continue;
                    }
                    
                    // now we are in the final stage 
                    // if the tip is covered very badly we delete it with no doubt
                    if (TipHasVeryLowRelativeCoverage(tip)){
					    TRACE("Edge "  << this->graph().str(tip) << " judged to be a tip with a very low coverage");
                        removed_with_check++;
                        RemoveTip(tip);
                        TRACE("Edge "  << tip << " removed as a tip");
                        continue;
                    }

                    // additional topology kind of check at the final stages
                    if (tipchecker.TipCanBeProjected(tip)){
                        TRACE("Edge "  << this->graph().str(tip) << " judged to be a tip");
                        removed_with_check++;
                        RemoveTip(tip);
                        TRACE("Edge "  << tip << " removed as a tip");
                        continue;
                    }
					TRACE("Edge "  << this->graph().str(tip) << " is not a tip");
				} else {
					TRACE("Edge "  << this->graph().str(tip) << " judged NOT to be a tip");
				}
			} else {
				TRACE("Edge "  << this->graph().str(tip) << " judged NOT to look like a tip topologically");
			}
			TRACE("Try to find next edge");
			TRACE("Use next edge");
		}
		TRACE("Tip clipping finished");
        
        DEBUG("REMOVED STATS " << removed_with_check << " " << removed);
        DEBUG("LOCKED " << locked);
        //DEBUG("REMOVED GOOD " << good_removed);
        //DEBUG("TOTAL GOOD " << good_total);

		Compressor<Graph> compressor(this->graph());
		compressor.CompressAllVertices();
	}

    void ClipTips(){
        ClipTips(false);   
    }


private:
	DECL_LOGGER("AdvancedTipClipper")
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
