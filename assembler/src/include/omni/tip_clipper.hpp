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

#include "omni_utils.hpp"
#include "sequential_algorithm.hpp"
#include "sequential_algorihtm_factory.hpp"
#include "xmath.h"
#include "concurrent_conjugate_graph_component.hpp"
#include "concurrent_edge_algorithm.hpp"
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
template<class Graph>
class AbstractTipClipper
		: public SequentialAlgorithm<typename Graph::EdgeId>,
		  private boost::noncopyable {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph &graph_;
	size_t removed_;

public:
	const size_t max_tip_length_;

private:	
    boost::function<void(EdgeId)> removal_handler_;

protected:

	/**
	 * Create TipClipper with specified parameters. Those parameters could probably be replaced later with
	 * certain generic checker class.
	 */
	AbstractTipClipper(
			Graph &graph,
			size_t max_tip_length,
			boost::function<void(EdgeId)> removal_handler = 0,
			boost::function<double(EdgeId)> qual_f = 0)
				: graph_(graph),
				  removed_(0),
				  max_tip_length_(max_tip_length),
				  removal_handler_(removal_handler) {

	}

	const Graph& graph() const{
		return graph_;
	}

    Graph& graph(){
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

			// for debug. It should happen only if CanCompressVertex
			// invalidate graph.
//			VERIFY(graph_.IsValid());

			graph_.CompressVertex(splitVertex);
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

	virtual void RemoveTip(EdgeId tip) {
		VertexId start = graph_.EdgeStart(tip);
		VertexId end = graph_.EdgeEnd(tip);
		if (removal_handler_) {
			removal_handler_(tip);
		}
		graph_.DeleteEdge(tip);
		ProcessVertex(start);
		ProcessVertex(end);
	}

	virtual bool TryToRemoveTip(EdgeId tip) {
//		if (!graph_.IsInternalSafe(tip)) {
//			// for algorithm to process this edge sequently.
//			graph_.InvalidateComponent();
//		}
//
//		if (graph_.IsValid()) {
			RemoveTip(tip);
			TRACE("Edge removed");
//			return true;
//		} else {
//			TRACE("Component is invalid. " << "Edge "  << graph_.str(tip) << " can not be removed in parallel.");
//			return false;
//		}
			return true;
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
    virtual void ProcessNext(const EdgeId& tip) {
		TRACE("Checking edge for being tip "  << this->graph().str(tip));
		if (this->IsTip(tip)) {
			TRACE("Edge "  << this->graph().str(tip) << " judged to look like a tip topologically");
			if (AdditionalCondition(tip)) {
				TRACE("Edge "  << this->graph().str(tip) << " judged to be a tip");

				if (this->TryToRemoveTip(tip)) {
					removed_++;
				}

			} else {
				TRACE("Edge "  << this->graph().str(tip) << " judged NOT to be tip");
			}
		} else {
			TRACE("Edge "  << this->graph().str(tip) << " judged NOT to look like tip topologically");
		}
    }

    virtual void Preprocessing() {
    	TRACE("Tip clipping (maximal corruption mode) started");
    }

    virtual void Postprocessing() {
    	TRACE("Tip clipping finished");
    	DEBUG("REMOVED " << removed_);
    }

private:
	DECL_LOGGER("AbstractTipClipper")
};

template<class Graph>
class DefaultTipClipper: public AbstractTipClipper<Graph> {

typedef typename Graph::EdgeId EdgeId;
typedef typename Graph::VertexId VertexId;
typedef AbstractTipClipper<Graph> base;
        
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

	DefaultTipClipper(
			Graph &graph,
			size_t max_tip_length,
			size_t max_coverage,
			double max_relative_coverage,
			boost::function<void(EdgeId)> removal_handler = 0,
			boost::function<double(EdgeId)> qual_f = 0)
				: base(graph, max_tip_length, removal_handler, qual_f),
				  max_coverage_(max_coverage),
				  max_relative_coverage_(max_relative_coverage) {
	}

private:
    DECL_LOGGER("DefaultTipClipper")
};

template<class Graph>
class DefaultTipClipperFactory : public SequentialAlgorihtmFactory<Graph, typename Graph::EdgeId> {

public:
	typedef typename Graph::EdgeId EdgeId;
	typedef SequentialAlgorihtmFactory<Graph, EdgeId> Base;
	typedef typename Base::AlgorithmPtr AlgorithmPtr;

	DefaultTipClipperFactory(
			size_t max_tip_length,
			size_t max_coverage,
			double max_relative_coverage,
			boost::function<void(EdgeId)> removal_handler = 0,
			boost::function<double(EdgeId)> qual_f = 0)
				: 	max_tip_length_(max_tip_length),
				  	max_coverage_(max_coverage),
				  	max_relative_coverage_(max_relative_coverage),
				  	removal_handler_(removal_handler),
				  	qual_f_(qual_f) {
	}

	virtual AlgorithmPtr CreateAlgorithm(Graph& graph) {
		AlgorithmPtr ptr(
				new DefaultTipClipper<Graph>(
						graph,
						max_tip_length_,
						max_coverage_,
						max_relative_coverage_,
						removal_handler_,
						qual_f_));
		return ptr;
	}

private:
	size_t max_tip_length_;
	size_t max_coverage_;
	double max_relative_coverage_;
	boost::function<void(EdgeId)> removal_handler_;
	boost::function<double(EdgeId)> qual_f_;
};



template<class Graph>
class AdvancedTipClipper: public AbstractTipClipper<Graph> {

private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef AbstractTipClipper<Graph> base;
    TipLock<EdgeId> tip_lock_;

    const size_t max_coverage_;
    const double max_relative_coverage_;
    const size_t max_iterations_;
    const size_t max_levenshtein_;
    const size_t max_ec_length_;

    size_t removed_;
    size_t removed_with_check_;
    size_t locked_;
    bool final_stage_;

    TipChecker<Graph> tipchecker_;

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
        bool backward = this->IsTip(this->graph().EdgeStart(tip));
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

    bool TipHasLowRelativeCoverage(EdgeId tip) const {
        double min_covlen = MinCompetitorCoverage(tip);  
        if (final_stage_ && this->graph().length(tip) < this->graph().k() / 2)
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

	AdvancedTipClipper(
			Graph &graph,
			size_t max_tip_length,
			size_t max_coverage,
			double max_relative_coverage,
			size_t max_iterations,
			size_t max_levenshtein,
			size_t max_ec_length,
			boost::function<void(EdgeId)> removal_handler = 0,
			bool final_stage = false)
				: base(graph, max_tip_length, removal_handler),
				  max_coverage_(max_coverage),
				  max_relative_coverage_(max_relative_coverage),
				  max_iterations_(max_iterations),
				  max_levenshtein_(max_levenshtein),
				  max_ec_length_(max_ec_length),
				  final_stage_(final_stage),
				  tipchecker_(graph, tip_lock_, max_iterations_, max_levenshtein_, max_tip_length, max_ec_length_) {

        removed_ = 0;
        removed_with_check_ = 0;
        locked_ = 0;
	}


	virtual void Preprocessing() {
		TRACE("Tip clipping started");
	}

	virtual void Postprocessing() {
		TRACE("Tip clipping finished");
        DEBUG("REMOVED STATS " << removed_with_check_ << " " << removed_);
        DEBUG("LOCKED " << locked_);
	}

    // Method deletes tips from the graph carefully, its work depends on the number of simplification iteration
	virtual void ProcessNext(const EdgeId& tip) {
        TRACE("Use next edge");
		TRACE("Checking edge for being a tip "  << this->graph().str(tip));

		if (this->IsTip(tip)) {
			TRACE("Edge "  << this->graph().str(tip) << " judged to look like tip topologically");
			if (AdditionalCondition(tip)) {
				TRACE("Additional checking");
				removed_++;

				// if tip was locked, we should not delete it
				if (tip_lock_.IsLocked(tip)){
					TRACE("Tip " << this->graph().str(tip) << " was locked => can not remove it");
					locked_++;
					return;
				}

				// removing only if stage is final
				if (final_stage_ && CheckUniqueExtension(tip)){
					TRACE("Edge " << this->graph().str(tip) << " has a unique extension");
					this->TryToRemoveTip(tip);
					return;
				}

				// tricky condition -- not removing short edges with high coverage until the final stage
				if (!TipHasLowRelativeCoverage(tip)){
					TRACE("Tip is covered well too much => not removing");
					return;
				}

				// now we delete tip if we are not in the final stage
				if (!final_stage_){
					TRACE("Edge "  << this->graph().str(tip) << " judged to be a tip with a very low coverage");
					if (this->TryToRemoveTip(tip)) {
						removed_with_check_++;
					}
					return;
				}

				// now we are in the final stage
				// if the tip is covered very badly we delete it with no doubt
				if (TipHasVeryLowRelativeCoverage(tip)){
					TRACE("Edge "  << this->graph().str(tip) << " judged to be a tip with a very low coverage");
					if (this->TryToRemoveTip(tip)) {
						removed_with_check_++;
					}
					return;
				}

				// additional topology kind of check at the final stages
				if (tipchecker_.TipCanBeProjected(tip)){
					TRACE("Edge "  << this->graph().str(tip) << " judged to be a tip");
					if (this->TryToRemoveTip(tip)) {
						removed_with_check_++;
					}
					return;
				}
				TRACE("Edge "  << this->graph().str(tip) << " is not a tip");
			} else {
				TRACE("Edge "  << this->graph().str(tip) << " judged NOT to be a tip");
			}
		} else {
			TRACE("Edge "  << this->graph().str(tip) << " judged NOT to look like a tip topologically");
		}
	}


private:
	DECL_LOGGER("AdvancedTipClipper")
};


template <class Graph>
class AdvancedTipClipperFactory : public SequentialAlgorihtmFactory<Graph, typename Graph::EdgeId> {

public:
	typedef typename Graph::EdgeId EdgeId;
	typedef SequentialAlgorihtmFactory<Graph, EdgeId> Base;
	typedef typename Base::AlgorithmPtr AlgorithmPtr;

	AdvancedTipClipperFactory (
			size_t max_tip_length,
			size_t max_coverage,
			double max_relative_coverage,
			size_t max_iterations,
			size_t max_levenshtein,
			size_t max_ec_length,
			boost::function<void(EdgeId)> removal_handler = 0,
			bool final_stage = false)
				: 	max_tip_length_(max_tip_length),
				  	max_coverage_(max_coverage),
				  	max_relative_coverage_(max_relative_coverage),
				  	max_iterations_(max_iterations),
				  	max_levenshtein_(max_levenshtein),
				  	max_ec_length_(max_ec_length),
				  	removal_handler_(removal_handler),
				  	final_stage_(final_stage) {
	}

	virtual AlgorithmPtr CreateAlgorithm(Graph& graph) {
		AlgorithmPtr ptr(
			new AdvancedTipClipper<Graph>(
					graph,
					max_tip_length_,
					max_coverage_,
				  	max_relative_coverage_,
				  	max_iterations_,
				  	max_levenshtein_,
				  	max_ec_length_,
				  	removal_handler_,
				  	final_stage_));

		return ptr;
	}

private:
	size_t max_tip_length_;
	size_t max_coverage_;
	double max_relative_coverage_;
	size_t max_iterations_;
	size_t max_levenshtein_;
	size_t max_ec_length_;
	boost::function<void(EdgeId)> removal_handler_;
	bool final_stage_;

};


template<class Graph>
class TopologyTipClipper : public AbstractTipClipper<Graph> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef AbstractTipClipper<Graph> base;

	size_t uniqueness_length_;
	size_t plausibility_length_;
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

	bool CheckPlausibleAlternative(const vector<EdgeId> &rivals) const {
		for(auto it = rivals.begin(); it != rivals.end(); ++it) {
			if(this->graph().length(*it) >= plausibility_length_)
				return true;
		}
		return false;
	}

protected:

	/*virtual*/
	bool AdditionalCondition(EdgeId tip) const {
		vector<EdgeId> unique_path;
		vector<EdgeId> rivals;
		if (this->graph().IsDeadEnd(this->graph().EdgeEnd(tip))
				&& this->graph().CheckUniqueIncomingEdge(this->graph().EdgeStart(tip))) {
			unique_path = unique_path_finder_.UniquePathBackward(
					this->graph().GetUniqueIncomingEdge(this->graph().EdgeStart(tip)));
			rivals = this->graph().OutgoingEdges(this->graph().EdgeStart(tip));
		} else if (this->graph().IsDeadStart(this->graph().EdgeStart(tip))
				&& this->graph().CheckUniqueOutgoingEdge(this->graph().EdgeEnd(tip))) {
			unique_path = unique_path_finder_.UniquePathForward(
					this->graph().GetUniqueOutgoingEdge(this->graph().EdgeEnd(tip)));
			rivals = this->graph().IncomingEdges(this->graph().EdgeEnd(tip));
		}
		return CummulativeLength(this->graph(), unique_path) >= uniqueness_length_ && CheckPlausibleAlternative(rivals);
	}

public:
	TopologyTipClipper(
			Graph &graph,
			size_t max_tip_length,
			size_t uniqueness_length,
			size_t plausibility_length,
			boost::function<void(EdgeId)> removal_handler = 0)
				: base(graph, max_tip_length, removal_handler),
				  uniqueness_length_(uniqueness_length),
				  plausibility_length_(plausibility_length),
				  unique_path_finder_(graph) {
	}

private:
	DECL_LOGGER("TopologyTipClipper")
};

template<class Graph>
void RunTopologyTipClipper(Graph &graph,
		size_t max_tip_length,
		size_t uniqueness_length,
		size_t plausibility_length,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0) {
	ConcurrentConjugateGraphComponent<Graph> all_graph_component (
			graph,
			restricted::PeriodicIdDistributor(
				restricted::GlobalIdDistributor::GetInstance()->GetId(),
				1
			),
			graph.begin(),
			graph.end()
	);
	omnigraph::TopologyTipClipper<ConcurrentConjugateGraphComponent<Graph>> tc(all_graph_component, max_tip_length, uniqueness_length, plausibility_length, removal_handler);
	SequentialEdgeAlgorithm<ConcurrentConjugateGraphComponent<Graph>, omnigraph::TopologyTipClipper<ConcurrentConjugateGraphComponent<Graph>>>(all_graph_component, tc).Run();
	all_graph_component.Synchronize();
}


} // namespace omnigraph

#endif /* TIP_CLIPPER_HPP_ */
