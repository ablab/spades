//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


/*
 * concurrent_edge_algorithm.hpp
 *
 *  Created on: Sep 7, 2012
 *      Author: Alexander Opeykin (alexander.opeykin@gmail.com)
 */


#ifndef CONCURRENT_EDGE_ALGORITHM_HPP_
#define CONCURRENT_EDGE_ALGORITHM_HPP_

#include "order_and_law.hpp"
#include "devisible_tree.hpp"
#include "omni_tools.hpp"
#include "sequential_algorihtm_factory.hpp"
#include "concurrent_graph_component.hpp"
#include "concurrent_conjugate_graph_component.hpp"
#include "conjugate_vertex_glued_graph.hpp"
#include "component_algorithm_runner.hpp"
#include "perfcounter.hpp"

#include <memory>

namespace omnigraph {

template <class Graph>
class ConcurrentEdgeAlgorithm {

public:

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	typedef ConjugateVertexGluedGraph<Graph> GluedVertexGraph;
	typedef ConcurrentConjugateGraphComponent<Graph> ConjugateComponent;

	typedef ConcurrentGraphComponent<Graph> Component;
	typedef std::shared_ptr<Component> ComponentPtr;

	typedef SequentialAlgorihtmFactory<Component, EdgeId> Factory;
	typedef std::shared_ptr<Factory> FactoryPtr;

	typedef ComponentAlgorithmRunner<Graph, EdgeId> Runner;
	typedef std::shared_ptr<Runner> RunnerPtr;

	ConcurrentEdgeAlgorithm(const size_t nthreads, Graph& graph, FactoryPtr factory)
			: nthreads_(nthreads), graph_(graph), factory_(factory) {

		TRACE("Run in " << nthreads_ << " threads")

		GluedVertexGraph glued_vertex_graph (graph);
		DevisibleTree<GluedVertexGraph> tree (glued_vertex_graph);

		const size_t component_size = tree.GetSize() / nthreads;

		for (size_t thread = 0; thread < nthreads_; ++thread) {
			vector<VertexId> vertices;
			if (thread == nthreads_ - 1) {
				tree.SeparateVertices(vertices, tree.GetSize());
			} else {
				tree.SeparateVertices(vertices, component_size);
			}

			size_t actual_size = vertices.size();
			for (size_t i = 0; i < actual_size; ++i) {
				vertices.push_back(graph.conjugate(vertices[i]));
			}

			ComponentPtr ptr (
					new ConjugateComponent(
							graph,
							restricted::PeriodicIdDistributor(graph.GetGraphIdDistributor(),
								graph.GetGraphIdDistributor()->GetId(),
								nthreads
							),
							vertices.begin(),
							vertices.end()
					)
			);

			components_.push_back(ptr);
		}

		for (size_t i = 0; i < nthreads_; ++i) {
			RunnerPtr ptr (new Runner(*components_[i], factory->CreateAlgorithm(*components_[i])));
			runners_.push_back(ptr);
		}
	}

	// Comparator is used to define edge processing order. Not important by default.
	template <class Comparator = std::less<EdgeId>>
	void Run(Comparator comparator = std::less<EdgeId>()) {

		if (nthreads_ > 1) {
			VERIFY(graph_.AllHandlersThreadSafe());
		}

		vector<EdgeId> not_processed_edges_with_duplicates;

		#pragma omp parallel for num_threads(nthreads_)
		for (size_t i = 0; i < nthreads_; ++i) {
			auto it = components_[i]->SmartEdgeBegin(comparator);
			runners_[i]->Run(it);
		}

		for (size_t i = 0; i < nthreads_; ++i) {
			components_[i]->Synchronize();
		}

		restricted::PeriodicIdDistributor id_distributor(graph.GetGraphIdDistributor(),
				graph_.GetGraphIdDistributor()->GetId(), 1);

		ConjugateComponent all_graph_component(
				graph_, id_distributor, graph_.begin(), graph_.end());

		for (size_t i = 0; i < nthreads_; ++i) {
			components_[i]->GetEdgesGoingOutOfComponent(not_processed_edges_with_duplicates);
		}

		for (size_t i = 0; i < nthreads_; ++i) {
			runners_[i]->GetNotProcessedArguments(not_processed_edges_with_duplicates);
		}

		Runner border_runner(all_graph_component, factory_->CreateAlgorithm(all_graph_component));

		auto border_edge_iterator =	all_graph_component.SmartEdgeBegin(
				comparator, &not_processed_edges_with_duplicates);

		border_runner.Run(border_edge_iterator);

		// TODO: for debug only. remove.
		vector<EdgeId> border_not_processed_edges; // test vector. should have size = 0.
		border_runner.GetNotProcessedArguments(border_not_processed_edges);
		if (border_not_processed_edges.size() != 0) {
			INFO("WARNING: there are " << border_not_processed_edges.size() << " not processed in parallel");
		}
//		VERIFY(border_not_processed_edges.size() == 0);

		all_graph_component.Synchronize();
	}

private:
	const size_t nthreads_;
	Graph& graph_;
	FactoryPtr factory_;
	vector<ComponentPtr> components_;
	vector<RunnerPtr> runners_;
};


template <class Graph, class Algorithm>
class SequentialEdgeAlgorithm {
private:
	Graph &graph_;
	Algorithm &algorithm_;

public:

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	SequentialEdgeAlgorithm(Graph& graph, Algorithm &algorithm)
			: graph_(graph), algorithm_(algorithm) {
	}

	// Comparator is used to define edge processing order. Not important by default.
	template <class Comparator = std::less<EdgeId>>
	void Run(Comparator comparator = std::less<EdgeId>()) {

		algorithm_.Preprocessing();

		for (auto it = graph_.SmartEdgeBegin(comparator); !it.IsEnd(); ++it) {
			algorithm_.ProcessNext(*it);
		}

		algorithm_.Postprocessing();
	}

private:
	DECL_LOGGER("ConcurrentEdgeAlgorithm")
};
} // namespace omnigraph


#endif /* CONCURRENT_EDGE_ALGORITHM_HPP_ */
