#pragma once

#include "logger/logger.hpp"
#include "graph_iterators.hpp"
#include "graph_processing_algorithm.hpp"

namespace omnigraph {

//todo add conjugate filtration
template<class Graph, class ElementType>
class AlgorithmRunner {
    const Graph& g_;

    template<class Algo, class It>
    bool ProcessBucket(Algo& algo, It begin, It end) {
        bool changed = false;
        for (auto it = begin; it != end; ++it) {
            changed |= algo.Process(*it);
        }
        return changed;
    }

public:

    const Graph& g() const {
        return g_;
    }

    AlgorithmRunner(Graph& g)
            : g_(g) {

    }

    template<class Algo, class ItVec>
    bool RunFromChunkIterators(Algo& algo, const ItVec& chunk_iterators) {
        DEBUG("Running from " << chunk_iterators.size() - 1 << "chunks");
        VERIFY(chunk_iterators.size() > 1);
        bool changed = false;
        #pragma omp parallel for schedule(guided) reduction(|:changed)
        for (size_t i = 0; i < chunk_iterators.size() - 1; ++i) {
            changed |= ProcessBucket(algo, chunk_iterators[i], chunk_iterators[i + 1]);
        }
        DEBUG("Finished");
        return changed;
    }
private:
    DECL_LOGGER("AlgorithmRunner")
    ;
};

template<class Graph, class ElementType>
class TwoStepAlgorithmRunner {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    const bool filter_conjugate_;
    std::vector<std::vector<ElementType>> elements_of_interest_;

    template<class Algo>
    bool ProcessBucket(Algo& algo, const std::vector<ElementType>& bucket, size_t idx_offset) const {
        bool changed = false;
        for (ElementType el : bucket) {
            changed |= algo.Process(el, idx_offset++);
        }
        return changed;
    }

    template<class Algo>
    bool Process(Algo& algo) const {
        std::vector<size_t> cumulative_bucket_sizes;
        cumulative_bucket_sizes.push_back(0);
        for (const auto& bucket : elements_of_interest_) {
            cumulative_bucket_sizes.push_back(cumulative_bucket_sizes.back() + bucket.size());
        }
        DEBUG("Preparing for processing");
        algo.PrepareForProcessing(cumulative_bucket_sizes.back());
        bool changed = false;
        DEBUG("Processing buckets");
        #pragma omp parallel for schedule(guided) reduction(|:changed)
        for (size_t i = 0; i < elements_of_interest_.size(); ++i) {
            changed |= ProcessBucket(algo, elements_of_interest_[i], cumulative_bucket_sizes[i]);
        }
        return changed;
    }

    template<class Algo>
    void CountElement(Algo& algo, ElementType el, size_t bucket) {
        if (filter_conjugate_ && g_.conjugate(el) < el)
            return;
        if (algo.IsOfInterest(el)) {
            INFO("Element " << g_.str(el) << " is of interest");
            elements_of_interest_[bucket].push_back(el);
        } else {
            INFO("Element " << g_.str(el) << " is not interesting");
        }
    }

    template<class Algo, class It>
    void CountAll(Algo& algo, It begin, It end, size_t bucket) {
        for (auto it = begin; !(it == end); ++it) {
            CountElement(algo, *it, bucket);
        }
    }

public:

    const Graph& g() const {
        return g_;
    }

    //conjugate elements are filtered based on ids
    //should be used only if both conjugate elements are simultaneously either interesting or not
    //fixme filter_conjugate is redundant
    TwoStepAlgorithmRunner(Graph& g, bool filter_conjugate)
            : g_(g),
              filter_conjugate_(filter_conjugate) {

    }

    template<class Algo, class ItVec>
    bool RunFromChunkIterators(Algo& algo, const ItVec& chunk_iterators) {
        DEBUG("Started running from " << chunk_iterators.size() - 1 << " chunks");
        VERIFY(algo.ShouldFilterConjugate() == filter_conjugate_);
        VERIFY(chunk_iterators.size() > 1);
        elements_of_interest_.clear();
        elements_of_interest_.resize(chunk_iterators.size() - 1);
        DEBUG("Searching elements of interest");
        #pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < chunk_iterators.size() - 1; ++i) {
            CountAll(algo, chunk_iterators[i], chunk_iterators[i + 1], i);
        }
        DEBUG("Processing");
        return Process(algo);
    }

//    template<class Algo, class It>
//    void RunFromIterator(Algo& algo, It begin, It end) {
//        RunFromChunkIterators(algo, std::vector<It> { begin, end });
//    }
private:
    DECL_LOGGER("TwoStepAlgorithmRunner")
    ;
};

template<class Graph, class ElementType>
class SemiParallelAlgorithmRunner {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;

public:

    const Graph& g() const {
        return g_;
    }

    SemiParallelAlgorithmRunner(Graph& g)
            : g_(g) {

    }

    template<class Algo, class ItVec>
    bool RunFromChunkIterators(Algo& algo, const ItVec& chunk_iterators) {
        VERIFY(chunk_iterators.size() > 1);
        std::vector<std::vector<ElementType>> of_interest(chunk_iterators.size() - 1);

        #pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < chunk_iterators.size() - 1; ++i) {
            for (auto it = chunk_iterators[i], end = chunk_iterators[i + 1]; it != end; ++it) {
                ElementType t = *it;
                if (algo.IsOfInterest(t)) {
                    of_interest[i].push_back(t);
                }
            }
        }

        auto it = SmartSetIterator<Graph, ElementType>(g_);
        for (auto& chunk : of_interest) {
            it.insert(chunk.begin(), chunk.end());
        }
        bool changed = false;
        for (; !it.IsEnd(); ++it) {
            changed |= algo.Process(*it);
        }
        return changed;
    }

private:
    DECL_LOGGER("SemiParallelAlgorithmRunner")
    ;
};

//todo generalize to use for other algorithms if needed
template<class Graph>
class SemiParallelEdgeRemovingAlgorithm {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    Graph& g_;
    shared_ptr<func::Predicate<EdgeId>> condition_;
    EdgeRemover<Graph> edge_remover_;

public:
    SemiParallelEdgeRemovingAlgorithm(Graph& g,
                                      shared_ptr<func::Predicate<EdgeId>> condition,
                                      boost::function<void(EdgeId)> removal_handler = 0) :
            g_(g), condition_(condition), edge_remover_(g, removal_handler) {
    }

    bool IsOfInterest(EdgeId e) const {
        return condition_->Check(e);
    }

    bool Process(EdgeId e) {
        edge_remover_.DeleteEdge(e);
        return true;
    }
};

template<class Graph, class AlgoRunner, class Algo>
bool RunVertexAlgorithm(Graph& g, AlgoRunner& runner, Algo& algo, size_t chunk_cnt) {
    return runner.RunFromChunkIterators(algo, ParallelIterationHelper<Graph>(g).VertexChunks(chunk_cnt));
}

template<class Graph, class AlgoRunner, class Algo>
bool RunEdgeAlgorithm(Graph& g, AlgoRunner& runner, Algo& algo, size_t chunk_cnt) {
    return runner.RunFromChunkIterators(algo, ParallelIterationHelper<Graph>(g).EdgeChunks(chunk_cnt));
}

}
