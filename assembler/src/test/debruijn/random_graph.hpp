//***************************************************************************
//* Copyright (c) 2015-2018 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace debruijn_graph {

template<class Graph>
class RandomGraphAccessor {

public:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    RandomGraphAccessor(const Graph &graph) :
            graph_(graph) {}

    VertexId GetRandomVertex() const {
        size_t num = rand() % graph_.size();
        auto it = graph_.SmartVertexBegin();
        while (num--)
            ++it;
        return *it;
    }

    EdgeId GetRandomEdge() const {
        //Reservoir sampling
        EdgeId result = *(graph_.SmartEdgeBegin());
        size_t cur = 0;
        for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            cur++;
            if (rand() % cur == 0) {
                result = *it;
            }
        }
        return result;
    }

private:
    const Graph &graph_;
};

template<class Graph>
class RandomGraphConstructor : private RandomGraphAccessor<Graph> {
private:
    Graph &graph_;

    size_t max_size_;
    unsigned rand_seed_;

    Sequence GenerateRandomSequence(size_t length) {
        std::string result(length, 'A');
        for (size_t i = 0; i < length; i++) {
            result[i] = nucl((char)rand() % 4);
        }
        return Sequence(result);
    }

    void AddRandomVertex() {
        graph_.AddVertex();
    }

    void AddRandomEdge() {
        static const size_t MAX_SEQ_LENGTH = 1000;
        graph_.AddEdge(this->GetRandomVertex(), this->GetRandomVertex(),
                       GenerateRandomSequence(rand() % MAX_SEQ_LENGTH + graph_.k() + 1));
    }

    void RemoveRandomVertex() {
        graph_.ForceDeleteVertex(this->GetRandomVertex());
    }

    void RemoveRandomEdge() {
        graph_.DeleteEdge(this->GetRandomEdge());
    }

    void PerformRandomOperation() {
        if (!graph_.size())
            AddRandomVertex();
        else if (graph_.SmartEdgeBegin().IsEnd()) {
            if (rand() % 2)
                AddRandomVertex();
            else
                AddRandomEdge();
        } else if (graph_.size() > max_size_) {
            RemoveRandomVertex();
        } else {
            size_t tmp = rand() % 9;
            if (tmp == 0)
                AddRandomVertex();
            else if (tmp <= 6)
                AddRandomEdge();
            else
                RemoveRandomEdge();
        }
    }

public:
    RandomGraphConstructor(Graph &graph, size_t max_size) :
            RandomGraphAccessor<Graph>(graph), graph_(graph), max_size_(max_size) {
    }

    void Generate(size_t iterations, unsigned rand_seed = 100) {
        srand(rand_seed);
        while (iterations--)
            PerformRandomOperation();
        DEBUG("Generated graph of size " << graph_.size());
    }
};

template<class Index>
class RandomPairedIndexConstructor: private RandomGraphAccessor<typename Index::Graph> {
private:
    Index &index_;

    typedef typename Index::Graph Graph;
    typedef typename Index::Point Point;

    size_t max_size_;

    void AddRandomPoint() {
        using namespace omnigraph::de;
        const size_t MAX_DIST = 100;
        auto point = Point(DEDistance(rand() % MAX_DIST), DEWeight(1));
        index_.Add(this->GetRandomEdge(), this->GetRandomEdge(), point);
    }

    void RemoveRandomPoint() {
        size_t num = rand() % index_.size();
        auto i = omnigraph::de::pair_begin(index_);
        do {
            for (auto j : *i)
                if (num--) {
                    index_.Remove(i.first(), i.second(), j);
                    return;
                }
        } while (++i != omnigraph::de::pair_end(index_));
    }

    void RemoveRandomHistogram() {
        index_.Remove(this->GetRandomEdge(), this->GetRandomEdge());
    }

    void RemoveRandomEdgeInfo() {
        index_.Remove(this->GetRandomEdge());
    }

    void PerformRandomOperation() {
        enum OpDistr {RemoveEdge = 0, RemoveHist = 1, RemovePoint = 4, AddPoint = 12};
        size_t dice = rand() % AddPoint;
        if (!index_.size())
            dice = AddPoint;
        else if (index_.size() >= max_size_)
            dice = rand() % RemovePoint;
        if (dice <= RemoveEdge)
            RemoveRandomEdgeInfo();
        else if (dice <= RemoveHist)
            RemoveRandomHistogram();
        else if (dice <= RemovePoint)
            RemoveRandomPoint();
        else
            AddRandomPoint();
    }

public:
    RandomPairedIndexConstructor(Index &index, size_t max_size) :
            RandomGraphAccessor<Graph>(index.graph()), index_(index), max_size_(max_size) {
    }

    void Generate(size_t iterations, unsigned rand_seed = 100) {
        srand(rand_seed);
        while (iterations--)
            PerformRandomOperation();
        DEBUG("Generated index of size " << index_.size());
    }
};

}
