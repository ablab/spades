//***************************************************************************
//* Copyright (c) 2015-2018 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/alignment/kmer_mapper.hpp"

namespace debruijn_graph {

const size_t MAX_SEQ_LENGTH = 1000;
inline Sequence RandomSequence(size_t length) {
    std::string result(length, 'A');
    for (size_t i = 0; i < length; i++) {
        result[i] = nucl((char)rand() % 4);
    }
    return Sequence(result);
}

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

template<typename T>
class RandomConstructor {

public:
    RandomConstructor(const char *name, T &value, size_t max_size)
            : name_(name), value_(value), max_size_(max_size) {
    }

    void Generate(size_t iterations, unsigned rand_seed = 100) {
        srand(rand_seed);
        while (iterations--)
            PerformRandomOperation();
        DEBUG("Generated random " << name_);
    }

private:
    virtual void PerformRandomOperation() = 0;
    const char *name_;

protected:
    T &value_;
    const size_t max_size_;
};

template<class Graph>
class RandomGraph : public RandomConstructor<Graph>, private RandomGraphAccessor<Graph> {

public:
    RandomGraph(Graph &graph, size_t max_size = 1000)
            : RandomConstructor<Graph>("graph", graph, max_size)
            , RandomGraphAccessor<Graph>(graph) {
    }

private:
    void AddRandomVertex() {
        this->value_.AddVertex();
    }

    void AddRandomEdge() {
        this->value_.AddEdge(this->GetRandomVertex(), this->GetRandomVertex(),
                             RandomSequence(rand() % MAX_SEQ_LENGTH + this->value_.k() + 1));
    }

    void RemoveRandomVertex() {
        this->value_.ForceDeleteVertex(this->GetRandomVertex());
    }

    void RemoveRandomEdge() {
        this->value_.DeleteEdge(this->GetRandomEdge());
    }

    void PerformRandomOperation() override {
        if (!this->value_.size())
            AddRandomVertex();
        else if (this->value_.SmartEdgeBegin().IsEnd()) {
            if (rand() % 2)
                AddRandomVertex();
            else
                AddRandomEdge();
        } else if (this->value_.size() > this->max_size_) {
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
};

template<class Index>
class RandomPairedIndex: public RandomConstructor<Index>, private RandomGraphAccessor<typename Index::Graph> {

public:
    RandomPairedIndex(Index &index, size_t max_size = 1000)
            : RandomConstructor<Index>("paired index", index, max_size)
            , RandomGraphAccessor<Graph>(index.graph()) {
    }

private:
    typedef typename Index::Graph Graph;
    typedef typename Index::Point Point;

    size_t max_size_;

    void AddRandomPoint() {
        using namespace omnigraph::de;
        const size_t MAX_DIST = 100;
        auto point = Point(DEDistance(rand() % MAX_DIST), DEWeight(1));
        this->value_.Add(this->GetRandomEdge(), this->GetRandomEdge(), point);
    }

    void RemoveRandomPoint() {
        size_t num = rand() % this->value_.size();
        auto i = omnigraph::de::pair_begin(this->value_);
        do {
            for (auto j : *i)
                if (num--) {
                    this->value_.Remove(i.first(), i.second(), j);
                    return;
                }
        } while (++i != omnigraph::de::pair_end(this->value_));
    }

    void RemoveRandomHistogram() {
        this->value_.Remove(this->GetRandomEdge(), this->GetRandomEdge());
    }

    void RemoveRandomEdgeInfo() {
        this->value_.Remove(this->GetRandomEdge());
    }

    void PerformRandomOperation() override {
        enum OpDistr {RemoveEdge = 0, RemoveHist = 1, RemovePoint = 4, AddPoint = 12};
        size_t dice = rand() % AddPoint;
        if (!this->value_.size())
            dice = AddPoint;
        else if (this->value_.size() >= max_size_)
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
};

template<typename Graph>
class RandomKmerMapper : public RandomConstructor<KmerMapper<Graph>> {

public:
    typedef KmerMapper<Graph> Type;

    RandomKmerMapper(Type &mapper, size_t max_size = 100)
            : RandomConstructor<Type>("kmer mapper", mapper, max_size) {
    }

private:
    void PerformRandomOperation() override {
        if (this->value_.size() < this->max_size_) {
            auto length = rand() % MAX_SEQ_LENGTH + this->value_.k() + 1;
            this->value_.RemapKmers(RandomSequence(length), RandomSequence(length));
        }
    }
};

}
