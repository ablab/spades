//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/mapping_path.hpp"

#include <vector>
#include <queue>
#include <map>

namespace cds_subgraphs {
using namespace debruijn_graph;

template<typename T>
void hash_combine(size_t &seed, T const &key) {
    std::hash<T> hasher;
    seed ^= hasher(key) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

using EdgePath = omnigraph::Path<EdgeId>;
using CodonSet = std::vector<Sequence>;

extern CodonSet STOP_CODONS;
extern CodonSet RC_STOP_CODONS;

//EdgeId + offset pair
using GraphPos = std::pair<EdgeId, size_t>;

struct FramedPos {
    EdgeId e;
    size_t offset;
    //"state" of the frame: 0 means that the codon is complete
    unsigned frame;

    FramedPos(EdgeId e_, size_t offset_, unsigned frame_) : e(e_), offset(offset_), frame(frame_) {}

    FramedPos(EdgeId e_, size_t offset_) : FramedPos(e_, offset_, 0) {}

    FramedPos() : e(EdgeId()), offset(size_t(-1)) {}

    bool good_frame() const {
        return frame == 0;
    }

    bool last(const Graph &g) const {
        VERIFY(offset < g.length(e));
        return offset == g.length(e) - 1;
    }

    FramedPos next() const {
        return FramedPos(e, offset + 1, (frame + 1) % 3);
    }

    FramedPos next(EdgeId neighbour) const {
        return FramedPos(neighbour, 0, (frame + 1) % 3);
    }

    bool operator!=(const FramedPos &other) const {
        return e != other.e || offset != other.offset || frame != other.frame;
    }

    bool operator==(const FramedPos &other) const {
        return !(*this != other);
    }

    size_t hash() const {
        size_t seed(0);
        hash_combine(seed, e);
        hash_combine(seed, offset);
        hash_combine(seed, frame);
        return seed;
    }
};

}

namespace std {
template<>
struct hash<cds_subgraphs::GraphPos> {
    size_t operator()(const cds_subgraphs::GraphPos &gpos) const {
        size_t seed(0);
        cds_subgraphs::hash_combine(seed, gpos.first);
        cds_subgraphs::hash_combine(seed, gpos.second);
        return seed;
    }
};

template<>
struct hash<cds_subgraphs::FramedPos> {
    size_t operator()(const cds_subgraphs::FramedPos &fp) const {
        return fp.hash();
    }
};
}

namespace cds_subgraphs {

class CodonFinder {
    const Graph &g_;
    GraphPos init_pos_;
    const CodonSet terminators_;
    std::queue<std::pair<FramedPos, FramedPos>> queue_;
    std::unordered_map<FramedPos, FramedPos> prev_;

    //Can be optimized
    //Last three nucleotides of a k+1-mer on position pos in the graph
    Sequence Codon(GraphPos pos) const {
        VERIFY(pos.second < g_.length(pos.first));
        return g_.EdgeNucls(pos.first).Subseq(pos.second + g_.k() - 2, pos.second + g_.k() + 1);
    }

    bool CodonInSet(GraphPos pos, const CodonSet& codons) const {
        return std::find(codons.begin(), codons.end(), Codon(pos)) != codons.end();
    }

    bool Terminate(FramedPos fpos) const {
        return fpos.good_frame() && CodonInSet(GraphPos(fpos.e, fpos.offset), terminators_);
    }

    void NextToQueue(FramedPos fpos);

public:
    CodonFinder(const Graph &g, GraphPos init_pos, const CodonSet& terminators):
            g_(g), init_pos_(init_pos), terminators_(terminators) {}

    std::vector<EdgePath> Go();

    //map from graph position to its corresponding path length
    //Note that those are "exclusive" coordinates of the ends of the path
    std::unordered_map<GraphPos, size_t> Terminates(const std::vector<EdgePath> &paths) const;

    std::set<EdgeId> Edges(const std::vector<EdgePath> &paths) const;

};

}
