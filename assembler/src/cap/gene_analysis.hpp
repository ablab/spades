#pragma once

#include "standard_base.hpp"
#include "comparison_utils.hpp"

namespace cap {

typedef Sequence Genome;
typedef map<Genome, Range> GeneCoordinates;
typedef vector<Range> Coordinates;

template<class gp_t>
class CoordinatesUpdater {
    typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const gp_t& gp_;

    Range NewCoords(MappingPath<EdgeId> path, Range coord) {

    }

 public:
    CoordinatesUpdater(const gp_t& gp) : gp_(gp) {

    }

    Coordinates Update(const Genome& genome, const Coordinates& coords) const {
        auto mapper = *MapperInstance(gp_);
        auto mapping_path = mapper.MapSequence(genome);
    }
};

template<class Graph>
void WriteGeneLocality(const GeneCoordinates& coords) {

}

}
