#pragma once

#include "standard_base.hpp"
#include "comparison_utils.hpp"

namespace cap {

typedef Sequence Genome;
typedef map<Genome, Range> GeneCoordinates;

//range of nucleotide positions; true if main strand
typedef pair<Range, bool> Pos;
typedef size_t GenomeId;

typedef map<GenomeId, Pos> GenePosition;

struct GeneCollection {
    vector<Genome> genomes;
    vector<GenePosition> gene_positions;
};

typedef vector<Range> Coordinates;

template<class gp_t>
class CoordinatesUpdater {
    typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef MappingPath<EdgeId> Path;

    const gp_t& gp_;

    Range NewCoords(const Path& path, Range coord) const {
        size_t cumm_length = 0;
        Range answer;
        int i = 0;
        for (; i < path.size() && path[i].second.initial_range.end_pos <= coord.start_pos; ++i) {
            cumm_length += gp_.g.length(path[i].first);
        }
        VERIFY(i < path.size());
        VERIFY(path[i].second.initial_range.end_pos > coord.start_pos);
        VERIFY(path[i].second.initial_range.start_pos <= coord.start_pos);
        answer.start_pos = cumm_length;
        for (; i < path.size() && path[i].second.initial_range.end_pos < coord.end_pos; ++i) {
            cumm_length += gp_.g.length(path[i].first);
        }
        VERIFY(i < path.size());
        VERIFY(path[i].second.initial_range.end_pos >= coord.end_pos);
        VERIFY(path[i].second.initial_range.start_pos < coord.end_pos);
        answer.end_pos = cumm_length + gp_.g.length(path[i].first);
        return answer;
    }

 public:
    CoordinatesUpdater(const gp_t& gp) : gp_(gp) {

    }

    Coordinates Update(const Genome& genome, const Coordinates& coords) const {
        Coordinates answer;
        auto mapper = *MapperInstance(gp_);
        auto mapping_path = mapper.MapSequence(genome);
        FOREACH(Range r, coords) {
            answer.push_back(NewCoords(mapping_path, r));
        }
        return answer;
    }
};

template<class Graph>
void WriteGeneLocality(const GeneCoordinates& coords) {

}

}
