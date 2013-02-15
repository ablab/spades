#pragma once

#include "standard_base.hpp"
#include "simple_tools.hpp"
#include "comparison_utils.hpp"

namespace cap {

typedef Sequence Genome;
typedef map<Genome, Range> GeneCoordinates;

//range of nucleotide positions; true if main strand
typedef pair<Range, bool> Pos;
typedef size_t GenomeId;

typedef multimap<GenomeId, Pos> GenePositions;
typedef vector<Range> Coordinates;

Range NuclToKmerRange(Range nucl_range, size_t k) {
    return Range(std::max(0, int(nucl_range.start_pos) - int(k)),
                 nucl_range.end_pos);
}

Range KmerToNuclRange(Range kmer_range, size_t k) {
    return Range(min(kmer_range.start_pos + k, kmer_range.end_pos - 1),
                 kmer_range.end_pos);
}

Range OppositeStrandNuclCoord(Range nucl_coord, size_t genome_length) {
    return Range(genome_length - 1 - nucl_coord.end_pos, genome_length - 1 -nucl_coord.start_pos);
}

struct GeneCollection {
    map<GenomeId, Genome> genomes;
    vector<GenePositions> gene_positions;

 private:

    void Normalize(GenomeId genome_id, Pos& pos) {
        if (!pos.second.second) {
//            pos.
        }
    }

 public:
    template<class gp_t>
    void Update(const gp_t& gp) {
        CoordinatesUpdater<gp_t> updater(gp);
        FOREACH(GenomeId genome_id, key_set(genomes)) {
            Coordinates gene_coords;
            vector<size_t> gene_ids;
            for (size_t i = 0; i < gene_positions.size(); ++i) {
                std::multimap m;

                FOREACH(Pos pos, get_all(gene_positions[i], genome_id)) {
                    gene_ids.push_back(i);
                    VERIFY(pos.second);
                    gene_coords.push_back(pos.first);
                }
            }
            auto updated = updater.Update(genomes[genome_id], Coordinates);
            genomes[genome_id] = updated.first;

            FOREACH(GenePositions gene_poss, gene_positions) {
                gene_poss.clear();
            }

            for (size_t j = 0; j < gene_ids.size(); ++j) {
                gene_positions[gene_ids[j]].insert(make_pair(genome_id, make_pair(updated.second[j], true)));
            }
        }
    }
};

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

    pair<Genome, Coordinates> Update(const Genome& genome, const Coordinates& coords) const {
        pair<Genome, Coordinates> answer;
        auto mapper = *MapperInstance(gp_);
        auto mapping_path = mapper.MapSequence(genome);
        FOREACH(Range r, coords) {
            answer.second.push_back(NewCoords(mapping_path, r));
        }
        auto read_corrector = GraphReadCorrectorInstance(gp_.g, mapper);
        answer.first = read_corrector.Modify(genome);
        return answer;
    }
};

template<class Graph>
void WriteGeneLocality(const GeneCoordinates& coords) {

}

}
