#pragma once
#include "common/assembly_graph/core/graph.hpp"

#include <boost/optional.hpp>
#include <string>
#include <vector>

namespace sequence_corrector {

using Graph = debruijn_graph::Graph;
using EdgeId = Graph::EdgeId;
using Path = std::vector<EdgeId>;
class FillerChooser {
public:
    virtual boost::optional<std::string> operator()(std::vector<Path> const & paths, std::string const & ref) const = 0;
    virtual ~FillerChooser() = default;
};

class AlignerFiller : public FillerChooser {
    Graph const & graph;
    double score_domination_coeff;
public:
    AlignerFiller(Graph const & graph, double score_domination_coeff) 
        : graph(graph)
        , score_domination_coeff(score_domination_coeff)
    {}

    boost::optional<std::string> operator()(std::vector<Path> const & paths, std::string const & ref) const override;

    bool IsDominantScore(size_t dominator, size_t other) const noexcept;
};

void CompressTest();
void ConsensusTest();

} // namespace sequence_corrector
