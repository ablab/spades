#include "distribution_extractor.hpp"

namespace path_extend {
    namespace fragment_model {

        vector<Fragment> LengthBasedFragmentExtractor::ExtractFragments() const {
            omnigraph::IterationHelper<Graph, EdgeId> edge_it_helper(g_);

            for (const auto edge: edge_it_helper) {

            }
        }

    }
}
