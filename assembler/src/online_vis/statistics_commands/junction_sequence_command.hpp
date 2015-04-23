#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "omni/omni_utils.hpp"
#include "omni/path_processor.hpp"
#include "utils.hpp"

namespace online_visualization {
class JunctionSequenceCommand : public LocalCommand<DebruijnEnvironment> {

protected:
    size_t MinArgNumber() const {
        return 3;
    }

    bool CheckCorrectness(const vector<string>& args) const {
        return CheckEnoughArguments(args);
    }

public:
    string Usage() const {
        return string() + "Command `junction_sequence` \n" + " Usage:\n"
                + "> junction_sequence <flank length> <int ids of edges> \n"
                + " Edges might be not consecutive, then will try to fix the path. \n"
                + " <flank length> specifies how many bp will be taken from first and last edges in the path.";
    }

    JunctionSequenceCommand()
            : LocalCommand<DebruijnEnvironment>("junction_sequence") {
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        size_t flank_legth = boost::lexical_cast<size_t>(args[1]);
        LOG("Flank length " << flank_legth);

        vector<EdgeId> edges;
        for (size_t i = 2; i < args.size(); ++i) {
            edges.push_back(curr_env.finder().ReturnEdgeId(boost::lexical_cast<size_t>(args[i])));
        }

        if (curr_env.graph().length(edges.front()) < flank_legth || curr_env.graph().length(edges.back()) < flank_legth) {
            LOG("Flanking path edges can not be shorter than flank length!");
            return;
        }

        Path<EdgeId> path(edges, curr_env.graph().length(edges.front()) - flank_legth, flank_legth);

        MappingPathFixer<Graph> path_fixer(curr_env.graph());
        path = path_fixer.TryFixPath(path);

        if (path_fixer.CheckContiguous(path.sequence())) {
            LOG("Successfully fixed path, outputting sequence:");
            LOG(PathSequence(curr_env.graph(), path));
        } else {
            LOG("Couldn't fix path!");
        }
    }

private:
    DECL_LOGGER("JunctionSequenceCommand")
    ;
};
}
