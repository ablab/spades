//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "assembly_graph/paths/path_utils.hpp"

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
                + "> junction_sequence <insert size> <int ids of edges> \n"
                + " Edges might be not consecutive, then will try to fix the path. \n"
                + " <insert size> specifies how many bp will be taken from first and last edges in the path. \n "
                + " flank size = IS - K - (cumulative length of intermediate edge).";
    }

    JunctionSequenceCommand()
            : LocalCommand<DebruijnEnvironment>("junction_sequence") {
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        size_t insert_size = std::stoll(args[1]);
        LOG("Insert size " << insert_size);

        stringstream ss;
        ss << "Provided path: ";
        vector<EdgeId> edges;
        for (size_t i = 2; i < args.size(); ++i) {
            EdgeId e = curr_env.finder().ReturnEdgeId(std::stoll(args[i]));
            edges.push_back(e);
            ss << curr_env.graph().str(e) << " ; ";
        }

        LOG(ss.str());

        debruijn_graph::MappingPathFixer<Graph> path_fixer(curr_env.graph());
        edges = path_fixer.TryFixPath(edges, insert_size);

        if (path_fixer.CheckContiguous(edges)) {
            LOG("Successfully fixed path");
        } else {
            LOG("Couldn't fix path!");
            return;
        }

        VERIFY(edges.size() > 1);
        size_t interm_length = CumulativeLength(curr_env.graph(), vector<EdgeId>(edges.begin() + 1, edges.end() - 1));
        if (insert_size < curr_env.graph().k() + interm_length) {
            LOG("Intermediate sequence too long");
            return;
        }

        size_t flank_length = insert_size - curr_env.graph().k() - interm_length;
        LOG("Flank length " << flank_length);

        if (curr_env.graph().length(edges.front()) < flank_length || curr_env.graph().length(edges.back()) < flank_length) {
            LOG("Flanking path edges can not be shorter than flank length!");
            return;
        }

        Path<EdgeId> path(edges, curr_env.graph().length(edges.front()) - flank_length, flank_length);

        LOG("Outputting sequence:");
        LOG(PathSequence(curr_env.graph(), path));
    }

private:
    DECL_LOGGER("JunctionSequenceCommand")
    ;
};

}
