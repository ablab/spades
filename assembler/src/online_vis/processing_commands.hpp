//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "environment.hpp"
#include "command.hpp"
#include "errors.hpp"
#include "argument_list.hpp"
#include "omni/tip_clipper.hpp"
#include "genomic_quality.hpp"

namespace online_visualization {

class ClipTipsCommand : public NewLocalCommand<DebruijnEnvironment> {

public:
    string Usage() const {
        string answer;
        answer = answer + "Command `clip_tips` \n" + "Usage:\n"
                + "> clip_tips <length>\n" + " This command clips tips.\n"
                + " If length is not specified, "
                + "it will be counted from global settings ";
        return answer;
    }

    ClipTipsCommand()
            : NewLocalCommand<DebruijnEnvironment>("clip_tips", 1) {
    }

private:
    /*virtual*/ void InnerExecute(DebruijnEnvironment& curr_env,
                 const vector<string>& args) const {
        size_t length = GetInt(args[1]);
        debruijn_graph::EdgeQuality<Graph, Index> edge_qual(curr_env.graph(),
                                                            curr_env.index(),
                                                            curr_env.kmer_mapper(),
                                                            curr_env.genome());

        shared_ptr<func::Predicate<EdgeId>> condition = make_shared<AlwaysTrue<EdgeId>>();
        if (args.size() > 1 && (args[1] == "Y" || args[1] == "y")) {
            cout << "Trying to activate genome quality condition" << endl;
            if (curr_env.genome().size() == 0) {
                cout << "No reference was provided!!!" << endl;
            } else {
                cout << "Genome quality condition will be used" << endl;
//                condition = make_shared<make_shared<debruijn_graph::ZeroQualityCondition<Graph, Index>>(curr_env.graph(), edge_qual);
                condition = make_shared<func::AdaptorPredicate<EdgeId>>(boost::bind(&debruijn_graph::EdgeQuality<Graph, Index>::IsZeroQuality, boost::ref(edge_qual), _1));
            }
        }
        omnigraph::ClipTips(curr_env.graph(), length, condition);
    }
};
}
