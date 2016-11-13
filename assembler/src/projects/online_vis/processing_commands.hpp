//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "environment.hpp"
#include "command.hpp"
#include "errors.hpp"
#include "argument_list.hpp"
#include "assembly_graph/graph_support/genomic_quality.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"

namespace online_visualization {

class ClipTipsCommand : public NewLocalCommand<DebruijnEnvironment> {

public:
    string Usage() const {
        string answer;
        answer = answer + "Command `clip_tips` \n" + "Usage:\n"
                + "> clip_tips <length> (Y/y)\n" + " This command clips tips.\n"
                + " If length is not specified, "
                + "it will be counted from global settings. "
                + "If second argument Y/y is specified then genomic edges will be retained.";
        return answer;
    }

    ClipTipsCommand()
            : NewLocalCommand<DebruijnEnvironment>("clip_tips", 0) {
    }

private:
    /*virtual*/ void InnerExecute(DebruijnEnvironment& curr_env,
                 const vector<string>& args) const {
        size_t length = 0;
        if(args.size() > 0) {
            length = GetInt(args[1]);
        } else {
            length = curr_env.edge_length_bound();
        }

        func::TypedPredicate<EdgeId> condition = LengthUpperBound<Graph>(curr_env.graph(), length);
        if (args.size() > 2 && (args[2] == "Y" || args[2] == "y")) {
            cout << "Trying to activate genome quality condition" << endl;
            if (curr_env.genome().size() == 0) {
                cout << "No reference was provided!!!" << endl;
            } else {
                cout << "Genome quality condition will be used" << endl;

                curr_env.graph_pack().ClearQuality();
                curr_env.graph_pack().FillQuality();
//                condition = make_shared<make_shared<debruijn_graph::ZeroQualityCondition<Graph, Index>>(curr_env.graph(), edge_qual);
                condition = std::bind(&debruijn_graph::EdgeQuality<Graph>::IsZeroQuality,
                                      std::ref(curr_env.graph_pack().edge_qual), std::placeholders::_1);
            }
        }
        debruijn::simplification::SimplifInfoContainer info(debruijn_graph::config::pipeline_type::base);
        info.set_chunk_cnt(10);
        debruijn::simplification::TipClipperInstance(curr_env.graph(), condition, info, (omnigraph::EdgeRemovalHandlerF<Graph>)nullptr)->Run();
    }
};
}
