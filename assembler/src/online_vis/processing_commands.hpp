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
        omnigraph::ClipTips(curr_env.graph(), length);
    }
};
}
