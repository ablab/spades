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

namespace online_visualization {
    class ClearPositionCommand : public LocalCommand<DebruijnEnvironment> {
        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `clear_pos` \n" +
                                "Usage:\n" +
                                "> clear_pos\n" +
                                " This command resets the graph and clears all the labels you previously filled in.\n";
                return answer;
            }

            ClearPositionCommand() : LocalCommand<DebruijnEnvironment>("clear_pos")
            {
            }

            void Execute(DebruijnEnvironment& curr_env, const ArgumentList&) const {
                curr_env.ResetPositions();
            }

    };
}
