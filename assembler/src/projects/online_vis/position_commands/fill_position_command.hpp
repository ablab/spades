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

    class FillPositionCommand : public LocalCommand<DebruijnEnvironment> {

        protected:
            size_t MinArgNumber() const {
                return 2;   
            }

            bool CheckCorrectness(const vector<string>& args) const {
                if (!CheckEnoughArguments(args))
                    return false;
                //const string& name = args[1];
                const string& file = args[2];
                return CheckFileExists(file);
            }
        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `fill_pos` \n" + 
                                "Usage:\n" + 
                                "> fill_pos <label> <path_to_contigs>\n" + 
                                " This command maps contigs you provide to the graph.\n" + 
                                " You should specify a label of this contigs, which you want to see at the edge in the DB graph.";
                return answer;
            }

            FillPositionCommand() : LocalCommand<DebruijnEnvironment>("fill_pos")
            {
            }

            void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();
                if (!CheckCorrectness(args))
                    return;

                string name = args[1];
                string file = args[2];

                visualization::position_filler::FillPos(curr_env.graph_pack(), file, name, true);
            }
    };
}
