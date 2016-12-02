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
#include "../argument_list.hpp"

#include "drawing_command.hpp"

namespace online_visualization {
    class DrawPositionCommand : public DrawingCommand {
        private:
            void DrawPicture(DebruijnEnvironment& curr_env, RtSeq kmer, string label = "") const {
                kmer = curr_env.kmer_mapper().Substitute(kmer);
                if (!curr_env.index().contains(kmer)) {
                    cout << "No corresponding graph location " << endl;
                    return;
                }
                pair<EdgeId, size_t> position = curr_env.index().get(kmer);
                if (position.second * 2 < curr_env.graph().length(position.first))
                    DrawingCommand::DrawPicture(curr_env, curr_env.graph().EdgeStart(position.first), label);
                else
                    DrawingCommand::DrawPicture(curr_env, curr_env.graph().EdgeEnd(position.first), label);
            }

        protected:
            size_t MinArgNumber() const {
                return 1;   
            }
            
            bool CheckCorrectness(const vector<string>& args) const {
                if (!CheckEnoughArguments(args))
                    return false;
                if (!CheckIsNumber(args[1]))
                    return false;

                return true;
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `draw_position` \n" + 
                                "Usage:\n" + 
                                "> position <position> [--rc] [-r]\n" + 
                                " You should specify an integer position in the genome, which location you want to look at. Optionally you can use a flag -r, whether you want the tool to invert the positions,\n" +
                                "and an option --rc, if you would like to see the pictures of the second strand.";
                return answer;
            }

            DrawPositionCommand() : DrawingCommand("draw_position")
            {
            }

            void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();
                if (!CheckCorrectness(args))
                    return;

                int position = GetInt(args[1]);
                Sequence genome = curr_env.genome();
                if (arg_list["rc"] == "true") {
                    cout << "Inverting genome...";
                    genome = !genome;
                }

                if (CheckPositionBounds(position, genome.size(), curr_env.k_value())) {
                    DrawPicture(curr_env, genome.Subseq(position).start<RtSeq>(curr_env.k_value() + 1), args[1]);
                }

            }
    };
}
