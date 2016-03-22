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
    class DrawConnectedCommand : public DrawingCommand {


        protected:
            size_t MinArgNumber() const {
                return 0;
            }

            bool CheckCorrectness(const vector<string>& args) const {
                if (!CheckEnoughArguments(args))
                    return false;
                if (args.size() > 1 && !CheckIsNumber(args[1]))
                    return false;
                if (args.size() > 2 && !CheckIsNumber(args[2]))
                    return false;
                return true;
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `draw_connected` \n" +
                                "Usage: draw_connected\n" +
                                "or draw_connected min_component_size \n" +
                                "or draw_connected min_component_size max_component_size\n" +
                                "Takes no arguments, draw connected components with given size of the graph";
                return answer;
            }

            DrawConnectedCommand() : DrawingCommand("draw_connected")
            {
            }

            void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();
                if (!CheckCorrectness(args))
                    return;
                int min_size = 1;
                int max_size = 1000000000;
                if (args.size() > 1) min_size = GetInt(args[1]);
                if (args.size() > 2) max_size = GetInt(args[2]);
                DrawConnectedComponents(curr_env, min_size, max_size);
            }
    };
}
