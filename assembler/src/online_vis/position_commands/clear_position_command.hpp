#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../command_type.hpp"
#include "../errors.hpp"

namespace online_visualization {
    class ClearPositionCommand : public LocalCommand {
        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `clear_pos` \n" + 
                                "Usage:\n" + 
                                "> clear_pos\n" + 
                                " This command resets the graph and clears all the labels you previously filled in.\n";
                return answer;
            }

            ClearPositionCommand() : LocalCommand(CommandType::clear_pos) 
            {
            }

            void Execute(Environment& curr_env, const ArgumentList& args) const {
                curr_env.ResetPositions();
            }
        
    };
}
