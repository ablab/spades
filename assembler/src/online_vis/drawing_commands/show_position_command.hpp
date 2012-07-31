#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../command_type.hpp"
#include "../errors.hpp"
#include "../argument_list.hpp"

#include "drawing_command.hpp"

namespace online_visualization {
    class ShowPositionCommand : public DrawingCommand {
        private:
            int ShowPicture(Environment& curr_env, runtime_k::RtSeq kmer, string label = "") const {
                kmer = curr_env.kmer_mapper().Substitute(kmer);
                if (!curr_env.index().contains(kmer)) {
                    FireNoCorrespondingGraphLocation(label);
                    return -1;
                }
                pair<EdgeId, size_t> position = curr_env.index().get(kmer);
                if (position.second * 2 < curr_env.graph().length(position.first))
                    return DrawingCommand::ShowPicture(curr_env, curr_env.graph().EdgeStart(position.first), label);
                else
                    return DrawingCommand::ShowPicture(curr_env, curr_env.graph().EdgeEnd(position.first), label);
            }

        protected:
            size_t MinArgNumber() const {
                return 1;   
            }
            
            bool CheckCorrectness(const vector<string>& args) const {
                if (!CheckEnoughArguments(args))
                    return false;
                bool result = true;
                result = result & CheckIsNumber(args[0]);
                return result;
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `show_position` \n" + 
                                "Usage:\n" + 
                                "> show_position <position>\n" + 
                                " This command prints pictures for a neigbourhood of an edge in the DB graph, which corresponds to a given genome position.\n" + 
                                " You should specify an integer position in the genome, which location you want to look at.";
                return answer;
            }

            ShowPositionCommand() : DrawingCommand(CommandType::show_position)
            {
            }

            void Execute(Environment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args_ = arg_list.GetAllArguments();
                if (!CheckCorrectness(args_))
                    return;

                int position = GetInt(args_[0]);
                Sequence genome = curr_env.genome();
                if (arg_list["--rc"] == "true") {
                    cout << "Inverting genome..." << endl;
                    genome = !genome;
                }
                if (CheckPositionBounds(position, genome.size())) {
                    int result = ShowPicture(curr_env, genome.Subseq(position).start<runtime_k::RtSeq::max_size>(cfg::get().K + 1), args_[0]);
                    if (result) 
                        FireGenericError("Something is wrong");
                }

            }
    };
}
