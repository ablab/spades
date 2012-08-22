#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../command_type.hpp"
#include "../errors.hpp"

namespace online_visualization {

    class FillPositionCommand : public LocalCommand {

        protected:
            size_t MinArgNumber() const {
                return 2;   
            }

            bool CheckCorrectness(const vector<string>& args) const {
                if (!CheckEnoughArguments(args))
                    return false;
                //const string& name = args[1];
                const string& file = args[2];

                bool result = true;

                result &= CheckFileExists(file);

                return result;
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

            FillPositionCommand() : LocalCommand(CommandType::fill_pos)
            {
            }

            void Execute(Environment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();
                if (!CheckCorrectness(args))
                    return;

                string name = args[1];
                string file = args[2];

                io::Reader irs(file);
                
                FillerClass& filler = curr_env.filler();
                while (!irs.eof()) {
                    io::SingleRead read;
                    irs >> read;
                    if (read.IsValid()) {
                        Sequence contig = read.sequence();
                        filler.Process(contig,  name + "_" + read.name());
                        filler.Process(!contig, name + "_" + read.name() + "_RC");
                    }
                }
            }
    };
}
