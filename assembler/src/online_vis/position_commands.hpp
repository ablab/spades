#pragma once

#include "environment.hpp"
#include "command.hpp"
#include "command_type.hpp"
#include "vis_utils.hpp"

namespace online_visualization {

    class FillPositionCommand : public Command {

        protected:
            size_t MinArgNumber() const {
                return 2;   
            }

            bool CheckCorrectness(const vector<string>& args) const {
                if (args.size() < MinArgNumber()) {
                    cout << "Not enough arguments" << endl;
                    cout << "Please try again" << endl;
                    return false;
                }
                const string& name = args[0];
                const string& file = args[1];

                bool correct = true;

                correct &= CheckFileExists(file);
                for (auto iterator = environments.begin(); iterator != environments.end(); ++iterator) {
                    if (name == iterator->first) {
                        cout << "Name " << name << " already exists" << endl;
                        return false;
                    }
                }
                return correct;
            }
        public:
            FillPositionCommand() : Command(CommandType::fill_pos)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                string name = args_[0];
                string file = args_[1];
                if (!CheckCorrectness(args_)) {
                    cout << "Please try again" << endl;
                    return;
                }

                io::Reader irs(file);
                
                FillerClass& filler = curr_env->filler();
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

    class ClearPositionCommand : public Command {

        public:
            ClearPositionCommand() : Command(CommandType::clear_pos) 
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                curr_env->ResetPositions();
            }
        
    };


}
