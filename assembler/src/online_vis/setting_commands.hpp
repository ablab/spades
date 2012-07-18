#pragma once

#include "environment.hpp"
#include "command.hpp"
#include "command_type.hpp"
#include "vis_utils.hpp"

namespace online_visualization {
    class SetMaxVertCommand : public Command {
        protected:
            size_t MinArgNumber() const {
                return 1;   
            }
            
            bool CheckCorrectness(vector<string> args) const {
                if (args.size() < MinArgNumber()) {
                    cout << "Not enough arguments" << endl;
                    cout << "Please try again" << endl;
                    return false;
                }
                return IsNumber(args[0]);
            }

        public:
            SetMaxVertCommand() : Command(CommandType::set_max_vertices)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                if (!CheckCorrectness(args_)) {
                    return;
                }
                size_t max_v = GetInt(args_[0]);
                curr_env->set_max_vertices(max_v);
            }
    };

    class SetFolderCommand : public Command {
        protected:
            size_t MinArgNumber() const {
                return 1;   
            }
            
            bool CheckCorrectness(const vector<string>& args) const {
                if (args.size() < MinArgNumber()) {
                    cout << "Not enough arguments" << endl;
                    cout << "Please try again" << endl;
                    return false;
                }
                return true;
            }

        public:
            SetFolderCommand() : Command(CommandType::set_folder)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                if (!CheckCorrectness(args_)) {
                    return;
                }
                string folder_name = args_[0];
                curr_env->set_folder(folder_name);
            }
    };

    class SetFileNameCommand : public Command {
        protected:
            size_t MinArgNumber() const {
                return 1;   
            }

            bool CheckCorrectness(const vector<string>& args) const {
                if (args.size() < MinArgNumber()) {
                    cout << "Not enough arguments" << endl;
                    cout << "Please try again" << endl;
                    return false;
                }
                return true;
            }

        public:
            SetFileNameCommand() : Command(CommandType::set_file_name)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                if (!CheckCorrectness(args_)) {
                    cout << "Please try again" << endl;
                    return;
                }
                string file_name = args_[0];
                curr_env->set_file_name(file_name);
            }
    };
}

