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
            string Usage() const {
                string answer;
                answer = answer + "Command `set_max_vertices` \n" + 
                                "Usage:\n" + 
                                "> set_max_vertices <max_vertices> \n" + 
                                " You should specify an integer, which is an upper bound for the number of vertices in the picture.";
                return answer;
            }

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
            string Usage() const {
                string answer;
                answer = answer + "Command `set_folder` \n" + 
                                "Usage:\n" + 
                                "> set_folder <folder_name> \n" + 
                                " You should specify a string, which is a new name for a pictures' folder.";
                return answer;
            }
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
            string Usage() const {
                string answer;
                answer = answer + "Command `set_file_name` \n" + 
                                "Usage:\n" + 
                                "> set_file_name <file_base_name>\n" + 
                                " You should specify a string, which is a new base_name for all the pictures, that you generate.";
                return answer;
            }
        
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

