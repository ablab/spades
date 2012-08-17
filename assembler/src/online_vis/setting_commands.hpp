#pragma once

#include "environment.hpp"
#include "command.hpp"
#include "command_type.hpp"
#include "errors.hpp"

namespace online_visualization {
    class SetMaxVertCommand : public LocalCommand {
        protected:
            size_t MinArgNumber() const {
                return 1;   
            }
            
            bool CheckCorrectness(const vector<string>& args) const {
                if (!CheckEnoughArguments(args))
                    return false;
                return CheckIsNumber(args[1]);
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

            SetMaxVertCommand() : LocalCommand(CommandType::set_max_vertices)
            {
            }

            void Execute(Environment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();
                if (!CheckCorrectness(args)) {
                    return;
                }
                size_t max_v = GetInt(args[1]);
                curr_env.set_max_vertices(max_v);
            }
    };

    class SetFolderCommand : public LocalCommand {
        protected:
            size_t MinArgNumber() const {
                return 1;   
            }
            
            bool CheckCorrectness(const vector<string>& args) const {
                return CheckEnoughArguments(args);
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
            SetFolderCommand() : LocalCommand(CommandType::set_folder)
            {
            }

            void Execute(Environment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();
                if (!CheckCorrectness(args)) 
                    return;
                string folder_name = args[1];
                curr_env.set_folder(folder_name);
            }
    };

    class SetFileNameCommand : public LocalCommand {
        protected:
            size_t MinArgNumber() const {
                return 1;   
            }

            bool CheckCorrectness(const vector<string>& args) const {
                return CheckEnoughArguments(args);
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
        
            SetFileNameCommand() : LocalCommand(CommandType::set_file_name)
            {
            }

            void Execute(Environment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();
                if (!CheckCorrectness(args))
                    return;
                string file_name = args[1];
                curr_env.set_file_name(file_name);
            }
    };
}

