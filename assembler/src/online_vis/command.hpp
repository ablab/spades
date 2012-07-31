#pragma once

#include "environment.hpp"
#include "command_type.hpp"
#include "loaded_environments.hpp"
#include "argument_list.hpp"
#include "errors.hpp"

namespace online_visualization {
    
    class Command {

        protected:
            CommandType command_id_;

            virtual size_t MinArgNumber() const {
                return 0;
            }

            virtual bool CheckCorrectness(const ArgumentList& arg_list) const {
                return false;
            }

            bool CheckEnoughArguments(const vector<string>& args) const {
                bool result = args.size() >= MinArgNumber();
                if (!result) 
                    FireNotEnoughArguments();
                return result;
            }

        public:
            virtual string Usage() const {
                string answer;
                answer = answer + " Welcome to GAF (Graph Analysis Framework). This framework allows to work with the de Bruijn Graph interactively.\n " +
                                " You can see the list of command names below. To see a command's help message just type\n" +
                                "> help <command_name>\n" +
                                " The list of command names : \n" + 
                                " exit\n" +
                                " help\n" +
                                " load\n" +
                                " list\n" +
                                " switch\n" +
                                " rep\n" +
                                " load_genome\n" +
                                " set_folder\n" +
                                " set_file_name\n" +
                                " set_max_vertices\n" +
                                " fill_pos\n" +
                                " clear_pos\n" +
                                " vertex\n" +
                                " edge\n" +
                                " position\n" +
                                " paths";
                return answer;
            }

            Command(CommandType command_id) : command_id_(command_id) 
            {
            }

            virtual ~Command() {
            }

            CommandType command_id() const {
                return command_id_;
            }

            // system command, curr_env can point to null
            virtual void Execute(EnvironmentPtr& curr_env, LoadedEnvironments& loaded_environments, const ArgumentList& arg_list) const = 0;

            // virtual void Execute(EnvironmentPtr& curr_env, const ArgumentList& arg_list) const = 0;

    };

    class LocalCommand : public Command {
        
        public:
            LocalCommand(CommandType command_id) : Command(command_id)
            {
            }
            
            // command for the current environment
            virtual void Execute(Environment& curr_env, const ArgumentList& arg_list) const = 0;

            // !!!! NO OVERRIDING !!!!
            virtual void Execute(EnvironmentPtr& curr_env, LoadedEnvironments& loaded_environments, const ArgumentList& arg_list) const {
                Execute(*curr_env, arg_list);
            }

    };


}
