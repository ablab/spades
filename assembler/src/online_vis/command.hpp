#pragma once

#include "environment.hpp"
#include "command_type.hpp"
#include "loaded_environments.hpp"


namespace online_visualization {
    
    typedef shared_ptr<Environment> EnvironmentPtr;

    typedef map<string, EnvironmentPtr > LoadedEnvironments;

    class Command {

        protected:
            CommandType command_id_;
            LoadedEnvironments& environments;

            //ArgumentList args;
            
            bool IsNumber(const string& s) const {
                 if (s.empty())
                     return false;
                 for  (auto iter = s.begin(); iter != s.end(); ++iter) {
                    if (!std::isdigit(*iter))
                        return false;
                 }
                 return true;
            }

            virtual size_t MinArgNumber() const {
                return 0;
            }

            // TODO: create ArgumentList class, which can return the value for an option you ask.
            //const set<string> GetAllTokens(stringstream& args) const { 
                //set<string> answer;
                //while (!args.eof()) {
                    //string arg;
                    //args >> arg;
                    //answer.insert(arg);
                //}
                //return answer;
            //}

            //const vector<string> SortArgsList(set<string> args) const {
                //return {};
            //}

            const vector<string> SplitInTokens(stringstream& args) const { 
                vector<string> answer;
                while (!args.eof()) {
                    string arg;
                    args >> arg;
                    answer.push_back(arg);
                }
                return answer;
            }

            int GetInt(string str) const {
                stringstream ss(str);
                int ans;
                ss >> ans;
                return ans;
            }
            
            virtual bool CheckCorrectness(vector<string> args) const {
                return false;
            }

        public:
            virtual string Usage() const {
                string answer;
                answer = answer + "Welcome to GAF (Graph Analysis Framework). This framework allows to work with the de Bruijn Graph interactively.\n " +
                                "You can see the list of command names below. To see a command's help message just type\n" +
                                "> help <command_name>\n" +
                                "The list of command names\n" + 
                                "exit\n" +
                                "help\n" +
                                "load\n" +
                                "list\n" +
                                "switch\n" +
                                "rep\n" +
                                "set_folder\n" +
                                "set_file_name\n" +
                                "set_max_vertices\n" +
                                "fill_pos\n" +
                                "clear_pos\n" +
                                "vertex\n" +
                                "edge\n" +
                                "position\n" +
                                "paths";
                return answer;
            }

            Command(CommandType command_id) : command_id_(command_id), environments(GetLoadedEnvironments()) 
            {
            }

            virtual ~Command() {
            }

            CommandType command_id() const {
                return command_id_;
            }

            virtual void Execute(EnvironmentPtr& curr_env, stringstream& args) const = 0;

    };

}
