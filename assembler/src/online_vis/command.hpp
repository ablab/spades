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
