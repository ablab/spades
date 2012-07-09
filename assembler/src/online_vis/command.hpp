#pragma once

#include "environment.hpp"
#include "command_type.hpp"
#include "loaded_environments.hpp"


namespace online_visualization {
    
    typedef shared_ptr<Environment> EnvironmentPtr;

    typedef map<string, EnvironmentPt > LoadedEnvironments;

    class Command {

        protected:
            CommandType command_id_;
            LoadedEnvironments& environments;
            
            virtual void PostError() const {
                
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
