#pragma once
#include "environment.hpp"
#include "command_struct.hpp"


namespace online_visualization {
    
    typedef map<string, shared_ptr<Environment>> LoadedEnvironments;

    //typedef typename online_visualization::CommandType CommandType;
    
    class Command {
        CommandType command_id_;
        size_t args_max_number;
    public:
        Command(CommandType command_id) : command_id_(command_id) {
        }

        CommandType command_id() const {
            return command_id_;
        }

        virtual void Execute(LoadedEnvironments& loaded_environments
                , shared_ptr<Environment>& curr_env, string args) {
            Execute(*curr_env, args);
        }

        virtual void Execute(Environment& curr_env, string args) = 0;

        virtual ~Command() {
        }
    };

    Command& GetCommand(CommandType command_id);

    namespace command_impl {

        typedef map<CommandType, shared_ptr<Command>> CommandMapping;   
     
        void AddMapping(CommandMapping& mapping, Command* command) {
            mapping.insert(make_pair(command->command_id(), shared_ptr<Command>(command)));
        }

        CommandMapping FillCommandMapping() {
            CommandMapping mapping;
            AddMapping(mapping, new LoadCommand());
            return mapping;
        }
    
        
        Command& GetCommand(CommandType command_id) {
            static CommandMapping mapping = command_impl::FillCommandMapping();
            auto it = mapping.find(command_id);
            VERIFY(it != mapping.end());
            return *(it->second);
        }
    }


    class LoadCommand : public Command {
        private:
            shared_ptr<Environment> MakeNewEnvironment(string args) {
                
            }

        public:
            LoadCommand() : 
                Command(CommandType::load)
            {
            }

            void Execute(LoadedEnvironments& loaded_environment, shared_ptr<Environment>& curr_env, string args) {
                MakeNewEnvironment(args); 
            }

            void Execute(Environment& curr_env, string args) {
                INFO("NOTHING TO SAY HERE");
            }
    }
}
