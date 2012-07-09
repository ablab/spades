#pragma once
#include "environment.hpp"
#include "command_struct.hpp"


namespace online_visualization {
    
    typedef map<string, shared_ptr<Environment>> LoadedEnvironments;

    //typedef typename online_visualization::CommandType CommandType;
    
    class Command {
        CommandType command_id_;
    public:
        Command(CommandType command_id) : command_id_(command_id) {
        }

        CommandType command_id() const {
            return command_id_;
        }

        virtual void Execute(LoadedEnvironments& loaded_environments
                , shared_ptr<Environment>& curr_env, string args) const {
            Execute(*curr_env, args);
        }

        virtual void Execute(Environment& curr_env, string args) const = 0;

        virtual ~Command() {
        }
    };

    Command& GetCommand(CommandType command_id);

    // loading new environment from folder with saves
    class LoadCommand : public Command {
        private:
            //TODO: IO check
            shared_ptr<Environment> MakeNewEnvironment(string saves, string name) const {
                fs::path path = fs::path(saves);
                
                Environment environment(name, path.string());
                return shared_ptr<Environment>(&environment);
            }

        public:
            LoadCommand() : 
                Command(CommandType::load)
            {
            }

            void Execute(LoadedEnvironments& loaded_environments, shared_ptr<Environment>& curr_env, string args) const {
                stringstream ss(args);
                string name; 
                ss >> name;
                string saves;
                ss >> saves;

                cout << name << " " << saves << endl;

                for (auto iterator = loaded_environments.begin(); iterator != loaded_environments.end(); ++iterator) {
                    if (name == iterator->first) {
                        cerr << "Name " << name << " already exists" << endl;
                        return;
                    }
                }

                const shared_ptr<Environment>& new_env = MakeNewEnvironment(saves, name);
                loaded_environments.insert(make_pair(name, new_env));
                curr_env = new_env;
            }

            void Execute(Environment& curr_env, string args) const {
                cout << "NOTHING TO SAY HERE" << endl;
                exit(1);
            }
    };

    // exit
    class ExitCommand : public Command {
        public:
            ExitCommand() : 
                Command(CommandType::_exit_)
            {
            }

            void Execute(Environment& curr_env, string args) const {
                cout << "exiting" << endl;
                exit(0);
            }
    };


    typedef map<CommandType, shared_ptr<Command>> CommandMapping;   
 
    void AddMapping(CommandMapping& mapping, Command* command) {
        mapping.insert(make_pair(command->command_id(), shared_ptr<Command>(command)));
    }

    CommandMapping FillCommandMapping() {
        CommandMapping mapping;
        AddMapping(mapping, new LoadCommand);
        AddMapping(mapping, new ExitCommand);
        return mapping;
    }

    Command& GetCommand(CommandType command_id) {
        static CommandMapping mapping = FillCommandMapping();
        auto it = mapping.find(command_id);
        VERIFY(it != mapping.end());
        return *(it->second);
    }
}
