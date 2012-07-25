#pragma once

#include "command_struct.hpp"
#include "environment.hpp"
#include "command.hpp"
#include "command_type.hpp"
#include "vis_utils.hpp"

namespace online_visualization {


    typedef online_visualization::Command Command;

    // null
    class NullCommand : public Command {
        public:
            NullCommand() : 
                Command(CommandType::_null_)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
            }
    };

    class HelpCommand : public Command {
        protected:
            bool CheckCorrectness(const vector<string>& args) const {
                 return true;
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "The command `help` allows you to see a help message for any command. \n " +
                                "Usage: \n" +
                                "> help <name_of_command> \n" +
                                " Running `help` without parameters yields a list of all commands.";
                return answer;
            }

            HelpCommand() : 
                Command(CommandType::help)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                if (args_.size() == 0) 
                    cout << Command::Usage() << endl;
                else {
                    if (!CheckCorrectness(args_))
                        return;
                    string command_name = args_[0];
                    Command& command = GetCommand(CommandId(command_name));
                    cout << command.Usage() << endl;
                }
            }
    };

    // exit
    class ExitCommand : public Command {
        public:
            string Usage() const {
                return "The command `exit` allows you to exit this application.";
            }

            ExitCommand() : 
                Command(CommandType::_exit_)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                cout << "Exitting" << endl;
                exit(0);
            }
    };
    
    // loading new environment from folder with saves
    class LoadCommand : public Command {
        private:
            EnvironmentPtr MakeNewEnvironment(string name, string saves) const {
                EnvironmentPtr EnvPointer(new Environment(name, saves));
                return EnvPointer;
            }


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
                const string& saves = args[1];
                for (auto iterator = environments.begin(); iterator != environments.end(); ++iterator) {
                    if (name == iterator->first) {
                        cout << "Name " << name << " already exists" << endl;
                        cout << "Maybe you want to switch to this environment? " << name << endl;
                        cout << "Please try again" << endl;
                        return false;
                    }
                }
                bool correct = true;

                correct &= CheckFileExists(saves + ".grp");
                correct &= CheckFileExists(saves + ".sqn");

                return correct;
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `load` \n" + 
                                "Usage:\n" + 
                                "> load <environment_name> <path_to_saves>\n" + 
                                " You should specify the name of the new environment as well as a path to the graph saves. For example:\n" +
                                "> load GraphSimplified data/saves/simplified_graph\n" + 
                                " would load a new environment with the name `GraphSimplified` from the files\n" + 
                                " in the folder `data/saves/` with the basename `simplified_graph` (simplified_graph.grp, simplified_graph.sqn, e.t.c).";
                return answer;
            }

            LoadCommand() : 
                Command(CommandType::load)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                
                if (!CheckCorrectness(args_))
                    return;

                string name = args_[0]; 
                string saves = args_[1];

                cout << "Loading " << name << " " << saves << endl;
                
                const EnvironmentPtr& new_env = MakeNewEnvironment(name, saves);
                environments.insert(make_pair(name, new_env));
                curr_env = new_env;
            }

    };

    // loading new environment from folder with saves
    class SwitchCommand : public Command {
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
                answer = answer + "Command `switch` \n" + 
                                "Usage:\n" + 
                                " switch <environment_name>\n" + 
                                " You should specify the name of the environment you want to switch to. For example:\n" +
                                "> switch GraphSimplified \n" + 
                                " would switch you to the environment with the name `GraphSimplified`.\n" + 
                                " Of course this environment must be loaded first. To see all loaded environments, run command `list`.";
                return answer;
            }

            SwitchCommand() : 
                Command(CommandType::switch_env)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                
                if (!CheckCorrectness(args_))
                    return;

                string name = args_[0]; 

                bool okay = false;
                for (auto iterator = environments.begin(); iterator != environments.end(); ++iterator) {
                    if (name == iterator->first) {
                        okay = true;
                        curr_env = iterator->second;
                        break;
                    }
                }
                if (!okay) {
                    cout << "Name " << name << " does not exist" << endl;
                    cout << "Please try again" << endl;
                    return;
                } else
                    cout << "Switching to " << name << endl;
            }

    };

    class LoadGenomeCommand : public Command {

        protected:
            size_t MinArgNumber() const {
                return 1;    
            }

            bool CheckCorrectness(const vector<string>& args) const {
                const string& file = args[0];
                if (!CheckFileExists(file))
                    return false;
                return true;
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `load_genome` \n" + 
                                "Usage:\n" + 
                                "> load_genome <path_to_genome>\n" + 
                                " You should specify a path to the genome you want to load from. For example:\n" +
                                " Previously loaded genomes would be lost.";
                return answer;
            }

            LoadGenomeCommand() : Command(CommandType::load_genome)
            {                   
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
	            const vector<string>& args_ = SplitInTokens(args);
                if (!CheckCorrectness(args_))
                    return;
                const string& file = args_[0];
                if (file == "") {
                    cout << "Warning: loading empty genome" << endl;
                    return;
                }
                io::Reader genome_stream(file);
                io::SingleRead genome;
                genome_stream >> genome;
                if (genome.IsValid()) {
                    curr_env->LoadNewGenome(genome.sequence());
                }
                else {
                    cout << "Reference genome (" + file + ") has non-ACGT characters. Skipping it" << endl;
                    cout << "Please try again" << endl;
                }
            }
    };


    class ListCommand : public Command {
        protected:
            bool CheckCorrectness() const {
                return true;   
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `list` \n" +
                                "Usage:\n" +
                                "> list\n" + 
                                " This command lists all loaded environments.";
                return answer;
            }

            ListCommand() : Command(CommandType::list)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                cout << "Environments :" << endl;
                for (auto iter = environments.begin(); iter != environments.end(); ++iter) {
                    cout << iter->first << endl;
                }
                if (curr_env) 
                    cout << "Current environment is " << curr_env->str() << endl;
                else 
                    cout << "Current environment was not set" << endl;
            }
    };

    class ReplayCommand : public Command {
        private:
            bool CheckCorrectness(const vector<string>& args) const {
                if (args.size() == 0)
                    return true;
                return IsNumber(args[0]);
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `replay` \n" + 
                                "Usage:\n" + 
                                "> rep <command_number>\n" + 
                                " Runs the command <command_number> commands before. For example:\n" +
                                "> rep 1 \n" + 
                                " would run the previous command.\n" + 
                                " It is still under development.";
                return answer;
            }
            
            ReplayCommand() : Command(CommandType::replay)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                
                if (!CheckCorrectness(args_))
                    return;

                size_t number = GetInt(args_[0]);

                const vector<string>& history = GetHistory();
                
                cout << "Executing the command " << number << " command(s) before... " << endl;
                const string& command_with_args = history[history.size() - 1 - number];
                cout << command_with_args << endl;

                stringstream ss(command_with_args);
                string command_string;
                ss >> command_string;

                Command& command = GetCommand(CommandId(command_string));
                command.Execute(curr_env, ss);
                
            }
    };


//**********************

    class PrintPathsCommand : public Command {
        
        typedef vector<EdgeId> Path;
        
        private:
            bool CheckVertexExists(const IdTrackHandler<Graph>& int_ids, size_t vertex_id) const {
                VertexId vertex = int_ids.ReturnVertexId(vertex_id);
                if ((vertex == VertexId(NULL))) {
                    cout << "Ignoring request. Vertex " << vertex_id << " does not exist" << endl;
                    cout << "Please try again" << endl;
                    return false;
                }
                return true;
            }

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
                return true;
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `print_paths` \n" + 
                                "Usage:\n" + 
                                "> paths <vertex_from> <vertex_to> [<max_length>] \n" + 
                                " This command prints all paths between two given vertices, that do not exceed `max_length` parameter.\n" +
                                " You should specify two integers (id of vertices), between which you want to find paths. Optionally you can provide `max_length` integer, \n" +
                                " so that tool does not consider paths longer than `max_length`.";
                return answer;
            }
            
            PrintPathsCommand() : Command(CommandType::print_paths)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                if (!CheckCorrectness(args_))
                    return;
                
                
                size_t from = GetInt(args_[0]);
                size_t to = GetInt(args_[1]);
                size_t max_length = 0;
                if (args_.size() > 2) 
                    max_length = GetInt(args_[2]);

                if (!CheckVertexExists(curr_env->int_ids(), from) || !CheckVertexExists(curr_env->int_ids(), to))
                    return;

                PathStorageCallback<Graph> callback(curr_env->graph());
                PathProcessor<Graph> pp(curr_env->graph(), 0, max_length,
                        curr_env->int_ids().ReturnVertexId(from),
                        curr_env->int_ids().ReturnVertexId(to), callback);
                pp.Process();
                const vector<Path>& paths = callback.paths();
                cout << paths.size() << " path(s) have been found : " << endl;
                for (size_t i = 0; i < paths.size(); ++i) {
                    cout << (i + 1) << "-th path  :::  ";
                    for (size_t j = 0; j < paths[i].size(); ++j) {
                        cout << curr_env->int_ids().ReturnIntId(paths[i][j]) << "(" << curr_env->graph().length(paths[i][j]) << ") ";       
                    }
                    cout << endl;
                }
                
            }
    };
        
}
