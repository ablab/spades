#pragma once

#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"
#include "standard_vis.hpp"
#include "command.hpp"
#include "loaded_environments.hpp"
#include "environment.hpp"
#include "command_struct.hpp"

#include "all_commands.hpp"

namespace online_visualization {


    //TODO : BatchCommand
    class OnlineVisualizer {
    private:

        shared_ptr<Environment> current_environment_;

        LoadedEnvironments loaded_environments;

        void AddAllCommands() {
            AddCommand(shared_ptr<Command>(new NullCommand));
            AddCommand(shared_ptr<Command>(new LoadCommand));
            AddCommand(shared_ptr<Command>(new ExitCommand));
            AddCommand(shared_ptr<Command>(new HelpCommand));
            AddCommand(shared_ptr<Command>(new ListCommand));
            AddCommand(shared_ptr<Command>(new SwitchCommand));
            AddCommand(shared_ptr<Command>(new ReplayCommand));
            AddCommand(shared_ptr<Command>(new LoadGenomeCommand));
            AddCommand(shared_ptr<Command>(new LogCommand));
            AddCommand(shared_ptr<Command>(new SaveBatchCommand));
            AddCommand(shared_ptr<Command>(new BatchCommand));

            AddCommand(shared_ptr<Command>(new SetMaxVertCommand));
            AddCommand(shared_ptr<Command>(new SetFolderCommand));
            AddCommand(shared_ptr<Command>(new SetFileNameCommand));

            AddCommand(shared_ptr<Command>(new FillPositionCommand));
            AddCommand(shared_ptr<Command>(new ClearPositionCommand));

            AddCommand(shared_ptr<Command>(new DrawVertexCommand));
            AddCommand(shared_ptr<Command>(new DrawEdgeCommand));
            AddCommand(shared_ptr<Command>(new DrawPositionCommand));
            AddCommand(shared_ptr<Command>(new DrawPartOfGenomeCommand));
            AddCommand(shared_ptr<Command>(new DrawContigCommand));
            AddCommand(shared_ptr<Command>(new ShowPositionCommand));

            AddCommand(shared_ptr<Command>(new PrintPathsCommand));
            AddCommand(shared_ptr<Command>(new PrintContigsStatsCommand));
        }

    public:
        OnlineVisualizer()
        {
	        fs::path p = fs::path(cfg::get().load_from) / "late_pair_info_counted";
            stringstream ss("load default " + p.string());
            AddAllCommands();
            Command& LoadCommand = GetCommand(CommandId("load"));
            LoadCommand.Execute(current_environment_, loaded_environments, ss);
        }

        void run() {

            vector<string>& history = GetHistory();
            //const size_t max_buffer_size = 10000;

            while (true) {
                cout << "GAF$> ";
                string command_with_args;
                getline(cin, command_with_args);
                stringstream ss(command_with_args);
                ArgumentList arg_list(ss);
                string processed_command = arg_list.Preprocess(history);

                DEBUG("processed string " << processed_command);

                const string& command_string = arg_list.GetAllArguments()[0];
                Command& command = GetCommand(CommandId(command_string));
                command.Execute(current_environment_, loaded_environments, arg_list);

                history.push_back(processed_command);
            }
        }
    private:
        DECL_LOGGER("OnlineVisualizer");
    };

}
