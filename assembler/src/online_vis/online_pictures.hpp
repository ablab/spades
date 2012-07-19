#pragma once

#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"
#include "standard_vis.hpp"
#include "command.hpp"
#include "environment.hpp"
#include "command_struct.hpp"

#include "system_commands.hpp"
#include "drawing_commands.hpp"
#include "position_commands.hpp"
#include "setting_commands.hpp"

namespace online_visualization {


    class OnlineVisualizer {
    private:
        shared_ptr<Environment> current_environment_;

        void AddAllCommands() {
            AddCommand(shared_ptr<Command>(new NullCommand));
            AddCommand(shared_ptr<Command>(new LoadCommand));
            AddCommand(shared_ptr<Command>(new ExitCommand));
            AddCommand(shared_ptr<Command>(new ListCommand));
            AddCommand(shared_ptr<Command>(new SwitchCommand));
            AddCommand(shared_ptr<Command>(new ReplayCommand));

            AddCommand(shared_ptr<Command>(new SetMaxVertCommand));
            AddCommand(shared_ptr<Command>(new SetFolderCommand));
            AddCommand(shared_ptr<Command>(new SetFileNameCommand));

            AddCommand(shared_ptr<Command>(new FillPositionCommand));
            AddCommand(shared_ptr<Command>(new ClearPositionCommand));

            AddCommand(shared_ptr<Command>(new DrawVertexCommand));
            AddCommand(shared_ptr<Command>(new DrawEdgeCommand));
            AddCommand(shared_ptr<Command>(new DrawPositionCommand));

            AddCommand(shared_ptr<Command>(new PrintPathsCommand));
        }

    public:
        OnlineVisualizer()
        {
	        fs::path p = fs::path(cfg::get().load_from) / "constructed_graph";
            stringstream ss("default " + p.string());
            AddAllCommands();
            Command& LoadCommand = GetCommand(CommandId("load"));
            LoadCommand.Execute(current_environment_, ss);
        }

        void run() {

            vector<string>& history = GetHistory();
            const size_t max_buffer_size = 10000;

            while (true) {
                cout << "> ";
                string command_with_args;
                getline(cin, command_with_args);
                while (history.size() >= max_buffer_size) 
                    history.pop_back();

                history.push_back(command_with_args);
                stringstream ss(command_with_args);
                string command_string;
                ss >> command_string;

                Command& command = GetCommand(CommandId(command_string));
                command.Execute(current_environment_, ss);
            }
        }
    };

}
