#pragma once

#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"
#include "standard_vis.hpp"
#include "command.hpp"
#include "environment.hpp"
#include "command_struct.hpp"

namespace online_visualization {


    class OnlineVisualizer {
    private:
        shared_ptr<Environment> current_environment_;

    public:
        OnlineVisualizer()
        {
	        fs::path p = fs::path(cfg::get().load_from) / "constructed_graph";
            //stringstream ss;
            //ss << "default ";
            //ss << p.string();
            //cout << ss.str() << endl;;
            stringstream ss("default " + p.string());
            Command& LoadCommand = GetCommand(CommandId("load"));
            LoadCommand.Execute(current_environment_, ss);
        }

        void run() {

            while (true) {
                cout << "> ";
                string command_with_args;
                getline(cin, command_with_args);
                stringstream ss(command_with_args);
                string command_string;
                ss >> command_string;

                Command& command = GetCommand(CommandId(command_string));
                command.Execute(current_environment_, ss);
                

                //} else if (com == "paths") {
                    //FindPaths(ss);
            }
        }


        //void SetMaxVertices(stringstream & ss) {
            //size_t max;
            //ss >> max;
            //max_vertices_ = max;
        //}


        //void SetFolder(stringstream & ss) {
            //ss >> folder_;
        //}

        //void SetFileName(stringstream & ss) {
            //ss >> file_name_base_;
        //}

    };

}
