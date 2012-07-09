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
                

                //else if (com == "set_folder") {
                    //SetFolder(ss);
                //} else if (com == "set_file_name") {
                    //SetFileName(ss);
                //} else if (com == "set_max_vertices") {
                    //SetMaxVertices(ss);
                //} else if (com == "clear_pos") {
                    //ClearPos(ss);
                //} else if (com == "vertex") {
                    //DrawVertex(ss);
                //} else if (com == "edge") {
                    //DrawEdge(ss);
                //} else if (com == "position") {
                    //DrawGenomePosition(ss);
                //} else if (com == "paths") {
                    //FindPaths(ss);
                //} else {
                    //cout << "ignoring command " << command_with_args << endl;
                //}
            }
        }

        //void FindPaths(stringstream & ss) {
            //size_t from;
            //size_t to;
            //size_t max_length;
            //ss >> from;
            //ss >> to;
            //ss >> max_length;
            //CountingCallback callback;
            //PathProcessor<Graph> pp(gp_.g, 0, max_length,
                    //gp_.g.int_ids().ReturnVertexId(from),
                    //gp_.g.int_ids().ReturnVertexId(to), callback);
            //pp.Process();
            //cout << callback.get_cnt() << " paths found" << endl;
        //}

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
