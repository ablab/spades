//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "drawing_commands/drawing_command.hpp"
#include "drawing_commands/draw_position_command.hpp"
#include "drawing_commands/show_position_command.hpp"
#include "drawing_commands/draw_part_of_genome_command.hpp"
#include "drawing_commands/draw_contig_command.hpp"
#include "drawing_commands/draw_poorly_assembled.hpp"
#include "drawing_commands/draw_connected_command.hpp"
#include "drawing_commands/draw_missasemblies.hpp"
#include "environment.hpp"
#include "command.hpp"
#include "errors.hpp"
#include "argument_list.hpp"
#include "drawing_commands/draw_polymorphic_regions.hpp"

namespace online_visualization {

    class DrawVertexCommand : public DrawingCommand {
        protected:
            size_t MinArgNumber() const {
                return 1;   
            }
            
            bool CheckCorrectness(const vector<string>& args) const {
                if (!CheckEnoughArguments(args))
                    return false;
                if (!CheckIsNumber(args[1]))
                    return false;

                return true;
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `draw_vertex` \n" + 
                                "Usage:\n" + 
                                "> draw_vertex <vertex_id>\n" +
                                " This command prints pictures for a neigbourhood of a vertex in the DB graph.\n" + 
                                " You should specify an id of the vertex in the DB graph, which neighbourhood you want to look at.";
                return answer;
            }
            
            DrawVertexCommand() : DrawingCommand("draw_vertex")
            {
            }

            void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();

                if (!CheckCorrectness(args))
                    return;
                size_t vertex_id = GetInt(args[1]);
                if (CheckVertexExists(curr_env.finder(), vertex_id))
                    DrawVertex(curr_env, vertex_id, args[1]);
            }
    };

    class DrawEdgeCommand : public DrawingCommand {
        protected:
            size_t MinArgNumber() const {
                return 1;   
            }
            
            bool CheckCorrectness(const vector<string>& args) const {
                return CheckEnoughArguments(args);
            }

            void DrawEdge(DebruijnEnvironment& curr_env, EdgeId edge, string label = "") const {
                DrawingCommand::DrawPicture(curr_env, curr_env.graph().EdgeStart(edge), label); 
            }

            void DrawEdge(DebruijnEnvironment& curr_env, size_t edge_id, string label = "") const {
                DrawEdge(curr_env, curr_env.finder().ReturnEdgeId(edge_id), label);
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `draw_edge` \n" + 
                                "Usage:\n" + 
                                "> draw_edge <edge_id>\n" +
                                " This command prints pictures for a neigbourhood of an edge in the DB graph.\n" + 
                                " You should specify an id of the edge in the DB graph, which location you want to look at.";
                return answer;
            }

            DrawEdgeCommand() : DrawingCommand("draw_edge")
            {
            }

            void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();

                if (!CheckCorrectness(args))
                     return;

                size_t edge_id = GetInt(args[1]);
                if (CheckEdgeExists(curr_env.finder(), edge_id)) {
                    DrawEdge(curr_env, edge_id, args[1]);
                }
            }
    };
}

