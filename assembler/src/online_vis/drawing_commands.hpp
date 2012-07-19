#pragma once

#include "environment.hpp"
#include "command.hpp"
#include "command_type.hpp"
#include "vis_utils.hpp"

namespace online_visualization {

    class DrawingCommand : public Command {
        protected:
            void DrawPicture(EnvironmentPtr& curr_env, VertexId v) const {
                make_dir(curr_env->folder_);
                stringstream namestream;
                namestream << curr_env->folder_ << "/" << curr_env->file_name_base_ << "_" << curr_env->picture_counter_ << ".dot";
                string file_name = namestream.str();
                stringstream linkstream;
                linkstream  << curr_env->folder_ << "/" << curr_env->file_name_base_ << "_latest.dot";
                VertexNeighborhoodFinder<Graph> splitter(curr_env->graph(), v, curr_env->max_vertices_, curr_env->edge_length_bound_);
                //EdgePosGraphLabeler<Graph> labeler(curr_env->graph(), gp_.edge_pos);
                WriteComponents <Graph> (curr_env->graph(), splitter, file_name, *DefaultColorer(curr_env->graph(), curr_env->coloring_), curr_env->tot_lab_);
                WriteComponents <Graph> (curr_env->graph(), splitter, linkstream.str(), *DefaultColorer(curr_env->graph(), curr_env->coloring_), curr_env->tot_lab_);
                cout << "Picture is written to " << file_name << endl;
                curr_env->picture_counter_++;
            }

            void DrawVertex(EnvironmentPtr& curr_env, size_t vertex_id) const {
                DrawPicture(curr_env, curr_env->int_ids().ReturnVertexId(vertex_id));
            }


        public:
            DrawingCommand(CommandType command_type) : Command(command_type)
            {
            }

            virtual ~DrawingCommand()
            {
            }
    };

    class DrawPositionCommand : public DrawingCommand {
        private:
            void DrawPicture(EnvironmentPtr curr_env, Seq<debruijn_graph::K + 1> kmer) const {
                kmer = curr_env->kmer_mapper().Substitute(kmer);
                if (!curr_env->index().contains(kmer)) {
                    cout << "No corresponding graph location " << endl;
                    return;
                }
                pair<EdgeId, size_t> position = curr_env->index().get(kmer);
                if (position.second * 2 < curr_env->graph().length(position.first))
                    DrawingCommand::DrawPicture(curr_env, curr_env->graph().EdgeStart(position.first));
                else
                    DrawingCommand::DrawPicture(curr_env, curr_env->graph().EdgeEnd(position.first));
            }

            bool CheckPositionBounds(int position, size_t total_size) const {
                bool result = !(position + debruijn_graph::K + 1 > total_size);
                if (!result) {
                    cout << "Ignoring request. Position is out of range : required position is " 
                         << position << " while length of the sequence is "
                         << total_size << endl;
                    cout << "Please try again" << endl;
                }
                return result;
            }

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

                bool result = IsNumber(args[0]);
                if (!result) {
                    cout << "The argument " << args[0] << " is not a number" << endl;
                    cout << "Please try again" << endl;
                }
                return result;
            }

        public:
            DrawPositionCommand() : DrawingCommand(CommandType::draw_position)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                if (!CheckCorrectness(args_))
                    return;

                int position = GetInt(args_[0]);
                Sequence genome = curr_env->genome();
                if (args_.size() > 1) {
                    if (args_[1] == "--rc") {
                        genome = !genome;
                    }
                }
                if (CheckPositionBounds(position, genome.size()))
                    DrawPicture(curr_env, genome.Subseq(position).start<debruijn_graph::K + 1>());

            }
    };

    class DrawVertexCommand : public DrawingCommand {
        private:
            bool CheckVertexExists(const IdTrackHandler<Graph>& int_ids, size_t vertex_id) const {
                VertexId vertex = int_ids.ReturnVertexId(vertex_id);
                if ((vertex == VertexId(NULL))) {
                    cout << "Ignoring request. Vertex " << vertex_id << " does not exist" << endl;
                    cout << "Please try again" << endl;
                    cout << vertex << endl;
                    return false;
                }
                return true;
            }

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
                bool result = IsNumber(args[0]);
                if (!result) {
                    cout << "The argument " << args[0] << " is not a number" << endl;
                    cout << "Please try again" << endl;   
                }
                return result;
            }

        public:
            DrawVertexCommand() : DrawingCommand(CommandType::draw_vertex)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);

                if (!CheckCorrectness(args_))
                    return;
                
                size_t vertex_id = GetInt(args_[0]);
                if (CheckVertexExists(curr_env->int_ids(), vertex_id)) 
                    DrawVertex(curr_env, vertex_id);
            }
    };

    class DrawEdgeCommand : public DrawingCommand {
        private:
            bool CheckEdgeExists(const IdTrackHandler<Graph>& int_ids, size_t edge_id) const {
                EdgeId edge = int_ids.ReturnEdgeId(edge_id);
                if ((edge == EdgeId(NULL))) {
                    cout << "Ignoring request. Edge " << edge_id << " does not exist" << endl;
                    cout << "Please try again" << endl;
                    return false;
                }
                return true;
            }

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
                bool result = IsNumber(args[0]);
                if (!result) {
                    cout << "The argument " << args[0] << " is not a number" << endl;
                    cout << "Please try again" << endl;   
                }
                return result;
            }

            void DrawEdge(EnvironmentPtr& curr_env, EdgeId edge) const {
                DrawingCommand::DrawPicture(curr_env, curr_env->graph().EdgeStart(edge)); 
            }

            void DrawEdge(EnvironmentPtr& curr_env, size_t edge_id) const {
                DrawEdge(curr_env, curr_env->int_ids().ReturnEdgeId(edge_id));   
            }

        public:
            DrawEdgeCommand() : DrawingCommand(CommandType::draw_edge)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);

                if (!CheckCorrectness(args_))
                     return;

                size_t edge_id = GetInt(args_[0]);
                if (CheckEdgeExists(curr_env->int_ids(), edge_id)) {
                    DrawEdge(curr_env, edge_id);
                }
            }
    };
}

