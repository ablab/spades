#pragma once

#include "environment.hpp"
#include "command.hpp"
#include "command_type.hpp"
#include "vis_utils.hpp"

namespace online_visualization {

    class DrawingCommand : public Command {
        protected:
            void DrawPicture(EnvironmentPtr& curr_env, VertexId vertex, string label = "") const {
                make_dir(curr_env->folder_);
                stringstream namestream;
                namestream << curr_env->folder_ << "/" << curr_env->file_name_base_ << "_" << label << "_" << curr_env->picture_counter_ << ".dot";
                string file_name = namestream.str();
                //stringstream linkstream;
                //linkstream  << curr_env->folder_ << "/" << curr_env->file_name_base_ << "_latest.dot";
                VertexNeighborhoodFinder<Graph> splitter(curr_env->graph(), vertex, curr_env->max_vertices_, curr_env->edge_length_bound_);
                //EdgePosGraphLabeler<Graph> labeler(curr_env->graph(), gp_.edge_pos);
                WriteComponents<Graph>(curr_env->graph(), splitter, file_name, *DefaultColorer(curr_env->graph(), curr_env->coloring_), curr_env->tot_lab_);
                //WriteComponents <Graph> (curr_env->graph(), splitter, linkstream.str(), *DefaultColorer(curr_env->graph(), curr_env->coloring_), curr_env->tot_lab_);
                cout << "Picture is written to " << file_name << endl;
                
                curr_env->picture_counter_++;
            }

              
            int ShowPicture(EnvironmentPtr& curr_env, VertexId vertex, string label = "") const {
                DrawPicture(curr_env, vertex, label);
                stringstream command_line_string;
                command_line_string << "gnome-open " << curr_env->folder_ << "/" << curr_env->file_name_base_ << "_" 
                                    << label << "_" << curr_env->picture_counter_ 
                                    << "_*_.dot & > /dev/null < /dev/null";
                int result = system(command_line_string.str().c_str());

                return result;
            }

            void DrawVertex(EnvironmentPtr& curr_env, size_t vertex_id, string label = "") const {
                DrawPicture(curr_env, curr_env->int_ids().ReturnVertexId(vertex_id), label);
            }


        public:
            DrawingCommand(CommandType command_type) : Command(command_type)
            {
            }

            virtual ~DrawingCommand()
            {
            }
    };

    class DrawVertexCommand : public DrawingCommand {
        private:
            bool CheckVertexExists(const IdTrackHandler<Graph>& int_ids, size_t vertex_id) const {
                VertexId vertex = int_ids.ReturnVertexId(vertex_id);
                if ((vertex == VertexId(NULL))) {
                    cout << "Ignoring the request. Vertex " << vertex_id << " does not exist" << endl;
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
            string Usage() const {
                string answer;
                answer = answer + "Command `draw_vertex` \n" + 
                                "Usage:\n" + 
                                "vertex <vertex_id>\n" + 
                                "This command prints pictures for a neigbourhood of a vertex in the DB graph.\n" + 
                                "You should specify an id of the vertex in the DB graph, which neighbourhood you want to look at.";
                return answer;
            }
            
            DrawVertexCommand() : DrawingCommand(CommandType::draw_vertex)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);

                if (!CheckCorrectness(args_))
                    return;
                
                size_t vertex_id = GetInt(args_[0]);
                if (CheckVertexExists(curr_env->int_ids(), vertex_id)) 
                    DrawVertex(curr_env, vertex_id, args_[0]);
            }
    };

    class DrawEdgeCommand : public DrawingCommand {
        private:
            bool CheckEdgeExists(const IdTrackHandler<Graph>& int_ids, size_t edge_id) const {
                EdgeId edge = int_ids.ReturnEdgeId(edge_id);
                if ((edge == EdgeId(NULL))) {
                    cout << "Ignoring the request. Edge " << edge_id << " does not exist" << endl;
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

            void DrawEdge(EnvironmentPtr& curr_env, EdgeId edge, string label = "") const {
                DrawingCommand::DrawPicture(curr_env, curr_env->graph().EdgeStart(edge), label); 
            }

            void DrawEdge(EnvironmentPtr& curr_env, size_t edge_id, string label = "") const {
                DrawEdge(curr_env, curr_env->int_ids().ReturnEdgeId(edge_id), label);
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `draw_edge` \n" + 
                                "Usage:\n" + 
                                "edge <edge_id>\n" + 
                                "This command prints pictures for a neigbourhood of an edge in the DB graph.\n" + 
                                "You should specify an id of the edge in the DB graph, which location you want to look at.";
                return answer;
            }

            DrawEdgeCommand() : DrawingCommand(CommandType::draw_edge)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);

                if (!CheckCorrectness(args_))
                     return;

                size_t edge_id = GetInt(args_[0]);
                if (CheckEdgeExists(curr_env->int_ids(), edge_id)) {
                    DrawEdge(curr_env, edge_id, args_[0]);
                }
            }
    };

    class ShowPositionCommand : public DrawingCommand {
        private:
            int ShowPicture(EnvironmentPtr curr_env, runtime_k::RtSeq kmer, string label = "") const {
                kmer = curr_env->kmer_mapper().Substitute(kmer);
                if (!curr_env->index().contains(kmer)) {
                    cout << "No corresponding graph location " << endl;
                    return -1;
                }
                pair<EdgeId, size_t> position = curr_env->index().get(kmer);
                if (position.second * 2 < curr_env->graph().length(position.first))
                    return DrawingCommand::ShowPicture(curr_env, curr_env->graph().EdgeStart(position.first), label);
                else
                    return DrawingCommand::ShowPicture(curr_env, curr_env->graph().EdgeEnd(position.first), label);
            }

            bool CheckPositionBounds(int position, size_t total_size) const {
                bool result = !(position + cfg::get().K + 1 > total_size);
                if (!result) {
                    cout << "Ignoring the request. Position is out of range : required position is " 
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
            string Usage() const {
                string answer;
                answer = answer + "Command `show_position` \n" + 
                                "Usage:\n" + 
                                "show_position <position>\n" + 
                                "This command prints pictures for a neigbourhood of an edge in the DB graph, which corresponds to a given genome position.\n" + 
                                "You should specify an integer position in the genome, which location you want to look at.";
                return answer;
            }

            ShowPositionCommand() : DrawingCommand(CommandType::show_position)
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
                if (CheckPositionBounds(position, genome.size())) {
                    int result = ShowPicture(curr_env, genome.Subseq(position).start<runtime_k::RtSeq::max_size>(cfg::get().K + 1), args_[0]);
                    if (result) {
                        cout << "Something went wrong" << endl;
                        cout << "Please try again" << endl;
                    }
                }

            }
    };

    //TODO: Finish!
    class DrawPartOfGenomeCommand : public DrawingCommand {
        private:
            void DrawPicture(EnvironmentPtr curr_env, runtime_k::RtSeq kmer, string label = "") const {
                kmer = curr_env->kmer_mapper().Substitute(kmer);
                if (!curr_env->index().contains(kmer)) {
                    cout << "No corresponding graph location " << endl;
                    return;
                }
                pair<EdgeId, size_t> position = curr_env->index().get(kmer);
                if (position.second * 2 < curr_env->graph().length(position.first))
                    DrawingCommand::DrawPicture(curr_env, curr_env->graph().EdgeStart(position.first), label);
                else
                    DrawingCommand::DrawPicture(curr_env, curr_env->graph().EdgeEnd(position.first), label);
            }

            bool CheckPositionBounds(int position, size_t total_size) const {
                bool result = !(position + cfg::get().K + 1 > total_size);
                if (!result) {
                    cout << "Ignoring the request. Position is out of range : required position is " 
                         << position << " while length of the sequence is "
                         << total_size << endl;
                    cout << "Please try again" << endl;
                }
                return result;
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

                for (size_t i = 0; i < 2; ++i) {
                    if (!IsNumber(args[i])) {
                        cout << "The argument " << args[i] << " is not a number" << endl;
                        cout << "Please try again" << endl;
                        return false;
                    }
                }
                return true;
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `draw_part_of_genome` \n" + 
                                "Usage:\n" + 
                                "position <position>\n" + 
                                "You should specify an integer position in the genome, which location you want to look at.";
                return answer;
            }

            DrawPartOfGenomeCommand() : DrawingCommand(CommandType::draw_part_of_genome)
            {
            }

            void Execute(EnvironmentPtr& curr_env, stringstream& args) const {
                const vector<string>& args_ = SplitInTokens(args);
                if (!CheckCorrectness(args_))
                    return;

                int first_position = GetInt(args_[0]);
                int second_position = GetInt(args_[1]);
                Sequence genome = curr_env->genome();
                if (args_.size() > 2) {
                    if (args_[2] == "--rc") {
                        genome = !genome;
                    }
                }
                if (CheckPositionBounds(first_position, genome.size()) && CheckPositionBounds(second_position, genome.size()))
                    DrawPicture(curr_env, genome.Subseq(first_position, second_position).start<runtime_k::RtSeq::max_size>(cfg::get().K + 1), args_[0]);

            }
    };

    class DrawPositionCommand : public DrawingCommand {
        private:
            void DrawPicture(EnvironmentPtr curr_env, runtime_k::RtSeq kmer, string label = "") const {
                kmer = curr_env->kmer_mapper().Substitute(kmer);
                if (!curr_env->index().contains(kmer)) {
                    cout << "No corresponding graph location " << endl;
                    return;
                }
                pair<EdgeId, size_t> position = curr_env->index().get(kmer);
                if (position.second * 2 < curr_env->graph().length(position.first))
                    DrawingCommand::DrawPicture(curr_env, curr_env->graph().EdgeStart(position.first), label);
                else
                    DrawingCommand::DrawPicture(curr_env, curr_env->graph().EdgeEnd(position.first), label);
            }

            bool CheckPositionBounds(int position, size_t total_size) const {
                bool result = !(position + cfg::get().K + 1 > total_size);
                if (!result) {
                    cout << "Ignoring the request. Position is out of range : required position is " 
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

                if (!IsNumber(args[0])) {
                    cout << "The argument " << args[0] << " is not a number" << endl;
                    cout << "Please try again" << endl;
                    return false;
                }
                return true;
            }

        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `draw_position` \n" + 
                                "Usage:\n" + 
                                "position <position>\n" + 
                                "You should specify an integer position in the genome, which location you want to look at.";
                return answer;
            }

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
                if (CheckPositionBounds(position, genome.size())) {
                    DrawPicture(curr_env, genome.Subseq(position).start<runtime_k::RtSeq::max_size>(cfg::get().K + 1), args_[0]);
                }

            }
    };
}

