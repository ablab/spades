//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "../argument_list.hpp"

#include "drawing_command.hpp"

namespace online_visualization {
    class DrawPartOfGenomeCommand : public DrawingCommand {
        private:
            void CheckPathIntegrity(const GraphDistanceFinder<Graph>& dist_finder, EdgeId first_edge, EdgeId second_edge) const {
                const vector<size_t>& distances = dist_finder.GetGraphDistancesLengths(first_edge, second_edge);
                if (distances[0] == 0) {
                    INFO("Edges " << first_edge << " and " << second_edge << " are neighbouring");
                } else 
                    INFO("Edges " << first_edge << " and " << second_edge << " are at distance of " << distances[0]);
            }

        private:
            void DrawPicturesAlongGenomePart(DebruijnEnvironment& curr_env, const Sequence& piece_of_genome, string label = "") const {
                const MappingPath<EdgeId>& mapping_path = curr_env.mapper().MapSequence(piece_of_genome);
                DrawingCommand::DrawPicturesAlongPath(curr_env, mapping_path, label);
            }
            
            void CountStatsAlongGenomePart(DebruijnEnvironment& curr_env, Sequence& piece_of_genome, string label = "") const {    
                GraphDistanceFinder<Graph> dist_finder(curr_env.graph(), cfg::get().ds.IS(), cfg::get().ds.RL(),
                        size_t(cfg::get().ds.is_var()));
                cout << "Statistics for the part of genome :" << endl;
                const MappingPath<EdgeId>& mapping_path = curr_env.mapper().MapSequence(piece_of_genome);
                for (size_t i = 0; i < mapping_path.size(); ++i) {
                    cout << "Edge # " << i << endl;
                    //const pair<EdgeId, MappingRange>& mapping_edge = mapping_path[i];
                    
                    if (i > 0) {
                        INFO("Checking connection between neighbouring edges");
                        CheckPathIntegrity(dist_finder, mapping_path[i - 1].first, mapping_path[i].first);
                    }
                }
            }

        protected:
            size_t MinArgNumber() const {
                return 2;   
            }
            
            bool CheckCorrectness(const vector<string>& args) const {
                if (!CheckEnoughArguments(args))
                    return false;
                if (!CheckIsNumber(args[1]))
                    return false;
                if (!CheckIsNumber(args[2]))
                    return false;

                return true;
            }
 
        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `draw_part_of_genome` \n" + 
                                " Usage:\n" + 
                                "> genome <first_pos> <second_pos> [--rc] [-r]\n" + 
                                " Prints a .dot picture of a substring [first_pos, second_pos] of the genome.\n" +
                                " Optionally you can use a flag -r, whether you want the tool to invert the positions,\n" +
                                " and an option --rc, if you would like to see the pictures of the second strand.";
                return answer;
            }

            DrawPartOfGenomeCommand() : DrawingCommand("draw_part_of_genome")
            {
            }

            void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();
                if (!CheckCorrectness(args))
                    return;

                size_t first_position = (size_t) GetInt(args[1]);
                size_t second_position = (size_t) GetInt(args[2]);
                Sequence genome = curr_env.genome();
                if (arg_list["rc"] == "true") {
                    cout << "Inverting genome..." << endl;
                    genome = !genome;
                }

                //experimental
                if (arg_list.contains("r")) {
                    cout << "Inverting positions..." << endl;
                    first_position = genome.size() - curr_env.k_value() - 1 - first_position;
                    second_position = genome.size() - curr_env.k_value() - 1 - first_position;
                }

                if (CheckPositionBounds(first_position, genome.size(), curr_env.k_value()) && 
                        CheckPositionBounds(second_position, genome.size(), curr_env.k_value())) 
                {
                    const Sequence& part_of_genome = genome.Subseq(first_position, second_position);
                    string label = args[1] + "_" + args[2];
                    DrawPicturesAlongGenomePart(curr_env, part_of_genome, label);
                }

            }
    };
}
