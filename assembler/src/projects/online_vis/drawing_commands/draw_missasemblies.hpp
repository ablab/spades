//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "io/reads/wrapper_collection.hpp"

namespace online_visualization {
class DrawMisassemblies : public DrawingCommand {

protected:
    size_t MinArgNumber() const {
        return 1;
    }

    bool CheckCorrectness(const vector<string>& args) const {
        if (!CheckEnoughArguments(args))
            return false;
        if(!CheckFileExists(args[1]))
            return false;
        return true;
    }

private:

    vector<EdgeId> FilterByLength(Graph& g, const vector<EdgeId>& edges) const {
        vector<EdgeId> filtered_edges;
        for(auto e : edges) {
            if(g.length(e) > 500) {
                filtered_edges.push_back(e);
            }
        }
        return filtered_edges;
    }

    vector<EdgeId> FilterNonUnique(Graph& g, const vector<EdgeId>& edges, const vector<EdgeId>& genome_edges) const {
        vector<EdgeId> filtered_edges;
        std::set<EdgeId> set_edges;
        std::set<EdgeId> non_unique;
        std::set<EdgeId> genome_set;
        std::set<EdgeId> non_unique_genome;


        for(auto e : edges) {
            if(set_edges.find(e) != set_edges.end()) {
                non_unique.insert(e);
            }
            set_edges.insert(e);
        }

        for(auto e : genome_edges) {
            if(genome_set.find(e) != genome_set.end()) {
                non_unique_genome.insert(e);
            }
            genome_set.insert(e);
        }



        for(auto e : edges) {
            if(non_unique.find(e) == non_unique.end() && non_unique_genome.find(e) == non_unique_genome.end()) {
                filtered_edges.push_back(e);
                INFO("Put " << g.int_id(e) << " into filtered set");
            }
        }
        return filtered_edges;
    }


    void ProcessContig(DebruijnEnvironment& curr_env, MappingPath<EdgeId>& genome_path, MappingPath<EdgeId>& reverse_genome_path, MappingPath<EdgeId>& path, string name = "") const {
        genome_path.join(reverse_genome_path);
        vector<EdgeId> genome_edges = curr_env.path_finder().FindReadPath(genome_path);
        vector<EdgeId> rc_genome_edges = curr_env.path_finder().FindReadPath(reverse_genome_path);
        vector<EdgeId> rc_and_usual_genome_edges(genome_edges);
        utils::push_back_all(rc_and_usual_genome_edges, rc_genome_edges);
        vector<EdgeId> edges = path.simple_path();
        auto filtered_edges = FilterNonUnique(curr_env.graph(), edges, rc_and_usual_genome_edges);
        if(filtered_edges.size() < 2)
            return;

        auto it_genome = find(rc_and_usual_genome_edges.begin(), rc_and_usual_genome_edges.end(), filtered_edges[0]);
        size_t index_genome = it_genome - rc_and_usual_genome_edges.begin();
        size_t i = 0;


        auto it_contig = find(edges.begin(), edges.end(), filtered_edges[i]);
        while(it_contig == edges.end()) {
            ++i;
            if(i > filtered_edges.size()) {
                return;
            }
            it_contig = find(edges.begin(), edges.end(), filtered_edges[i]);
        }
        size_t index_contig = it_contig - edges.begin();
        INFO("Now at edge " << curr_env.graph().int_id(filtered_edges[i]));
        const int allowed_error = 3000;
        int real_difference = (int)genome_path[index_genome].second.initial_range.start_pos - (int)path[index_contig].second.initial_range.start_pos;
        INFO("Diff is set to " << real_difference);


        while(i < filtered_edges.size()) {
            INFO("Now at edge " << curr_env.graph().int_id(filtered_edges[i]));
            it_genome = find(rc_and_usual_genome_edges.begin(), rc_and_usual_genome_edges.end(), filtered_edges[i]);
            it_contig = find(edges.begin(), edges.end(), filtered_edges[i]);

            size_t index_genome = it_genome - rc_and_usual_genome_edges.begin();
            size_t index_contig = it_contig - edges.begin();

            if(it_genome == rc_and_usual_genome_edges.end()) {
                vector<EdgeId> path_to_draw;

                while(it_genome == rc_and_usual_genome_edges.end()) {
                    ++i;
                    if(i == filtered_edges.size())
                    {
                        break;
                    }
                    it_genome = find(rc_and_usual_genome_edges.begin(), rc_and_usual_genome_edges.end(), filtered_edges[i]);
                }

                auto new_it_contig = find(edges.begin(), edges.end(), filtered_edges[i]);
                size_t new_index_contig = new_it_contig - edges.begin();

                for(size_t z = index_contig; z <= new_index_contig ; ++z) {
                    path_to_draw.push_back(edges[z]);
                }


                DrawPicturesAlongPath(curr_env, path_to_draw, name + "_" + std::to_string(curr_env.graph().int_id(filtered_edges[i])));
                real_difference = (int)genome_path[index_genome].second.initial_range.start_pos - (int)path[index_contig].second.initial_range.start_pos;
                INFO("Diff is set to " << real_difference);
                continue;
            }

            int difference = (int)genome_path[index_genome].second.initial_range.start_pos - (int)path[index_contig].second.initial_range.start_pos;
            if(abs(difference - real_difference) > allowed_error) {
                real_difference = (int)genome_path[index_genome].second.initial_range.start_pos - (int)path[index_contig].second.initial_range.start_pos;
                vector<EdgeId> path_to_draw;
                path_to_draw.push_back(genome_path[index_genome].first);
                DrawPicturesAlongPath(curr_env, path_to_draw, name + "_" + std::to_string(curr_env.graph().int_id(filtered_edges[i])));
                INFO("Diff is set to " << real_difference);
            }
            ++i;

        }
    }

public:
    DrawMisassemblies() : DrawingCommand("draw_misassemblies") {

    }

    string Usage() const {
        string answer;
        answer = answer + "Command `draw_misassemblies` \n" +
                        "Usage:\n" +
                        "> draw_misassemblies <file with missasembled quast contigs>\n" +
                        "Reference genome should be loaded to use this command.\n" +
                        "This command tries to draw exact places of misassembles.";
        return answer;
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args)) {
            return;
        }

        if(curr_env.genome() == Sequence()) {
            cout << "Reference should be loaded. Command will not be executed" << endl;
            return;
        }

        string file = args[1];
        
        visualization::position_filler::FillPos(curr_env.graph_pack(), file, "miss", true);
        cout << "All contigs are mapped" << endl;


        auto genome_mapping_path = curr_env.mapper().MapSequence(curr_env.genome());
        auto rc_genome_mapping_path = curr_env.mapper().MapSequence(!curr_env.genome());

        cout << "Genome is mapped" << endl;

        io::FileReadStream reader(file);
        while(reader.eof()) {
            io::SingleRead read;
            reader >> read;
            auto mapping_path = curr_env.mapper().MapRead(read);
            ProcessContig(curr_env, genome_mapping_path, rc_genome_mapping_path, mapping_path, read.name());
            cout << "Read " << read.name() << " is processed." << endl;
        }
    }

};
}
