#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "io/wrapper_collection.hpp"

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
        	if(non_unique.find(e) == non_unique.end() && non_unique_genome.find(e) != non_unique_genome.end()) {
        		filtered_edges.push_back(e);
        	}
        }
        return filtered_edges;
    }


    void ProcessContig(DebruijnEnvironment& curr_env, MappingPath<EdgeId>& genome_path, MappingPath<EdgeId>& path, string name = "") const {
        vector<EdgeId> genome_edges = genome_path.simple_path();
        vector<EdgeId> edges = path.simple_path();
        auto filtered_edges = FilterNonUnique(curr_env.graph(), edges, genome_edges);
        if(filtered_edges.size() < 2)
            return;

        auto it_genome = find(genome_edges.begin(), genome_edges.end(), filtered_edges[0]);
        size_t index_genome = it_genome - genome_edges.begin();
        size_t i = 0;


        auto it_contig = find(edges.begin(), edges.end(), filtered_edges[i]);
        while(it_contig != edges.end()) {
        	++i;
        	if(i > filtered_edges.size()) {
        		return;
        	}
        	it_contig = find(edges.begin(), edges.end(), filtered_edges[i]);
        }
        size_t index_contig = it_contig - edges.begin();

        const int allowed_error = 3000;
        int real_difference = (int)genome_path[index_genome].second.initial_range.start_pos - (int)path[index_contig].second.initial_range.start_pos;


        while(i < filtered_edges.size() && it_genome != genome_edges.end()) {
            it_genome = find(genome_edges.begin(), genome_edges.end(), filtered_edges[i]);
            it_contig = find(edges.begin(), edges.end(), filtered_edges[i]);

            size_t index_genome = it_genome - genome_edges.begin();
            size_t index_contig = it_contig - edges.begin();
            if(index_genome == genome_edges.size() || index_contig == edges.size()) {
            	++i;
                continue;
            }
            int difference = (int)genome_path[index_genome].second.initial_range.start_pos - (int)path[index_contig].second.initial_range.start_pos;
            if(abs(difference - real_difference) > allowed_error) {
                DrawVertex(curr_env, curr_env.graph().EdgeStart(filtered_edges[i]).int_id(), name);
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
        auto reader = make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(file));
        FillerClass& filler = curr_env.filler();
        while (!reader->eof()) {
            io::SingleRead read;
            (*reader) >> read;
            Sequence contig = read.sequence();
            filler.Process(contig,  "miss_" + read.name());
            filler.Process(!contig, "miss_" + read.name() + "_RC");
        }
        reader->close();
        cout << "All contigs are mapped" << endl;
        reader = make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(file));

        auto genome_mapping_path = curr_env.mapper().MapSequence(curr_env.genome());
        cout << "Genome is mapped" << endl;

        while(!reader->eof()) {
            io::SingleRead read;
            (*reader) >> read;
            Sequence contig = read.sequence();
            cout << "Read " << read.name() << " is processed." << endl;

            auto mapping_path = curr_env.mapper().MapSequence(contig);
            ProcessContig(curr_env, genome_mapping_path, mapping_path);
        }
    }

};
}
