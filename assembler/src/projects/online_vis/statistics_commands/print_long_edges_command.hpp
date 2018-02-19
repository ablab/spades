#pragma once

#include "../debruijn_environment.hpp"
#include "../command.hpp"

namespace online_visualization {

class PrintLongEdgesCommand : public LocalCommand<DebruijnEnvironment> {

 protected:
    size_t MinArgNumber() const {
        return 3;
    }

    bool CheckCorrectness(const vector<string>& args) const {
        return CheckEnoughArguments(args);
    }

 public:
    string Usage() const {
        return string() + "Command `print_long_edges` \n" + " Usage:\n"
            + "> print_long_edges <contig_name> <path_to_contigs> <min length of the edge> \n"
            + " Prints list of long edges in the contig path";
    }

    PrintLongEdgesCommand()
        : LocalCommand<DebruijnEnvironment>("print_long_edges") {
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        string contig_name = args[1];
        LOG("Printing edges in contigs " << contig_name);

        bool starts_with = false;
        if (contig_name[contig_name.size() - 1] == '*') {
            starts_with = true;
            contig_name = contig_name.substr(0, contig_name.size() - 1);
        }
        string contigs_file = args[2];
        if (!CheckFileExists(contigs_file))
            return;
        unsigned long length_threshold = std::stoul(args[3]);

        io::FileReadStream reader(contigs_file);

        while (!reader.eof()) {
            io::SingleRead read;
            reader >> read;
            //LOG("Contig " << read.name() << " is being processed now");

            // if the name contains a given string <contig_name> as a substring.
            if((starts_with && read.name().find(contig_name) != string::npos) || contig_name == read.name()) {
                auto mapping_path = curr_env.mapper().MapRead(read).simple_path();
                set<EdgeId> long_edges;
                std::copy_if(mapping_path.begin(), mapping_path.end(), std::inserter(long_edges, long_edges.begin()),
                             [&curr_env, length_threshold](const EdgeId& edge) {
                               return curr_env.graph().length(edge) >= length_threshold;
                             });
                cout << long_edges.size() << " long edges" << endl;
                auto add_edge_length = [&curr_env](size_t current_sum, const EdgeId &edge) {
                  return current_sum + curr_env.graph().length(edge);
                };

                size_t long_edge_length = std::accumulate(long_edges.begin(), long_edges.end(), 0, add_edge_length);

                size_t total_length = std::accumulate(mapping_path.begin(), mapping_path.end(), 0, add_edge_length);
                cout << "Total length: " << total_length << endl;
                cout << "Long edge total length: " << long_edge_length << endl;

                for (const auto& edge: long_edges) {
                    cout << edge.int_id() << std::endl;
                }
            }
        }
    }

 private:

    DECL_LOGGER("JunctionSequenceCommand");
};


}
