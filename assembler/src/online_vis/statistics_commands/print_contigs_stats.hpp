#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../command_type.hpp"
#include "../errors.hpp"
#include "omni/omni_utils.hpp"

namespace online_visualization {
    class PrintContigsStatsCommand : public LocalCommand {
        //typedef vector<EdgeId> Path;
        
        private:
            vector<EdgeId> TryCloseGap(const Graph& graph, VertexId v1, VertexId v2) const {
                if (v1 == v2)
                    return vector<EdgeId>();
                TRACE("Trying to close gap between v1 =" << graph.int_id(v1) << " and v2 =" << graph.int_id(v2));
                PathStorageCallback<Graph> path_storage(graph);
            
                //  todo reduce value after investigation
                PathProcessor<Graph> path_processor(graph, 0, 50, v1, v2, path_storage);
                path_processor.Process();

                if (path_storage.size() == 0) {
                    TRACE("Failed to find closing path");
                    return vector<EdgeId>();
                } else if (path_storage.size() == 1) {
                    TRACE("Unique closing path found");
                } else {
                    TRACE("Several closing paths found, first chosen");
                }
                vector<EdgeId> answer = path_storage.paths().front();
                TRACE("Gap closed");
                TRACE("Cumulative closure length is " 
                        << CummulativeLength(graph, answer));
                return answer;
            }

            vector<EdgeId> TryFixPath(Environment& curr_env, const vector<EdgeId>& edges) const {
                vector<EdgeId> answer;
                if (edges.empty()) {
                    //  WARN("Mapping path was empty");
                    return vector<EdgeId>();
                }
                //  VERIFY(edges.size() > 0);
                answer.push_back(edges[0]);
                for (size_t i = 1; i < edges.size(); ++i) {
                    vector<EdgeId> closure = TryCloseGap(curr_env.graph(), curr_env.graph().EdgeEnd(edges[i - 1]), curr_env.graph().EdgeStart(edges[i]));
                    answer.insert(answer.end(), closure.begin(), closure.end());
                    answer.push_back(edges[i]);
                }
                return answer;
            }

            Path<EdgeId> TryFixPath(Environment& curr_env, const Path<EdgeId>& path) const {
                return Path<EdgeId>(TryFixPath(curr_env, path.sequence()), path.start_pos(), path.end_pos());
            }

        private:

            bool ProcessContig(Environment& curr_env, const Sequence& contig, const MappingPath<EdgeId>& genome_path, const string& contig_name) const {
                cout << " Checking the contig " << contig_name << endl;
                cout << " Length " << contig.size() << endl;
                const Path<EdgeId>& genome_path_completed = TryFixPath(curr_env, genome_path.simple_path());
                const MappingPath<EdgeId>& contig_path = curr_env.mapper().MapSequence(contig);
                bool found = false;
                EdgeId first_edge = contig_path[0].first;
                for (size_t i = 0; i < genome_path_completed.size(); ++i) {
                    TRACE("i-th edge of the genome " << genome_path_completed[i]);
                    if (genome_path_completed[i] == first_edge) {
                        found = true;
                        for (size_t j = 1; j < contig_path.size(); ++j) {
                            if (genome_path_completed[i + j] != contig_path[j].first) {
                                cout << " Break in the edge " << contig_path[j] << endl;
                                return false;
                            }
                        }
                    }
                }
                if (!found) {
                    cout << " First edge was not found" << endl;   
                    return false;
                } else {
                    cout << " No misassemblies" << endl;
                    return true;
                }
            }

        protected:
            size_t MinArgNumber() const {
                return 1;   
            }
            
            bool CheckCorrectness(const vector<string>& args) const {
                if (!CheckEnoughArguments(args))
                    return false;

                const string& file = args[0];
                if (!CheckFileExists(file))
                    return false;

                return true;
            }
 
        public:
            //TODO : REDO
            string Usage() const {
                string answer;
                answer = answer + "Command `print_contigs_stats` \n" + 
                                "Usage:\n" + 
                                "> print_contigs_stats <position> [--rc] [-r]\n" + 
                                " You should specify an integer position in the genome, which location you want to look at. Optionally you can use a flag -r, whether you want the tool to invert the positions,\n" +
                                "and an option --rc, if you would like to see the pictures of the second strand.";
                return answer;
            }

            PrintContigsStatsCommand() : LocalCommand(CommandType::print_contigs_stats)
            {
            }

            void Execute(Environment& curr_env, const ArgumentList& args) const {
                const vector<string>& args_ = args.GetAllArguments();
                if (!CheckCorrectness(args_))
                    return;

                string file = args_[0];
                if (!CheckCorrectness(args_))
                    return;

                io::Reader irs(file);
                
                const Sequence& genome = curr_env.genome();

                const MappingPath<EdgeId>& genome_path = curr_env.mapper().MapSequence(genome);

                while (!irs.eof()) {
                    io::SingleRead read;
                    irs >> read;
                    if (read.IsValid()) {
                        const Sequence& contig = read.sequence();
                        bool result = true;
                        result = result & ProcessContig(curr_env, contig, genome_path, "CONTIG_" + read.name());
                        result = result & ProcessContig(curr_env, !contig, genome_path, "CONTIG_" + read.name() + "_RC");
                        if (result) {
                            cout << " contig " << read.name() << " is okay" << endl;   
                        }
                        else 
                            cout << " contig " << read.name() << " is MISASSEMBLED" << endl;
                    }
                }
            }
        private:
            DECL_LOGGER("PrintContigsStats");
    };
}
