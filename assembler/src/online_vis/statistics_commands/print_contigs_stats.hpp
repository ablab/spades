#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../command_type.hpp"
#include "../errors.hpp"

namespace online_visualization {
    class PrintContigsStatsCommand : public LocalCommand {
        private:
            void CheckPathIntegrity(EdgeId first_edge, EdgeId second_edge) const {
                            
            }

            void CheckCoverage(const MappingRange& first_range, const MappingRange& second_range) const {
                   
            }

        private:
            //void CountStatsAlongGenomePart(Environment& curr_env, Sequence& piece_of_genome, string label = "") const {    
                //cout << "Statistics for the part of genome :" << endl;
                //const MappingPath<EdgeId>& mapping_path = curr_env.mapper().MapSequence(piece_of_genome);
                //for (size_t i = 0; i < mapping_path.size(); ++i) {
                    //cout << "Edge # " << i << endl;
                    //const pair<EdgeId, MappingRange>& mapping_edge = mapping_path[i];
                    
                    //if (i > 0) {
                        //LOG("Checking connection between neighbouring edges");
                        //CheckPathIntegrity(mapping_path[i - 1].first, mapping_path[i].first);
                        //CheckCoverage(mapping_path[i - 1].second, mapping_path[i].second);
                    //}
                //}
            //}

            void ProcessContig(Environment& curr_env, const Sequence& contig, const MappingPath<EdgeId>& genome_path, const string& contig_name) const {
                LOG("Checking the contig " << contig_name);
                LOG("Length " << contig.size());
                const MappingPath<EdgeId>& contig_path = curr_env.mapper().MapSequence(contig);
                bool found = false;
                EdgeId first_edge = contig_path[0].first;
                for (size_t i = 0; i < genome_path.size(); ++i) {
                    TRACE("i-th edge of the genome " << genome_path[i].first);
                    if (genome_path[i].first == first_edge) {
                        found = true;
                        for (size_t j = 1; j < contig_path.size(); ++j) {
                            if (genome_path[i + j].first != contig_path[j].first) {
                                LOG("Break in the edge " << contig_path[j]);
                                return;
                            }
                        }
                    }
                }
                if (!found) {
                    LOG("First edge was not found");   
                } else {
                    LOG("No misassemblies");   
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
                answer = answer + "Command `draw_part_of_genome` \n" + 
                                "Usage:\n" + 
                                "> position <position> [--rc] [-r]\n" + 
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
                        ProcessContig(curr_env, contig, genome_path, "CONTIG_" + read.name());
                        ProcessContig(curr_env, !contig, genome_path, "CONTIG_" + read.name() + "_RC");
                    }
                }
            }
        private:
            DECL_LOGGER("PrintContigsStats");
    };
}
