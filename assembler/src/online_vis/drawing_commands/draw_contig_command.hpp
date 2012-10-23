#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../command_type.hpp"
#include "../errors.hpp"

namespace online_visualization {
    class DrawContigCommand : public DrawingCommand {

        private:
            void DrawPicturesAlongGenomePart(Environment& curr_env, const Sequence& piece_of_genome, string label = "") const {
                const MappingPath<EdgeId>& mapping_path = curr_env.mapper().MapSequence(piece_of_genome);
                DrawingCommand::DrawPicturesAlongPath(curr_env, mapping_path, label);
            }
            
        protected:
            size_t MinArgNumber() const {
                return 2;   
            }
            
            bool CheckCorrectness(const vector<string>& args) const {
                if (!CheckEnoughArguments(args))
                    return false;

                return true;
            }
 
        public:
            string Usage() const {
                string answer;
                answer = answer + "Command `draw_contig` \n" + 
                                "Usage:\n" + 
                                "> contig <name_of_contig> <contigs_file>\n" + 
                                " Draws graph pictures for a contig.";
                return answer;
            }

            DrawContigCommand() : DrawingCommand(CommandType::draw_contig)
            {
            }

            void Execute(Environment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();
                if (!CheckCorrectness(args))
                    return;

                string contig_name = args[1];
                string contigs_file = args[2];
                if (!CheckFileExists(contigs_file))
                    return;

                io::Reader irs(contigs_file);

                while (!irs.eof()) {
                    io::SingleRead read;
                    irs >> read;
                    LOG("Contig " << read.name() << " is being processed now");

                    // if read is valid and also the name contains a given string <contig_name> as a substring.
                    if (read.IsValid() && read.name().find(contig_name) != string::npos) {
                        const Sequence& contig = read.sequence();
                        const string& label = read.name();
                        DrawPicturesAlongGenomePart(curr_env, contig, label);
                        LOG("Contig " << read.name() << " has been drawn");
                    }                        
                }

            }
          
    };
}
