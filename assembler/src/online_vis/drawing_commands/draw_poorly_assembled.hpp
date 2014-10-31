//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "sequence_mapper.hpp"

namespace online_visualization {
    class DrawPoorlyAssembledCommand : public DrawingCommand {
        const double WELL_ASSEMBLED_CONSTANT = 0.7;
        private:
            void DrawPicturesAlongGenomePart(DebruijnEnvironment& curr_env,
                                             const Sequence& piece_of_genome, string label = "") const {
                const MappingPath<EdgeId>& mapping_path = curr_env.mapper().MapSequence(piece_of_genome);
                DrawingCommand::DrawPicturesAlongPath(curr_env, mapping_path, label);
            }

            void DrawContig(DebruijnEnvironment& curr_env, io::SingleRead contig) const {
                Sequence seq = contig.sequence();
                string label = contig.name();
                DrawPicturesAlongGenomePart(curr_env, seq, label);
                LOG("Contig " << contig.name() << " has been drawn");
            }

            io::SingleRead MakeValid(const io::SingleRead& contig) const {
                std::string str = contig.GetSequenceString();
                for (size_t i = 0; i < str.length(); ++i) {
                    if (str[i] == 'N')
                        str[i] = nucl(char(i % 4));
                }
                return io::SingleRead(contig.name(), str);
            }

            bool IsPoorlyAssembled(const GraphPack& gp, io::SingleRead contig, string base_assembly_prefix) const {
                MappingPath<EdgeId> mapping_path = debruijn_graph::MapperInstance(gp)->MapRead(contig);
                auto pos_handler = gp.edge_pos;
                map<string, size_t> base_ctg_2_len;
                for (EdgeId e : mapping_path.simple_path()) {
                    auto positions = pos_handler.GetEdgePositions(e);
                    for (EdgePosition pos : positions) {
                        if (boost::starts_with(pos.contigId, base_assembly_prefix)) {
                            base_ctg_2_len[pos.contigId] += pos.mr.mapped_range.size();
                        }
                    }
                }
                for (pair<string, size_t> entry : base_ctg_2_len)
                    if (double(entry.second) > double(contig.size()) * WELL_ASSEMBLED_CONSTANT)
                        return false;

                return true;
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
                answer = answer + "Command `draw_poorly_assembled` \n" +
                                "Usage:\n" +
                                "> draw_poorly_assembled <contigs_file> <prefix_of_base_assembly> [first N contigs to analyze]\n" +
                                " Draws pictures of contigs that are not well covered with any contig in base assembly.";
                return answer;
            }

            DrawPoorlyAssembledCommand() : DrawingCommand("draw_poorly_assembled")
            {
            }

            void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
                const vector<string>& args = arg_list.GetAllArguments();
                if (!CheckCorrectness(args))
                    return;

                std::string contigs_file = args[1];
                string base_assembly_prefix = args[2];

                if (!CheckFileExists(contigs_file)) {
                    LOG("File with contigs " << contigs_file << " not found");
                }

                size_t contig_cnt = -1u;
                if (args.size() > 3) {
                    LOG("Will analyze first " << args[3] << " contigs");
                    contig_cnt = lexical_cast<size_t>(args[3]);
                }

                io::FileReadStream irs(contigs_file);

                size_t i = 0;
                while (!irs.eof() && i < contig_cnt) {
                    io::SingleRead contig;
                    irs >> contig;
                    contig = MakeValid(contig);
                    LOG("Considering contig " << contig.name());

                    // if read is valid and also the name contains a given string <contig_name> as a substring.
                    VERIFY(contig.IsValid());

                    if (IsPoorlyAssembled(curr_env.graph_pack(), contig, base_assembly_prefix)) {
                        LOG("Was poorly assembled, drawing");
                        DrawContig(curr_env, contig);
                    } else {
                        LOG("Was well assembled");
                    }

                    ++i;
                }

            }

    };
}
