//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "io/reads/wrapper_collection.hpp"

namespace online_visualization {
    class DrawGapWithinContigCommand : public DrawingCommand {

    protected:
        size_t MinArgNumber() const {
            return 4;
        }

        bool CheckCorrectness(const vector<string>& args) const {
            if (!CheckEnoughArguments(args))
                return false;

            return true;
        }

    public:
        string Usage() const {
            string answer;
            answer = answer + "Command `draw_gap_within_contig` \n" +
                     "Usage:\n" +
                     "> draw_gap_within_contig <name_of_contig> <contig_file> <edge1> <edge2> <barcode_threshold>\n" +
                     " Draws graph pictures between two edges of the contig. Edges that contain at least \n" +
                     " barcode_threshold barcodes shared by edge1 and edge2 are highlighted in blue.";
            return answer;
        }

        DrawGapWithinContigCommand() : DrawingCommand("draw_gap_within_contig")
        {
        }

        void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
            const vector<string>& args = arg_list.GetAllArguments();
            if (!CheckCorrectness(args))
                return;

            string contig_name = args[1];
            EdgeId edge1 = curr_env.finder().ReturnEdgeId(GetInt(args[3]));
            EdgeId edge2 = curr_env.finder().ReturnEdgeId(GetInt(args[4]));

            LOG("Trying to draw gap within contig " << contig_name << " between edges "
                                                    << edge1.int_id() << " and " << edge2.int_id());

            bool starts_with = false;
            if (contig_name[contig_name.size() - 1] == '*') {
                starts_with = true;
                contig_name = contig_name.substr(0, contig_name.size() - 1);
            }
            string contigs_file = args[2];
            if (!CheckFileExists(contigs_file))
                return;

            io::FileReadStream reader(contigs_file);

            size_t barcode_threshold = 1;
            if (args.size() == MinArgNumber() + 2) {
                barcode_threshold = static_cast<size_t>(GetInt(args[5]));
            }

            while (!reader.eof()) {
                io::SingleRead read;
                reader >> read;
                if((starts_with && read.name().find(contig_name) != string::npos) || contig_name == read.name()) {
                    DrawPicturesInContigGap(curr_env, read, edge1, edge2, barcode_threshold);
                }
            }
        }
    };
}
