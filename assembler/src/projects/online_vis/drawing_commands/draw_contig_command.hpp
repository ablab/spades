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
class DrawContigCommand : public DrawingCommand {

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
                        "> draw_contig <name_of_contig> <contigs_file>\n" + 
                        " Draws graph pictures for a contig.";
        return answer;
    }

    DrawContigCommand() : DrawingCommand("draw_contig")
    {
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        string contig_name = args[1];
        LOG("Trying to draw contig " << contig_name);

        bool starts_with = false;
        if (contig_name[contig_name.size() - 1] == '*') {
            starts_with = true;
            contig_name = contig_name.substr(0, contig_name.size() - 1);
        }
        string contigs_file = args[2];
        if (!CheckFileExists(contigs_file))
            return;

        io::FileReadStream reader(contigs_file);

        while (!reader.eof()) {
            io::SingleRead read;
            reader >> read;
            //LOG("Contig " << read.name() << " is being processed now");

            // if the name contains a given string <contig_name> as a substring.
            if((starts_with && read.name().find(contig_name) != string::npos) || contig_name == read.name()) {
                DrawPicturesAlongContig(curr_env, read);
            }
        }
    }
};

class DrawContigsCommand : public DrawingCommand {

protected:
    size_t MinArgNumber() const {
        return 1;   
    }
    
    bool CheckCorrectness(const vector<string>& args) const {
        if (!CheckEnoughArguments(args))
            return false;

        return true;
    }

public:
    string Usage() const {
        string answer;
        answer = answer + "Command `draw_contigs` \n" + 
                        "Usage:\n" + 
                        "> draw_contigs <contigs_file>\n" + 
                        " Draws graph pictures for contigs.";
        return answer;
    }

    DrawContigsCommand() : DrawingCommand("draw_contigs")
    {
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        string contigs_file = args[1];

        LOG("Drawing contigs from " << contigs_file);
        if (!CheckFileExists(contigs_file))
            return;

        io::FileReadStream reader(contigs_file);

        while (!reader.eof()) {
            io::SingleRead read;
            reader >> read;
            //LOG("Contig " << read.name() << " is being processed now");

            DrawPicturesAlongContig(curr_env, read);
        }
    }
};
}
