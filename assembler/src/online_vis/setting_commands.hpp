//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "environment.hpp"
#include "command.hpp"
#include "errors.hpp"

namespace online_visualization {

class LoadGenomeCommand : public LocalCommand<DebruijnEnvironment> {

    protected:
        size_t MinArgNumber() const {
            return 1;
        }

        bool CheckCorrectness(const vector<string>& args) const {
            if (!CheckEnoughArguments(args))
                return false;
            const string& file = args[1];
            if (!CheckFileExists(file))
                return false;
            return true;
        }

    public:
        string Usage() const {
            string answer;
            answer = answer + "Command `load_genome` \n" +
                            " Usage:\n" +
                            " > load_genome <path_to_genome>\n" +
                            " You should specify a path to the genome you want to load from.\n" +
                            " Previously loaded genomes would be lost.";
            return answer;
        }

        LoadGenomeCommand() : LocalCommand<DebruijnEnvironment>("load_genome")
        {
        }

        void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
          const vector<string>& args = arg_list.GetAllArguments();
            if (!CheckCorrectness(args))
                return;
            const string& file = args[1];
            auto genome_reader = make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(file));
            io::SingleRead genome;
            (*genome_reader) >> genome;
            curr_env.LoadNewGenome(genome.sequence());
        }
};

class SetMaxVertCommand : public LocalCommand<DebruijnEnvironment> {
    protected:
        size_t MinArgNumber() const {
            return 1;
        }

        bool CheckCorrectness(const vector<string>& args) const {
            if (!CheckEnoughArguments(args))
                return false;
            return CheckIsNumber(args[1]);
        }

    public:
        string Usage() const {
            string answer;
            answer = answer + "Command `set_max_vertices` \n" +
                            "Usage:\n" +
                            "> set_max_vertices <max_vertices> \n" +
                            " You should specify an integer, which is an upper bound for the number of vertices in the picture.";
            return answer;
        }

        SetMaxVertCommand() : LocalCommand<DebruijnEnvironment>("set_max_vertices")
        {
        }

        void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
            const vector<string>& args = arg_list.GetAllArguments();
            if (!CheckCorrectness(args)) {
                return;
            }
            size_t max_v = GetInt(args[1]);
            curr_env.set_max_vertices(max_v);
        }
};

class SetFolderCommand : public LocalCommand<DebruijnEnvironment> {
    protected:
        size_t MinArgNumber() const {
            return 1;
        }

        bool CheckCorrectness(const vector<string>& args) const {
            return CheckEnoughArguments(args);
        }

    public:
        string Usage() const {
            string answer;
            answer = answer + "Command `set_folder` \n" +
                            "Usage:\n" +
                            "> set_folder <folder_name> \n" +
                            " You should specify a string, which is a new name for a pictures' folder.";
            return answer;
        }
        SetFolderCommand() : LocalCommand<DebruijnEnvironment>("set_folder")
        {
        }

        void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
            const vector<string>& args = arg_list.GetAllArguments();
            if (!CheckCorrectness(args))
                return;
            string folder_name = args[1];
            path::make_dirs(folder_name);
            curr_env.set_folder(folder_name);
        }
};

class SetFileNameCommand : public LocalCommand<DebruijnEnvironment> {
    protected:
        size_t MinArgNumber() const {
            return 1;
        }

        bool CheckCorrectness(const vector<string>& args) const {
            return CheckEnoughArguments(args);
        }

    public:
        string Usage() const {
            string answer;
            answer = answer + "Command `set_file_name` \n" +
                            "Usage:\n" +
                            "> set_file_name <file_base_name>\n" +
                            " You should specify a string, which is a new base_name for all the pictures, that you generate.";
            return answer;
        }

        SetFileNameCommand() : LocalCommand<DebruijnEnvironment>("set_file_name")
        {
        }

        void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
            const vector<string>& args = arg_list.GetAllArguments();
            if (!CheckCorrectness(args))
                return;
            string file_name = args[1];
            curr_env.set_file_name(file_name);
        }
};
}

