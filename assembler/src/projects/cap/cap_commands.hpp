//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "cap_environment.hpp"
#include "cap_environment_manager.hpp"
#include "mosaic.hpp"
#include "io/reads/sequence_reader.hpp"
#include "utils/path/path_helper.hpp"

namespace online_visualization {

class AddGenomeCommand : public LocalCommand<CapEnvironment> {
 public:
  AddGenomeCommand() : LocalCommand<CapEnvironment>("add_genome") {
  }

  virtual std::string Usage() const {
    return "Command `add_genome`\n"
           " Adds genome contained in file specified to the environment. This involves:\n"
           "  * further tracking of genome (incl. refinement, searching for diffs, drawing pics, etc.\n"
           "  * storing all genome sequence in RAM throughout usage of current environment\n"
           "Usage:\n"
           "> add_genome <path_to_file> <genome_name> [<crop_repeats> (Yy)] \n"
           " You should specify path to file in which genome data is stored "
           "(.fasta, .gb, etc.). Also you should provide name for genome"
           "to display in future output.\n"
           "For example:\n"
           "> add_genome /home/puperuser/genomes/my_genome.fasta my_genome\n"
           " would add to the environment genome stored in file "
           "`my_genome.fasta` located in folder `/home/puperuser/genomes`\n"
           "Optionally N's other strange symbols and repeat families marked by programs s.a. RepeatMasker (written in small letters)"
           " can be ommited without loss of original coordinates\n";
  }

  virtual void Execute(CapEnvironment& curr_env, const ArgumentList& arg_list) const {
    if (!CheckCorrectness(arg_list)) {
      return;
    }
    const vector<string>& args = arg_list.GetAllArguments();
    const std::string &filename = args[1];
    std::string name = filename;
    name = args[2];
    bool crop_repeats = false;
    if (args.size() > 3) {
        VERIFY(args[3] == "Y" || args[3] == "y");
        crop_repeats = true;
        std::cout << "Repeat crop enabled! All small letters will be ignored with coordinated preserved\n";
    }

    bool success = curr_env.manager().AddGenomeFromFile(filename, name, crop_repeats);
    if (!success) {
      std::cout << "Failed. Genome is not valid. Please check input.\n";
    }
  }

 protected:
  virtual size_t MinArgNumber() const {
    return 2;
  }

  virtual bool CheckCorrectness(const ArgumentList& arg_list) const {
    const vector<std::string> &args = arg_list.GetAllArguments();
    if (!CheckEnoughArguments(args)) {
      std::cout << "Command takes one or more arguments. Aborting.\n";
      return false;
    }

    const std::string &filename = args[1];
    if (!CheckFileExists(filename)) {
      std::cout << "There is no file linked to the path given. Aborting.\n";
      return false;
    }

    return true;
  }
};

class SaveGenomesCommand : public LocalCommand<CapEnvironment> {
 public:
  SaveGenomesCommand() : LocalCommand<CapEnvironment>("save_genomes") {
  }

  virtual std::string Usage() const {
    return "Command `save_genomes`\n"
           " Saves all progress of refining the genomes.\n"
           " Namely, stores refined (modified) sequences on hard drive.\n"
           "Usage:\n"
           "> save_genomes [force]\n"
           " `force` is optional parameter. If `force` is /y|Y|(force)/ then\n"
           " genomes will be written even if this combination of genomes and Ks\n"
           " was stored before/\n";
  }

  virtual void Execute(CapEnvironment& curr_env, const ArgumentList& arg_list) const {
    const vector<string> &args = arg_list.GetAllArguments();

    bool force = false;
    if (args.size() > 1 && (args[1] == "y" || args[1] == "Y" || args[1] == "force")) {
      force = true;
    }

    std::string dir = curr_env.manager().GetDirForCurrentState();

    std::cout << "Saving genomes to " << dir << " ...";
    if (cap::utils::DirExist(dir)) {
      std::cout << "Looks like current state was already stored. ";
      if (force) {
        std::cout << "(!) FORCED WRITE..";
      } else {
        std::cout << "Omitting this stage.\n";
        return;
      }
    }
    curr_env.manager().SaveGenomesToDisk(force);

    cout << "Done.\n";
  }

 protected:
  virtual size_t MinArgNumber() const {
    return 0;
  }

};

class RefineCommand : public LocalCommand<CapEnvironment> {
 public:
  RefineCommand() : LocalCommand<CapEnvironment>("refine") {
  }

  virtual std::string Usage() const {
    return "Command `refine`\n"
           " Refines multicolored De Bruijn graph built from all genomes in environment with last chosen K.\n"
           " Some K should be selected and graph built before running this command (see `build_graph`)\n"
           "Usage:\n"
           "> refine\n";
  }

  virtual void Execute(CapEnvironment& curr_env, const ArgumentList& /* arg_list */) const {
    if (curr_env.GetGraphK() == CapEnvironment::kNoGraphK) {
      cout << "Graph has not yet been constructed, aborting.\n";
      return;
    }

    cout << "Refining graph..";
    curr_env.manager().Refine();
    cout << " Done.\n";
  }

 protected:
  virtual size_t MinArgNumber() const {
    return 0;
  }

};

class BuildGraphCommand : public LocalCommand<CapEnvironment> {
 public:
  BuildGraphCommand() : LocalCommand<CapEnvironment>("build_graph") {
  }

  virtual std::string Usage() const {
    return "Command `build_graph`\n"
           " Sets K for multicolored De Bruijn graph and builds graph from genomes previously added to environment (see `add_genome`)\n"
           " K should be odd.\n"
           "Usage:\n"
           "> build_graph <k>\n";
  }

  virtual void Execute(CapEnvironment& curr_env, const ArgumentList& arg_list) const {
    const vector<string> &args = arg_list.GetAllArguments();

    if (!CheckEnoughArguments(args)) {
      return;
    }
    unsigned k;

    std::stringstream ss(args[1]);
    ss >> k;

    if (k % 2 == 0) {
      cout << "K should be odd. Aborting.\n";
      return;
    }

    cout << "Building graph..";
    curr_env.manager().ConstructGraph(k);
    cout << " Done.\n";
  }

 protected:
  virtual size_t MinArgNumber() const {
    return 1;
  }
};

template <>
class LoadCommand<CapEnvironment> : public Command<CapEnvironment> {
 private:
  typedef CapEnvironment Env;
  shared_ptr<Env> MakeNewEnvironment(const string& name, const string &desc) const {
    DEBUG("Making new environment " << name);
    shared_ptr<Env> EnvPointer(new Env(name, desc));
    DEBUG("Done");
    return EnvPointer;
  }

 protected:
  size_t MinArgNumber() const {
    return 1;
  }

  bool CheckCorrectness(const vector<string>& args, LoadedEnvironments<Env>& loaded_environments) const
  {
    if (!this->CheckEnoughArguments(args))
      return false;

    string name = args[1];
    for (auto iterator = loaded_environments.begin(); iterator != loaded_environments.end(); ++iterator) {
      if (name == iterator->first) {
        cout << "Name " << name << " already exists" << endl;
        cout << "Maybe you want to switch to this environment? " << name << endl;
        cout << "Please try again" << endl;
        return false;
      }
    }
    return true;
  }

 public:
  string Usage() const {
    string answer;
    answer = answer + "Command `load` \n" +
      "Usage:\n" +
      "> load <environment_name> [description]\n" +
      " You should specify the name of the new environment. All data and cache concerning \n" +
      " this environment will be stored in " + cap_cfg::get().cache_root + "/<environment_name>/\n" +
      " See cap_config for changing cache root folder.\n";
    return answer;
  }

  LoadCommand() : Command<Env>("load")
  {
  }

  void Execute(shared_ptr<Env>& curr_env,
               LoadedEnvironments<Env>& loaded_environments,
               const ArgumentList& arg_list) const
  {
    vector<string> args = arg_list.GetAllArguments();
    string name  = args[1];
    string desc = "";
    for (size_t i = 2; i < args.size(); ++i) {
      if (i > 2) {
        desc += " ";
      }
      desc += args[i];
    }

    cout << "Loading " << name << endl;
    if (!CheckCorrectness(args, loaded_environments))
      return;

    shared_ptr<Env> new_env = MakeNewEnvironment(name, desc);
    loaded_environments.insert(make_pair(name, new_env));
    curr_env = new_env;
  }

};

class SaveEnvCommand : public NewLocalCommand<CapEnvironment> {
 public:
  SaveEnvCommand() : NewLocalCommand<CapEnvironment>("save_env", 1) {
  }

  virtual std::string Usage() const {
    return "Command `save_env`\n"
           "Usage:\n"
           "> save_graph <directory_to_save_to>\n";
  }

 private:
  virtual void InnerExecute(CapEnvironment& curr_env, const vector<string>& args) const {
    std::string folder;
    folder = args[1] + "/";

    cout << "Saving env in " << folder << " ...";

    cap::utils::MakeDirPath(folder);
    VERIFY(cap::utils::DirExist(folder));

    std::ofstream write_stream(folder + "environment");
    curr_env.WriteToStream(write_stream);
    write_stream.close();
    cout << " Done.\n";
  }

};

class LoadEnvCommand : public NewLocalCommand<CapEnvironment> {
 public:
  LoadEnvCommand() : NewLocalCommand<CapEnvironment>("load_env", 1) {
  }

  virtual std::string Usage() const {
    return "Command `load_env`\n"
           "Usage:\n"
           "> load_env <directory with save>\n";
  }

private:
  virtual void InnerExecute(CapEnvironment& curr_env, const vector<string>& args) const {
    std::string folder;
    VERIFY(args.size() > 1);
    folder = args[1] + "/";

    cout << "Load env from " << folder << " ...";

    std::ifstream read_stream(folder + "environment");
    curr_env.ReadFromStream(read_stream);
    read_stream.close();
    cout << " Done.\n";
  }

};

class SaveGraphCommand : public NewLocalCommand<CapEnvironment> {
 public:
  SaveGraphCommand() : NewLocalCommand<CapEnvironment>("save_graph", 0) {
  }

  virtual std::string Usage() const {
    return "Command `save_graph`\n"
           " Saves graph in common spades format in specified directory.\n"
           " If no directory is specified then default cache directory for current state is used.\n"
           "Usage:\n"
           "> save_graph <directory_to_save_to>\n";
  }

  virtual void InnerExecute(CapEnvironment& curr_env, const vector<string>& args) const {
    if (curr_env.GetGraphK() == CapEnvironment::kNoGraphK) {
      cout << "You should build graph prior to saving it. Aborting.\n";
      return;
    }

    string folder = TryFetchFolder(curr_env, args);

    cout << "Saving graph in " << folder << " ...";
    curr_env.manager().SaveGraph(folder + "saves/");
    cout << " Done.\n";
  }

};

class DrawPicsCommand : public LocalCommand<CapEnvironment> {
 public:
  DrawPicsCommand() : LocalCommand<CapEnvironment>("draw_pics") {
  }

  virtual std::string Usage() const {
    return "Command `draw_pics`\n"
           " Draws colored graph components in in specified directory.\n"
           " If no directory is specified then default cache directory for current state is used.\n"
           "Usage:\n"
           "> draw_pics <directory_to_save_to>\n";
  }

  virtual void Execute(CapEnvironment& curr_env, const ArgumentList& arg_list) const {
    if (curr_env.GetGraphK() == CapEnvironment::kNoGraphK) {
      cout << "You should build graph prior to saving it. Aborting.\n";
      return;
    }

    std::string folder = TryFetchFolder(curr_env, arg_list);

    cout << "Drawing pics in " << folder << " ...";
    curr_env.manager().DrawPics(folder + "pics/");
    cout << " Done.\n";
  }

};

class FindIndelsCommand : public LocalCommand<CapEnvironment> {
 public:
  FindIndelsCommand() : LocalCommand<CapEnvironment>("find_indels") {
  }

  virtual std::string Usage() const {
    return "Command `find_indels`\n"
           " Finds common in-del events that transform genomes into each other and writes out them.\n"
           " If no output file is specified, the results are written to the default file (see `log_file`)\n"
           " Also there is a feature to mask found indels in graph (!!! this does not affect sequences)\n"
           " Note that graph should be built prior to finding indel events\n"
           "Usage:\n"
           "> find_indels [<mask>=N [<filename>]]\n"
           "Where\n"
           " <mask> is either Y or N\n"
           "For example:\n"
           "> find_indels Y ./indels/log5.txt\n"
           "NOTE: when output file is specified it is overwritten if exists\n";
  }

  virtual void Execute(CapEnvironment &curr_env, const ArgumentList &arg_list) const {
    if (curr_env.GetGraphK() == CapEnvironment::kNoGraphK) {
      cout << "You should build graph prior to saving it. Aborting.\n";
      return;
    }

    const vector<string> &args = arg_list.GetAllArguments();

    bool mask_indels = false;
    std::string filename = curr_env.event_log_path();
    std::string mode = curr_env.event_log_file_mode();

    if (args.size() > 1) {
      VERIFY(args[1].size());
      if (args[1][0] == 'Y' || args[1][0] == 'y') {
        mask_indels = true;
      }
    }
    if (args.size() > 2) {
      filename = args[2];
      mode = "w";
    }

    int code = curr_env.manager().FindIndels(mask_indels, filename, mode);
    if (code == 1) {
      cout << "Output file could not be opened for writing. Aborting.\n";
    }
  }

};

class FindInversionsCommand : public LocalCommand<CapEnvironment> {
 public:
  FindInversionsCommand() : LocalCommand<CapEnvironment>("find_inversions") {
  }

  virtual std::string Usage() const {
    return "Command `find_inversions`\n"
           " Finds common inversion events that transform genomes into each other and writes out them (actually, NO).\n"
           " Note that graph should be built prior to finding inversion events\n"
           "Usage:\n"
           "> find_inversions\n";
  }

  virtual void Execute(CapEnvironment &curr_env, const ArgumentList &/* arg_list */) const {
    if (curr_env.GetGraphK() == CapEnvironment::kNoGraphK) {
      cout << "You should build graph prior to saving it. Aborting.\n";
      return;
    }

//    const vector<string> &args = arg_list.GetAllArguments();?

    bool mask_inversions = false;
    std::string filename = curr_env.event_log_path();
    std::string mode = curr_env.event_log_file_mode();

    /*
    if (args.size() > 1) {
      VERIFY(args[1].size());
      if (args[1][0] == 'Y' || args[1][0] == 'y') {
        mask_indels = true;
      }
    }
    if (args.size() > 2) {
      filename = args[2];
      mode = "w";
    }
    */

    /*int code = */curr_env.manager().FindInversions(mask_inversions, filename, mode);
    /*
    if (code == 1) {
      cout << "Output file could not be opened for writing. Aborting.\n";
    }
    */
  }

};

class BlocksToGRIMMFormat : public LocalCommand<CapEnvironment> {
  public:
    BlocksToGRIMMFormat() : LocalCommand<CapEnvironment>("blocks_to_grimm") {
    }

    virtual std::string Usage() const {
        return "Command `block_to_grimm`\n"
               " Converts blocks output by `save_blocks' to GRIMM format.\n"
               "Usage:\n"
               "> block_to_grimm <blocks_file> <grimm_file>\n";
    }

    virtual void Execute(CapEnvironment &/* curr_env */, const ArgumentList &arg_list) const {
      const vector<string> &args = arg_list.GetAllArguments();

      if (args.size() <= 2) {
        cerr << "Not emough arguments" << endl;
        return;
      }

      std::string file_from = args[1],
                  file_to = args[2];

      file_from = fs::make_full_path(file_from);
      file_to = fs::make_full_path(file_to);

      std::string dir = fs::parent_path(file_to);
      cap::utils::MakeDirPath(dir);

      BlockPrinter<Graph>::ConvertBlocksToGRIMM(file_from, file_to);
    }
};

class SaveBlocksCommand : public LocalCommand<CapEnvironment> {
 public:
    SaveBlocksCommand() : LocalCommand<CapEnvironment>("save_blocks") {
  }

  virtual std::string Usage() const {
    return "Command `save_blocks`\n"
           " Saves all trivial synteny blocks (aka graph edges).\n"
           " Synteny blocks are given new ids (with edge ids also in the file).\n"
           " All the coordinates ()\n"
           "Usage:\n"
           "> save_blocks <file_to_save_to> [unique]\n"
           "Where\n"
           " [unique] if set and equals to (unique|Y|y) then only blocks\n"
           " that appear exactly once in the contigs will be reported.\n";
  }

  virtual void Execute(CapEnvironment& curr_env, const ArgumentList& arg_list) const {
      const vector<string> &args = arg_list.GetAllArguments();
      const std::string folder = TryFetchFolder(curr_env, arg_list);

      bool unique = false;
      if (args.size() > 2 && (args[2] == "y" || args[2] == "Y" || args[2] == "unique")) {
          unique = true;
      }
      INFO("unique = " << unique << ", args[2] = " << args[2]);

      BlockPrinter<Graph> *printer;

      if (!unique) {
        printer = new BlockPrinter<Graph>(curr_env.graph(),
          curr_env.coordinates_handler(), folder + "blocks.txt");
      } else {
        vector<pair<size_t, size_t>> rc_pairs = PrepareRCContigPairs(curr_env);
        printer = new UniqueBlockPrinter<Graph>(curr_env.graph(),
          curr_env.coordinates_handler(), folder + "blocks.txt", rc_pairs);
      }

      for (unsigned i = 0; i < curr_env.genome_cnt(); ++i) {
          printer->ProcessContig(i, 2*i, curr_env.genome_names()[i]);
      }

      delete printer;
  }

 protected:
  virtual size_t MinArgNumber() const {
    return 1;
  }

 private:
  vector<pair<size_t, size_t>> PrepareRCContigPairs(const CapEnvironment &curr_env) const {
    size_t num_contigs = curr_env.genome_cnt();

    vector<pair<size_t, size_t>> res;
    for (size_t i = 0; i < num_contigs; ++i) {
      res.push_back(make_pair(2 * i, 2 * i + 1));
    }

    return res;
  }

};

class LoadGraphCommand : public LocalCommand<CapEnvironment> {
 public:
  LoadGraphCommand() : LocalCommand<CapEnvironment>("load_graph") {
  }

  virtual std::string Usage() const {
    return "Command `load_graph`\n"
           " Loads graph from previously saved saves\n"
           "Usage:\n"
           "> load_graph <K> <path>\n"
           "For example:\n"
           "> find_indels 55 ./masked/graph\n";
  }

  virtual void Execute(CapEnvironment &curr_env, const ArgumentList &arg_list) const {
    const vector<string> &args = arg_list.GetAllArguments();

    uint K = 21;
    stringstream ss(args[1]);
    ss >> K;
    const std::string &path = args[2];

    curr_env.manager().LoadGraphFromSaves(K, path);
  }

 protected:
  virtual size_t MinArgNumber() const {
    return 2;
  }

};

class MosaicAnalysisCommand : public NewLocalCommand<CapEnvironment> {
 public:
    MosaicAnalysisCommand() : NewLocalCommand<CapEnvironment>("mosaic", 0) {
  }

  virtual std::string Usage() const {
    return "Command `mosaic`";
  }

 private:
  virtual void InnerExecute(CapEnvironment& curr_env, const vector<string>& args) const {
      VERIFY(curr_env.genome_cnt() == 1);
//      const Sequence& genome = curr_env.genomes()[1];
      const Sequence& genome = curr_env.genomes()[0];
      size_t min_support_length = 100;
      size_t max_support_mult = 10;
      size_t max_inter_length = 1000;
      size_t min_reportable_mosaic_length = 500;
      size_t min_reportable_submosaic_length = 100;
      std::string folder = TryFetchFolder(curr_env, args);
      cout << "Mosaic analysis triggered" << endl;
      cout << "Min support block length " << min_support_length << endl;
      cout << "Max support block multiplicity " << max_support_mult << endl;
      cout << "Max inter-block length " << max_inter_length << endl;
      if (curr_env.LSeqIsUsed()) {
          VERIFY(false);
//          mosaic::PerformMosaicAnalysis(curr_env.l_seq_gp(), curr_env.coordinates_handler().AsMappingPath(0),
//                                        genome, min_support_length, max_support_mult, max_inter_length,
//                                        min_reportable_mosaic_length,
//                                        min_reportable_submosaic_length, out);
      } else {
          mosaic::PerformMosaicAnalysis(curr_env.rt_seq_gp(), curr_env.coordinates_handler().AsMappingPath(0),
                                        genome, min_support_length, max_support_mult, max_inter_length,
                                        min_reportable_mosaic_length,
                                        min_reportable_submosaic_length, folder);
      }
  }
};

//todo works for finished genomes, not contigs!!!
ContigStreams ConvertRefsToStreams(const vector<Sequence>& ss, const vector<string>& names) {
    ContigStreams answer;
    VERIFY(ss.size() == names.size());
    for (size_t i = 0; i < ss.size(); ++i) {
        answer.push_back(make_shared<io::SequenceReadStream<Contig>>(ss[i], names[i]));
    }
    return answer;
}

class MaskRepeatsCommand : public NewLocalCommand<CapEnvironment> {
public:
    MaskRepeatsCommand()
            : NewLocalCommand<CapEnvironment>("mask_repeats", 2) {
    }

    virtual std::string Usage() const {
        return "Command `mask_repeats <k> <max_iter_count>`";
    }

private:

    vector<string> AppendFasta(const vector<string>& files) const {
        vector<string> answer;
        for (string s : files) {
            answer.push_back(s + ".fasta");
        }
        return answer;
    }

    Sequence ReadSequence(ContigStream& reader) const {
        VERIFY(!reader.eof());
        io::SingleRead read;
        reader >> read;
        return read.sequence();
    }

    void UpdateGenomes(ContigStreams streams, CapEnvironment& curr_env) const {
        vector<Sequence>& genomes = curr_env.mutable_genomes();
        VERIFY(streams.size() == genomes.size());
        for (size_t i = 0; i < streams.size(); ++i) {
            genomes[i] = ReadSequence(streams[i]);
        }
    }

    /*virtual*/
    void InnerExecute(CapEnvironment& curr_env,
                      const vector<string>& args) const {
        size_t k = GetInt(args[1]);
        size_t iteration_cnt = GetInt(args[2]);

        cout << "Masking repeats for k=" << k << " in " << iteration_cnt << " iterations" << endl;

        ContigStreams streams = ConvertRefsToStreams(
                curr_env.genomes(), curr_env.genome_names());

        //todo temporary hack
        curr_env.manager().SaveGenomesToDisk(false);

        string folder = this->CurrentFolder(curr_env) + "masking/";
        make_dir(folder);
        bool success = MaskRepeats(k, streams, AppendFasta(curr_env.genome_names()),
                                   iteration_cnt, folder);
        if (!success) {
            cout << "Failed to mask repeats in " << iteration_cnt
                    << " iterations" << endl;
        } else {
            cout << "Repeats successfully masked" << endl;
            cout << "Updating genomes in environment" << endl;
            UpdateGenomes(OpenStreams(CurrentFolder(curr_env) + "masking/masked/", AppendFasta(curr_env.genome_names()), false), curr_env);
        }
    }

};

}
