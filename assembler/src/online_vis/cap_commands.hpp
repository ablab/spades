#pragma once

#include "cap_environment.hpp"
#include "cap_environment_manager.hpp"

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
           "> add_genome <path_to_file> [<genome_name> = <path_to_file>]\n"
           " You should specify path to file in which genome data is stored "
           "(.fasta, .gb, etc.). Also you can provide optional name for genome"
           "to display in future output.\n"
           "For example:\n"
           "> add_genome /home/puperuser/genomes/my_genome.fasta my_genome\n"
           " would add to the environment genome stored in file "
           "`my_genome.fasta` located in folder `/home/puperuser/genomes`\n";
  }

  virtual void Execute(CapEnvironment& curr_env, const ArgumentList& arg_list) const {
    if (!CheckCorrectness(arg_list)) {
      return;
    }
    const vector<string>& args = arg_list.GetAllArguments();
    const std::string &filename = args[1];
    std::string name = filename;
    if (args.size() > 2) {
      name = args[2];
    }

    bool success = curr_env.manager().AddGenomeFromFile(filename, name);
    if (!success) {
      std::cout << "Failed. Genome is not valid. Please check input.\n";
    }
  }

 protected:
  virtual size_t MinArgNumber() const {
    return 1;
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

    std::cout << "Saving genomes to disk.. ";

    std::string dir = curr_env.manager().GetDirForCurrentState();
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

  virtual void Execute(CapEnvironment& curr_env, const ArgumentList& arg_list) const {
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
    std::stringstream ss(args[1]);
    size_t k;
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

class SaveGraphCommand : public LocalCommand<CapEnvironment> {
 public:
  SaveGraphCommand() : LocalCommand<CapEnvironment>("save_graph") {
  }

  virtual std::string Usage() const {
    return "Command `save_graph`\n"
           " Saves graph in common spades format in specified directory.\n"
           " If no directory is specified then default cache directory for current state is used.\n"
           "Usage:\n"
           "> save_graph <directory_to_save_to>\n";
  }
  
  virtual void Execute(CapEnvironment& curr_env, const ArgumentList& arg_list) const {
    if (curr_env.GetGraphK() == CapEnvironment::kNoGraphK) {
      cout << "You should build graph prior to saving it. Aborting.\n";
      return;
    }

    const vector<string> &args = arg_list.GetAllArguments();

    std::string folder;
    if (args.size() > 1) {
      folder = args[1];
    } else {
      folder = curr_env.manager().GetDirForCurrentState();
    }

    cout << "Saving graph in " << folder << " ...";
    curr_env.manager().SaveGraph(folder);
    cout << " Done.\n";
  }

};

}
