#pragma once

#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"
#include "standard_vis.hpp"
#include "command.hpp"
#include "loaded_environments.hpp"
#include "environment.hpp"

//#include "all_commands.hpp"
#include "base_commands.hpp"

namespace online_visualization {

template <class Env = Environment>
class OnlineVisualizer {
 private:

  shared_ptr<Env> current_environment_;
  LoadedEnvironments<Env> loaded_environments_;
  CommandMapping<Env> command_mapping_;

  void AddBaseCommands() {
    AddCommand(shared_ptr<Command<Env> >(new NullCommand<Env>));
    AddCommand(shared_ptr<Command<Env> >(new ExitCommand<Env>));
    AddCommand(shared_ptr<Command<Env> >(new ListCommand<Env>));
    AddCommand(shared_ptr<Command<Env> >(new HelpCommand<Env>(&command_mapping_)));
    
    AddCommand(shared_ptr<Command<Env> >(new LogCommand<Env>));
    AddCommand(shared_ptr<Command<Env> >(new SaveBatchCommand<Env>));
    AddCommand(shared_ptr<Command<Env> >(new BatchCommand<Env>(&command_mapping_)));

    AddCommand(shared_ptr<Command<Env> >(new SwitchCommand<Env>));
    AddCommand(shared_ptr<Command<Env> >(new ReplayCommand<Env>(&command_mapping_)));
  }

 protected:
  void AddCommand(shared_ptr<Command<Env> > command) {
    command_mapping_.AddCommand(command);
  }

  virtual void AddSpecificCommands() {
    AddCommand(shared_ptr<Command<Env> >(new LoadCommand<Env>));
  }

 private:

 public:
  OnlineVisualizer() : command_mapping_() {
  }

  virtual ~OnlineVisualizer() {
  }

  inline void init() {
    string p = path::append_path(cfg::get().load_from, "late_pair_info_counted");

    path::make_dir("tmp");
    stringstream ss("load default " + p);
    DEBUG("Adding Commands");
    AddBaseCommands();
    AddSpecificCommands();
    DEBUG("Commands added");
    const Command<Env>& load_command = command_mapping_.GetCommand("load");
    DEBUG("Loading current environment");
    load_command.Execute(current_environment_, loaded_environments_, ss);
    DEBUG("Environment loaded");
  }

  void run() {
    vector<string>& history = GetHistory();
    //const size_t max_buffer_size = 10000;

    while (true) {
      cout << "GAF$> ";
      string command_with_args;
      getline(cin, command_with_args);
      stringstream ss(command_with_args);
      TRACE("Delegating to the ArgumentList class");
      ArgumentList arg_list(ss);
      string processed_command = arg_list.Preprocess(history);

      DEBUG("processed string " << processed_command);

      const string& command_string = arg_list.GetAllArguments()[0];
      const Command<Env>& command = command_mapping_.GetCommand(command_string);
      command.Execute(current_environment_, loaded_environments_, arg_list);

      history.push_back(processed_command);
    }
  }

 private:
  DECL_LOGGER("OnlineVisualizer");
};

}
