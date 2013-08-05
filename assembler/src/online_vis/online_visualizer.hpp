//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"
#include "standard_vis.hpp"
#include "command.hpp"
#include "loaded_environments.hpp"
#include "environment.hpp"
#include "autocompletion.hpp"

//#include "all_commands.hpp"
#include "base_commands.hpp"

namespace online_visualization {

template <class Env = Environment>
class OnlineVisualizer {

 public:
  OnlineVisualizer() : command_mapping_() {
  }

  virtual ~OnlineVisualizer() {
  }

  inline void init() {
    string p = path::append_path(cfg::get().load_from, "simplified_graph"); // just for default

    path::make_dir("tmp");
    stringstream ss("load default " + p);
    DEBUG("Adding Commands");
    AddBaseCommands();
    AddSpecificCommands();
    DEBUG("Commands added");
    DEBUG("Adding auto-completion option");
    online_vis_autocompletion::Init(command_mapping_.GetCommandNamesList());
    const Command<Env>& load_command = command_mapping_.GetCommand("load");
    DEBUG("Loading current environment");
    load_command.Execute(current_environment_, loaded_environments_, ss);
    DEBUG("Environment loaded");
  }

  void run() {
    History& history = History::GetHistory();
    bool done = false;

    while (!done) {
      char* line = readline(prompt);
      if (!line)
        exit(1);

      if (*line) {
        string command_with_args(line);
        stringstream ss(command_with_args);
        TRACE("Delegating to the ArgumentList class");
        ArgumentList arg_list(ss);
        string processed_command = arg_list.Preprocess(history);
        if (processed_command == "")
          continue;
        //DEBUG("Processed string " << processed_command);
        string command_string = arg_list.GetAllArguments()[0];
        const Command<Env>& command = command_mapping_.GetCommand(command_string);
        command.Execute(current_environment_, loaded_environments_, arg_list);
        history.AddEntry(processed_command);
        free(line);
      }
    }
  }

 protected:
  void AddCommand(shared_ptr<Command<Env>> command) {
    command_mapping_.AddCommand(command);
  }

  virtual void AddSpecificCommands() {
    AddCommand(shared_ptr<Command<Env> >(new LoadCommand<Env>));
  }


 private:
  static const char* prompt;

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

  shared_ptr<Env> current_environment_;
  LoadedEnvironments<Env> loaded_environments_;
  CommandMapping<Env> command_mapping_;

  DECL_LOGGER("OnlineVisualizer");
};

  template<class Env>
  const char* OnlineVisualizer<Env>::prompt = "GAF$> ";

}
