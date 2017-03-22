//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/graph_pack.hpp"
#include "visualization/visualization_utils.hpp"
#include "standard_vis.hpp"
#include "command.hpp"
#include "loaded_environments.hpp"
#include "environment.hpp"
#include "utils/autocompletion.hpp"

//#include "all_commands.hpp"
#include "base_commands.hpp"

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <pthread.h>
#include <atomic>


namespace online_visualization {

std::atomic<bool> ctrlc_handler_;

inline void * wait_for_second_ctrlc(void *) {
    ctrlc_handler_ = true;
    cerr << endl << "Hit Ctrl+C within 1 second once more to exit" << endl;
    sleep(1);
    ctrlc_handler_ = false;
    return NULL;
}

inline void ctrlc_handler(int /*s*/) {
    if (!ctrlc_handler_) {
        pthread_t thread;
        pthread_create( &thread, NULL, wait_for_second_ctrlc, NULL);
    }
    else {
        exit(-1);
    }
}

template <class Env = Environment>
class OnlineVisualizer {
 public:
  OnlineVisualizer() : command_mapping_() {
  }

  virtual ~OnlineVisualizer() {
  }

  inline void init() {
    string p = fs::append_path(cfg::get().load_from, "simplification"); // just for default

    fs::make_dir("tmp");
    DEBUG("Adding Commands");
    AddBaseCommands();
    AddSpecificCommands();
    DEBUG("Commands added");
    DEBUG("Adding auto-completion option");
    utils::InitAutocompletion(command_mapping_.GetCommandNamesList());
    //stringstream ss("load default " + p);
    //const Command<Env>& load_command = command_mapping_.GetCommand("load");
    //DEBUG("Loading current environment");
    //load_command.Execute(current_environment_, loaded_environments_, ss);
    //DEBUG("Environment loaded");
  }

  string read_line() {
    cout << "[end]" << endl;
      char* line = readline(prompt);
      if (!line)
        exit(1);
      string answer(line);
      free(line);
      return answer;
  }

  void run(const string& batch_file = "") {
    ctrlc_handler_ = false;
    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = ctrlc_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, NULL);

    History& history = History::GetHistory();

    string command_with_args;
    if (batch_file != "") {
        command_with_args = "batch " + batch_file;
    } else {
        command_with_args = read_line();
    }
    bool done = false;

    while (!done) {
      if (!command_with_args.empty()) {
        stringstream ss(command_with_args);
        TRACE("Delegating to the ArgumentList class");
        ArgumentList arg_list(ss);
        string processed_command = arg_list.Preprocess(history);
        if (processed_command == "")
          continue;
        //DEBUG("Processed string " << processed_command);
        string command_string = arg_list.GetAllArguments()[0];
        const Command<Env>& command = command_mapping_.GetCommand(command_string);
        DEBUG("Command " << processed_command << " starting to execute");
        command.Execute(current_environment_, loaded_environments_, arg_list);
        DEBUG("Command " << processed_command << " executed");

        history.AddEntry(command_with_args);
        DEBUG("Command " << processed_command << " added to history");

      }
      command_with_args = read_line();
    }
  }

 protected:
  void AddCommand(shared_ptr<Command<Env>> command) {
    command_mapping_.AddCommand(command);
  }

  virtual void AddSpecificCommands() {
  }


 private:
  static const char* prompt;

  void AddBaseCommands() {
    AddCommand(make_shared<NullCommand<Env>>());
    AddCommand(make_shared<ExitCommand<Env>>());
    AddCommand(make_shared<ListCommand<Env>>());
    AddCommand(make_shared<HelpCommand<Env>>(&command_mapping_));

    AddCommand(make_shared<LogCommand<Env>>());
    AddCommand(make_shared<SaveBatchCommand<Env>>());
    AddCommand(make_shared<BatchCommand<Env>>(&command_mapping_));

    AddCommand(make_shared<SwitchCommand<Env>>());
    AddCommand(make_shared<ReplayCommand<Env>>(&command_mapping_));

    //todo think about why it was in the specific commands
    AddCommand(make_shared<LoadCommand<Env>>());
  }

  shared_ptr<Env> current_environment_;
  LoadedEnvironments<Env> loaded_environments_;
  CommandMapping<Env> command_mapping_;

  DECL_LOGGER("OnlineVisualizer");
};

  template<class Env>
  const char* OnlineVisualizer<Env>::prompt = "GAF$> ";

}
