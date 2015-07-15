//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "command.hpp"


namespace online_visualization {

template <class Env>
class CommandMapping {
  map<string, shared_ptr<Command<Env>>> command_map_;

 public:
  CommandMapping() : command_map_() {
  }

  const Command<Env>& GetCommand(string name) const {
    auto it = command_map_.find(name);
    if (it == command_map_.end()) {
      cout << "No such command `" << name << "`, try again" << endl;
      it = command_map_.find("null");
      VERIFY(it != command_map_.end());
    }
    return *(it->second);
  }

  void AddCommand(shared_ptr<Command<Env>> command) {
    string command_invocation_string = command->invocation_string();
    auto it = command_map_.find(command_invocation_string);
    VERIFY_MSG(it == command_map_.end(),
               "Cannot add a command with existing name `"
            << command_invocation_string << "'. Program exits.");

    command_map_[command_invocation_string] = command;
  }

  vector<string> GetCommandNamesList() const {
    vector<string> result;
    result.reserve(command_map_.size());
    for (auto it = command_map_.begin(); it != command_map_.end(); ++it) {
      if (it->first != "null")
        result.push_back(it->first);
    }
    return result;
  }
};

}
