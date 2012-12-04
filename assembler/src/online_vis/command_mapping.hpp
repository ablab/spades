#pragma once

#include "command.hpp"


namespace online_visualization {

template <class Env>
class CommandMapping {
  map<string, shared_ptr<Command<Env> > > command_map_;

 public:
  CommandMapping() : command_map_() {
  }

  const Command<Env> &GetCommand(string name) const {
    auto it = command_map_.find(name);
    if (it == command_map_.end()) {
      cout << "No such command `" << name << "`, try again" << endl;
      it = command_map_.find("null");
      VERIFY(it != command_map_.end());
    }
    return *(it->second);
  }

  void AddCommand(shared_ptr<Command<Env> > command) {
    string command_invocation_string = command->invocation_string();
    auto it = command_map_.find(command_invocation_string);
    if (it != command_map_.end()) {
      cout << "You are trying to add two commands with the same invocation"
              " string! Check for typos." << endl;
      VERIFY(false);
    }

    command_map_[command_invocation_string] = command;
  }

  vector<string> GetCommandNamesList() const {
    vector<string> result;
    result.reserve(command_map_.size());
    for (auto it = command_map_.begin(); it != command_map_.end(); ++it) {
      result.push_back(it->first);
    }
    return result;
  }
};

}
