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
  std::map<std::string, std::shared_ptr<Command<Env>>> command_map_;

 public:
  CommandMapping() : command_map_() {
  }

  const Command<Env> &GetCommand(const std::string &name) const {
    auto it = command_map_.find(name);
    if (it == command_map_.end()) {
      cout << "No such command `" << name << "`, try again" << endl;
      it = command_map_.find("null");
      VERIFY(it != command_map_.end());
    }
    return *(it->second);
  }

  void AddCommand(std::shared_ptr<Command<Env>> command) {
    auto command_invocation_string = command->invocation_string();
    auto it = command_map_.find(command_invocation_string);
    CHECK_FATAL_ERROR(it == command_map_.end(),
               "Cannot add a command with existing name `"
            << command_invocation_string << "'. Program exits.");

    command_map_[command_invocation_string] = command;
  }

  std::vector<std::string> GetCommandNamesList() const {
    std::vector<std::string> result;
    result.reserve(command_map_.size());
    for (const auto &i : command_map_) {
      if (i.first != "null")
        result.push_back(i.first);
    }
    return result;
  }
};

}
