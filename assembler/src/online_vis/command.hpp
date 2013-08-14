//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "environment.hpp"
#include "loaded_environments.hpp"
#include "argument_list.hpp"
#include "errors.hpp"

namespace online_visualization {

template <class Env>
class CommandMapping;

template <class Env>
class Command {
 protected:
  string invocation_string_;

  virtual size_t MinArgNumber() const {
    return 0;
  }

  bool CheckEnoughArguments(const vector<string>& args) const {
    bool result = (args.size() > MinArgNumber());
    if (!result)
      FireNotEnoughArguments();
    return result;
  }

 public:
  virtual string Usage() const = 0;

  Command(string invocation_string)
      : invocation_string_(invocation_string) {
  }

  virtual ~Command() {
  }

  string invocation_string() const {
    return invocation_string_;
  }

  // system command, curr_env can point to null
  virtual void Execute(shared_ptr<Env>& curr_env, LoadedEnvironments<Env>& loaded_environments, const ArgumentList& arg_list) const = 0;

  // virtual void Execute(shared_ptr<Env>& curr_env, const ArgumentList& arg_list) const = 0;

};

template <class Env>
class LocalCommand : public Command<Env> {

 public:
  LocalCommand(string invocation_string) : Command<Env>(invocation_string) {
  }

  // command for the current environment
  virtual void Execute(Env& curr_env, const ArgumentList& arg_list) const = 0;

  // !!!! NO OVERRIDING !!!!
  virtual void Execute(shared_ptr<Env>& curr_env, LoadedEnvironments<Env>& loaded_environments, const ArgumentList& arg_list) const {
    if (arg_list["all"] == "true")
      for (auto iter = loaded_environments.begin(); iter != loaded_environments.end(); ++iter)
        Execute(*(iter->second), arg_list);
    else if (curr_env) {
      Execute(*curr_env, arg_list);
    }
    else
      cout << "The environment is not loaded" << endl;
  }

 protected:

  string TryFetchFolder(Env& curr_env, const ArgumentList& arg_list, size_t arg_nmb = 1) const {
    const vector<string> &args = arg_list.GetAllArguments();

    if (args.size() > arg_nmb) {
      return args[arg_nmb] + "/";
    } else {
      return curr_env.manager().GetDirForCurrentState();
    }
  }
};

template <class Env>
class CommandServingCommand : public Command<Env> {
 protected:
  CommandMapping<Env> *command_container_;

 public:
  CommandServingCommand(string invocation_string, CommandMapping<Env> *command_mapper)
      : Command<Env>(invocation_string),
        command_container_(command_mapper)
  {
  }
};

}
