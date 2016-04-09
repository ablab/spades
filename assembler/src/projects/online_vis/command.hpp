//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "standard_vis.hpp"
#include "environment.hpp"
#include "loaded_environments.hpp"
#include "argument_list.hpp"
#include "errors.hpp"

namespace online_visualization {

template <class Env>
class CommandMapping;

template <class Env>
class Command {
 private:
    virtual size_t MinArgNumber() const {
      return 0;
    }

 //todo fix modifier
 protected:
  string invocation_string_;


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

//todo reduce code duplication in cap's test_utils
void MakeDirPath(const std::string& path) {
  if (path.size() == 0) {
    TRACE("Somewhat delirium: trying to create directory ``");
    return;
  }

  size_t slash_pos = 0;
  while ((slash_pos = path.find_first_of('/', slash_pos + 1)) != std::string::npos) {
    make_dir(path.substr(0, slash_pos));
  }
  if (path[path.size() - 1] != '/') {
    make_dir(path);
  }
}

bool DirExist(std::string path) {
  struct stat st;
  return (stat(path.c_str(), &st) == 0) && (S_ISDIR(st.st_mode));
}

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

  string TryFetchFolder(Env& curr_env, const vector<string>& args, size_t arg_nmb = 1) const {
    if (args.size() > arg_nmb) {
      return MakeDirIfAbsent(args[arg_nmb] + "/");
    } else {
      return MakeDirIfAbsent(CurrentFolder(curr_env));
    }
  }

  string TryFetchFolder(Env& curr_env, const ArgumentList& arg_list, size_t arg_nmb = 1) const {
      const vector<string>& args = arg_list.GetAllArguments();
      return TryFetchFolder(curr_env, args, arg_nmb);
  }

  string CurrentFolder(Env& curr_env) const {
    return curr_env.manager().GetDirForCurrentState();
  }

private:
  string MakeDirIfAbsent(const string& folder) const {
      if (!DirExist(folder))
          MakeDirPath(folder);
      return folder;
  }

};

//todo integrate into basic LocalCommand (after iteratively switching to it in all commands)
template <class Env>
class NewLocalCommand : public LocalCommand<Env> {
    size_t min_arg_num_;

public:
    NewLocalCommand(string invocation_string, size_t min_arg_num)
            : LocalCommand<Env>(invocation_string), min_arg_num_(min_arg_num) {
    }

    // command for the current environment
    /*virtual*/ void Execute(Env& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!this->CheckEnoughArguments(args))
            return;
        InnerExecute(curr_env, args);
    }

private:

    virtual size_t MinArgNumber() const {
      return min_arg_num_;
    }

    virtual void InnerExecute(Env& curr_env, const vector<string>& args) const = 0;

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
