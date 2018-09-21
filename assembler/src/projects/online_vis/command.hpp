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

typedef std::vector<std::string> Args;

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
    std::string invocation_string_;


  bool CheckEnoughArguments(const Args &args) const {
    bool result = (args.size() > MinArgNumber());
    if (!result)
      FireNotEnoughArguments();
    return result;
  }

 public:
  virtual std::string Usage() const = 0;

  Command(std::string invocation_string)
      : invocation_string_(std::move(invocation_string)) {
  }

  virtual ~Command() {
  }

  const std::string &invocation_string() const {
    return invocation_string_;
  }

  // system command, curr_env can point to null
  virtual void Execute(std::shared_ptr<Env> &curr_env,
                       LoadedEnvironments<Env> &loaded_environments,
                       const ArgumentList &arg_list) const = 0;

  // virtual void Execute(shared_ptr<Env>& curr_env, const ArgumentList& arg_list) const = 0;

};

//todo reduce code duplication in cap's test_utils
void MakeDirPath(const std::string &path) {
  if (path.size() == 0) {
    TRACE("Somewhat delirium: trying to create directory ``");
    return;
  }

  size_t slash_pos = 0;
  while ((slash_pos = path.find_first_of('/', slash_pos + 1)) != std::string::npos) {
    fs::make_dir(path.substr(0, slash_pos));
  }
  if (path[path.size() - 1] != '/') {
    fs::make_dir(path);
  }
}

inline bool DirExist(const std::string &path) {
  struct stat st;
  return (stat(path.c_str(), &st) == 0) && (S_ISDIR(st.st_mode));
}

template <class Env>
class LocalCommand : public Command<Env> {

 public:
  LocalCommand(std::string invocation_string) :
          Command<Env>(std::move(invocation_string)) {
  }

  // command for the current environment
  virtual void Execute(Env& curr_env, const ArgumentList& arg_list) const = 0;

  // !!!! NO OVERRIDING !!!!
  virtual void Execute(std::shared_ptr<Env> &curr_env,
                       LoadedEnvironments<Env> &loaded_environments,
                       const ArgumentList &arg_list) const {
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

  std::string TryFetchFolder(Env &curr_env, const Args &args, size_t arg_nmb = 1) const {
    if (args.size() > arg_nmb) {
      return MakeDirIfAbsent(args[arg_nmb] + "/");
    } else {
      return MakeDirIfAbsent(CurrentFolder(curr_env));
    }
  }

  std::string TryFetchFolder(Env &curr_env, const ArgumentList &arg_list, size_t arg_nmb = 1) const {
      const auto &args = arg_list.GetAllArguments();
      return TryFetchFolder(curr_env, args, arg_nmb);
  }

  std::string CurrentFolder(Env& curr_env) const {
    return curr_env.manager().GetDirForCurrentState();
  }

private:
  std::string MakeDirIfAbsent(const string& folder) const {
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
    NewLocalCommand(std::string invocation_string, size_t min_arg_num)
            : LocalCommand<Env>(std::move(invocation_string)), min_arg_num_(min_arg_num) {
    }

    // command for the current environment
    /*virtual*/ void Execute(Env& curr_env, const ArgumentList& arg_list) const {
        const auto &args = arg_list.GetAllArguments();
        if (!this->CheckEnoughArguments(args))
            return;
        InnerExecute(curr_env, args);
    }

private:

    virtual size_t MinArgNumber() const {
      return min_arg_num_;
    }

    virtual void InnerExecute(Env &curr_env, const Args &args) const = 0;

};

template <class Env>
class CommandServingCommand : public Command<Env> {
 protected:
  CommandMapping<Env> *command_container_;

 public:
  CommandServingCommand(std::string invocation_string, CommandMapping<Env> *command_mapper)
      : Command<Env>(std::move(invocation_string)),
        command_container_(command_mapper)
  {
  }
};

}
