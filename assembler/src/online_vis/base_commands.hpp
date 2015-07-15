//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "environment.hpp"
#include "command.hpp"
#include "errors.hpp"
#include "command_mapping.hpp"

namespace online_visualization {

// null
  template <class Env>
    class NullCommand : public Command<Env> {

     public:
      NullCommand() : Command<Env>("null")
      {
      }

      void Execute(shared_ptr<Env>& /* curr_env */,
                   LoadedEnvironments<Env>& /* loaded_environments */,
                   const ArgumentList& /* args */) const
      {
      }

      string Usage() const {
        return "Nothing to do here";
      }
    };

  template <class Env>
    class HelpCommand : public CommandServingCommand<Env> {
      std::string GetCommonUsageString() const {
        std::string answer =
          " Welcome to GAF (Graph Analysis Framework). This framework allows to work with the de Bruijn Graph interactively.\n "
          " You can see the list of command names below. To see a command's help message just type\n"
          "> help <command_name>\n"
          " The list of command names : \n";

        vector<string> command_names = this->command_container_->GetCommandNamesList();
        for (auto it = command_names.begin(); it != command_names.end(); ++it) {
          answer += *it;
          answer += '\n';
        }
        return answer;
      }

      public:
      string Usage() const {
        string answer;
        answer = answer + "The command `help` allows you to see a help message for any command. \n " +
          "Usage: \n" +
          "> help <name_of_command> \n" +
          " Running `help` without parameters yields a list of all commands.";
        return answer;
      }

      HelpCommand(CommandMapping<Env> *command_mapping)
        : CommandServingCommand<Env>("help", command_mapping) {
        }

      void Execute(shared_ptr<Env>& /* curr_env */, LoadedEnvironments<Env>& /* loaded_environments */, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (args.size() == 1) {
          cout << GetCommonUsageString() << endl;
        } else {
          string command_name = args[1];
          const Command<Env>& command = this->command_container_->GetCommand(command_name);
          if (command.invocation_string() == "null")
            return;
          cout << command.Usage() << endl;
        }
      }
    };

  // exit
  template <class Env>
    class ExitCommand : public Command<Env> {
      public:
        string Usage() const {
          return "The command `exit` allows you to exit this application.";
        }

        ExitCommand() :
          Command<Env>("exit")
      {
      }

        void Execute(shared_ptr<Env>& /* curr_env */, LoadedEnvironments<Env>& /* loaded_environments */, const ArgumentList& /* args */) const {
          cout << "Exiting" << endl;
          exit(0);
        }
    };

  // loading new environment from folder with saves
  template <class Env>
    class LoadCommand : public Command<Env> {
      private:
        shared_ptr<Env> MakeNewEnvironment(const string& name, const string& saves, size_t K) const {
          DEBUG("Making new environment " << name);
          shared_ptr<Env> EnvPointer(new Env(name, saves, K));
          DEBUG("Done");
          return EnvPointer;
        }

      protected:
        size_t MinArgNumber() const {
          return 2;
        }

        virtual bool CheckCorrectness(const vector<string>& args, LoadedEnvironments<Env>& loaded_environments) const
        {
          if (!this->CheckEnoughArguments(args))
            return false;

          string path = args[2];
          size_t K;
          if (args.size() > 3) {
            if (!CheckIsNumber(args[3]))
              return false;
            K = GetInt(args[3]);
          } else {
            K = cfg::get().K;
          }
          if (!CheckEnvIsCorrect(path, K))
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

        virtual bool CheckCorrectness(const vector<string>& args) const {
          return this->CheckEnoughArguments(args);
        }

      public:
        string Usage() const {
          string answer;
          answer = answer + "Command `load` \n" +
            "Usage:\n" +
            "> load <environment_name> <path_to_saves> [<k-value>]\n" +
            " You should specify the name of the new environment as well as a path to the graph saves. Optionally, \n" +
            " you can provide a k-value for these saves. \n: " +
            " For example:\n" +
            "> load GraphSimplified data/saves/simplification\n" +
            " would load a new environment with the name `GraphSimplified` from the files\n" +
            " in the folder `data/saves/` with the basename `simplification` (simplification.grp, simplification.sqn, e.t.c).";
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
          if (!CheckCorrectness(args))
            return;

          string name  = args[1];
          string saves = args[2];
          size_t K;
          if (args.size() > 3) {
            K = GetInt(args[3]);
          } else {
            K = cfg::get().K;
          }

          cout << "Loading " << name << " " << saves << endl;
          if (!CheckCorrectness(args, loaded_environments))
            return;

          shared_ptr<Env> new_env = MakeNewEnvironment(name, saves, K);
          loaded_environments.insert(make_pair(name, new_env));
          curr_env = new_env;
        }

    };

  // loading new environment from folder with saves
  template <class Env>
    class SwitchCommand : public Command<Env> {
      protected:
        size_t MinArgNumber() const {
          return 1;
        }

        virtual bool CheckCorrectness(const vector<string>& args) const {
          return this->CheckEnoughArguments(args);
        }

      public:
        string Usage() const {
          string answer;
          answer = answer + "Command `switch` \n" +
            "Usage:\n" +
            " switch <environment_name>\n" +
            " You should specify the name of the environment you want to switch to. For example:\n" +
            "> switch GraphSimplified \n" +
            " would switch you to the environment with the name `GraphSimplified`.\n" +
            " Of course this environment must be loaded first. To see all loaded environments, run command `list`.";
          return answer;
        }

        SwitchCommand() :
          Command<Env>("switch_env")
      {
      }

        void Execute(shared_ptr<Env>& curr_env, LoadedEnvironments<Env>& loaded_environments, const ArgumentList& arg_list) const {
          const vector<string>& args = arg_list.GetAllArguments();

          if (!CheckCorrectness(args))
            return;

          string name = args[1];

          bool okay = false;
          for (auto iterator = loaded_environments.begin(); iterator != loaded_environments.end(); ++iterator) {
            if (name == iterator->first) {
              okay = true;
              curr_env = iterator->second;
              break;
            }
          }
          if (!okay) {
            cout << "Name " << name << " does not exist" << endl;
            cout << "Please try again" << endl;
            return;
          } else
            cout << "Switching to " << name << endl;
        }

    };

  template <class Env>
    class ListCommand : public Command<Env> {
      protected:
        virtual bool CheckCorrectness() const {
          return true;
        }

      public:
        string Usage() const {
          string answer;
          answer = answer + "Command `list` \n" +
            "Usage:\n" +
            "> list\n" +
            " This command lists all loaded environments.";
          return answer;
        }

        ListCommand() : Command<Env>("list")
      {
      }

        void Execute(shared_ptr<Env>& curr_env, LoadedEnvironments<Env>& loaded_environments, const ArgumentList& /* arg_list */) const {
          cout << "Environments :" << endl;
          for (auto iter = loaded_environments.begin(); iter != loaded_environments.end(); ++iter) {
            cout << " " << iter->first << endl;
          }
          if (curr_env)
            cout << "Current environment is " << curr_env->str() << endl;
          else
            cout << "Current environment was not set" << endl;
        }
    };

  template <class Env>
    class ReplayCommand : public CommandServingCommand<Env> {
      private:
        virtual bool CheckCorrectness(const vector<string>& args) const {
          if (args.size() == 1)
            return true;
          return CheckIsNumber(args[1]);
        }

      public:
        string Usage() const {
          string answer;
          answer = answer + "Command `replay` \n" +
            " Usage:\n" +
            " > rep <command_number>\n" +
            " Runs the command <command_number> commands before. For example:\n" +
            " > rep 1 \n" +
            " would run the previous command.\n" +
            " It is still under development.";
          return answer;
        }

        ReplayCommand(CommandMapping<Env> *command_mapping) : CommandServingCommand<Env>("replay", command_mapping)
      {
      }

        void Execute(shared_ptr<Env>& curr_env, LoadedEnvironments<Env>& loaded_environments, const ArgumentList& arg_list) const {
          const vector<string>& args = arg_list.GetAllArguments();

          if (!CheckCorrectness(args))
            return;

          size_t number = GetInt(args[1]);
          if (number == 0 || number > 100000) {
            LOG(number << " is not in the range");
          }

          History& history = History::GetHistory();

          cout << "Executing the command " << number << " command(s) before... " << endl;
          string command_with_args = history[int(history.size() - number)];
          cout << command_with_args << endl;
          //inserting a command, which is to be repeated
          history.SetEntry(int(history.size() - 1), command_with_args);

          stringstream ss(command_with_args);
          TRACE("Delegating to the ArgumentList class");
          ArgumentList tmp_arg_list(ss);
          //inserting a command, which is to be repeated
          string processed_command = tmp_arg_list.Preprocess(history);
          DEBUG("processed string " << processed_command);
          const string& command_string = tmp_arg_list.GetAllArguments()[0];
          const Command<Env>& command = this->command_container_->GetCommand(command_string);
          command.Execute(curr_env, loaded_environments, tmp_arg_list);
          history.AddEntry(command_with_args);
        }
    };

  template <class Env>
    class LogCommand : public Command<Env> {
      private:
        size_t MinArgNumber() const {
          return 0;
        }

        virtual bool CheckCorrectness(const vector<string>& args) const {
          if (args.size() > 1)
            return CheckIsNumber(args[1]);
          return true;
        }
      public:
        string Usage() const {
          string answer;
          answer = answer + "Command `log` \n" +
            " Usage:\n" +
            " > log [<number_of_commands>]\n" +
            " Shows last <number_of_commands> in the history. Shows the whole log by default.";
          return answer;
        }

        LogCommand() : Command<Env>("log")
      {
      }

        void Execute(shared_ptr<Env>& /* curr_env */, LoadedEnvironments<Env>& /* loaded_environments */, const ArgumentList& arg_list) const {
          vector<string> args = arg_list.GetAllArguments();
          if (!CheckCorrectness(args))
            return;

          History& history = History::GetHistory();
          if (args.size() > 1) {
            size_t number = GetInt(args[1]);
            if (number > history.size())
              number = history.size();
            for (size_t i = 0; i < number; ++i)
              cout << " " << history[history.size() - int(number) + i] << endl;
          }
          else {
            for (size_t i = 0; i < history.size(); ++i) {
              cout << history[i] << endl;
            }
          }
        }
    };

  template <class Env>
    class SaveBatchCommand : public Command<Env> {
      private:
        size_t MinArgNumber() const {
          return 2;
        }

        virtual bool CheckCorrectness(const vector<string>& args) const {
          if (!this->CheckEnoughArguments(args))
            return false;
          return CheckIsNumber(args[1]);
        }

      public:
        string Usage() const {
          string answer;
          answer = answer + "Command `save` \n" +
            " Usage:\n" +
            " > save <number_of_commands> <file_name>\n" +
            " Saves last <number_of_commands> of the history in the file filename.";
          return answer;
        }

        SaveBatchCommand() : Command<Env>("save")
      {
      }

        void Execute(shared_ptr<Env>& /* curr_env */, LoadedEnvironments<Env>& /* loaded_environments */, const ArgumentList& arg_list) const {
          const vector<string>& args = arg_list.GetAllArguments();

          if (!CheckCorrectness(args))
            return;

          size_t number = GetInt(args[1]);
          const string& file = args[2];

          ofstream outfile;
          outfile.open(file);
          History& history = History::GetHistory();

          if (number > history.size())
            number = history.size();

          for (size_t i = 0; i < number; ++i) {
            outfile << history[int(history.size() - number + i)];
            if (i < number - 1)
              outfile << endl;
          }
          outfile.close();
        }
    };

  template <class Env>
    class BatchCommand : public CommandServingCommand<Env> {
      private:
        size_t MinArgNumber() const {
          return 1;
        }

        virtual bool CheckCorrectness(const vector<string>& args) const {
          if (!this->CheckEnoughArguments(args))
            return false;
          return true;
        }

      public:
        string Usage() const {
          string answer;
          answer = answer + "Command `batch` \n" +
            "Usage:\n" +
            "> batch <batch_filename>\n" +
            " Runs the commands from the file <batch_filename>.";
          return answer;
        }

        BatchCommand(CommandMapping<Env> *command_mapping) : CommandServingCommand<Env>("batch", command_mapping)
      {
      }

        void Execute(shared_ptr<Env>& curr_env, LoadedEnvironments<Env>& loaded_environments, const ArgumentList& arg_list) const {
          const vector<string>& args = arg_list.GetAllArguments();

          if (!CheckCorrectness(args))
            return;

          const string& file = args[1];
          if (!CheckFileExists(file))
            return;

          ifstream infile;
          infile.open(file);
          History& history = History::GetHistory();
          while (!infile.eof()) {
            string command_with_args;
            getline(infile, command_with_args);
            if (command_with_args == "")
              continue;
            cout << "> " << command_with_args << endl;
            stringstream ss(command_with_args);
            ArgumentList arg_list(ss);
            string processed_command = arg_list.Preprocess(history);

            const string& command_string = arg_list.GetAllArguments()[0];
            const Command<Env>& command = this->command_container_->GetCommand(command_string);
            command.Execute(curr_env, loaded_environments, arg_list);

            history.AddEntry(processed_command);
          }
          infile.close();
        }
    };

}
