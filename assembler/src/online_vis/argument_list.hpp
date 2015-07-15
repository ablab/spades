//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "vis_utils.hpp"
#include "history.hpp"
#include "errors.hpp"

#include <boost/tokenizer.hpp>

namespace online_visualization {

  class ArgumentList {

    public:
      ArgumentList() {
      }

      ArgumentList(stringstream& stream) {
        DEBUG("Splitting in tokens");
        vector<string> args = SplitInTokens(stream);
        DEBUG("Parsing args");
        ParseArguments(args);
      }

      string operator[](const string& option_name) const {
        // usual option
        if (options.find(option_name) == options.end()) {
          return "null";
        }
        string result = options.find(option_name)->second;
        return result;
      }

      bool contains(string short_opt) const {
        return (short_options.count(short_opt) > 0);
      }

      const vector<string>& GetAllArguments() const {
        return arguments;
      }

      string Preprocess(const History& history) {
        vector<string> new_arguments;

        for (auto iter = arguments.begin(); iter != arguments.end(); ++iter) {
          string arg = *iter;
          TRACE("Argument " << arg);
          if (arg == "!$") {
            TRACE("!$");
            stringstream ss(history.back());
            TRACE("Last command " << ss.str());
            ArgumentList arg_list(ss);
            const vector<string>& args = arg_list.GetAllArguments();
            string new_arg = args[args.size() - 1];
            TRACE("All args " << args);
            TRACE("New arg " << new_arg);
            new_arguments.push_back(new_arg);
          }
          else if (arg[0] == '!') {
            stringstream ss(history.back());
            size_t i = 1;
            if (arg[1] == '-')
              i = 2;
            string num_of_command = "";
            while (i < arg.size() && arg[i] != ':') {
              num_of_command = num_of_command + arg[i];
              ++i;
            }

            if (num_of_command == "")
              num_of_command = "1";
            TRACE("Number of the command " << num_of_command);

            if (IsNumber(num_of_command) && arg[i] == ':') {
              ++i;
              string num_of_arg = "";
              while (i < arg.size()) {
                num_of_arg = num_of_arg + arg[i];
                ++i;
              }
              TRACE("Number of the argument " << num_of_arg);
              if (num_of_arg == "$" || IsNumber(num_of_arg)) {
                int command_num = GetInt(num_of_command);
                if (command_num <= 0 || command_num > int(history.size())) {
                  FireNumberOutOfBounds(command_num);
                  return "";
                }
                stringstream ss(history[int(history.size()) - GetInt(num_of_command)]);
                TRACE("Got the command " << ss.str());
                ArgumentList arg_list(ss);
                string new_arg;
                // $ means the last one
                if (num_of_arg == "$")
                  new_arg = arg_list.GetAllArguments()[arg_list.GetAllArguments().size() - 1];
                else {
                  int arg_num = GetInt(num_of_arg);
                  if (0 <= arg_num && arg_num < int(arg_list.GetAllArguments().size())) {
                    TRACE("Got the argument " << arg_num);
                    new_arg = arg_list.GetAllArguments()[arg_num];
                  } else {
                    FireBadArgument(arg);
                    return "";
                  }
                }
                TRACE("New arg " << new_arg);
                new_arguments.push_back(new_arg);
              }
              else {
                FireBadArgument(arg);
                return "";
              }
            }
            else {
              new_arguments.push_back(arg);
            }
          }
          else {
            new_arguments.push_back(arg);
          }
        }
        arguments = new_arguments;
        stringstream result;
        for (auto iter = options.begin(); iter != options.end(); ++iter)
          result << iter->first + "=" + iter->second + " ";
        for (auto iter = short_options.begin(); iter != short_options.end(); ++iter)
          result << *iter + " ";
        for (size_t i = 0; i < arguments.size(); ++i) {
          result << arguments[i];
          if (i < arguments.size() - 1)
            result << " ";
        }

        return result.str();
      }

    private:
      const vector<string> SplitInTokens(stringstream& args) const {
        vector<string> answer;
        while (!args.eof()) {
          string arg;
          args >> arg;
          boost::char_separator<char> sep (" ,;");
          boost::tokenizer<boost::char_separator<char>> tokens(arg, sep);
          for (auto I = tokens.begin(); I != tokens.end(); ++I) {
            TRACE("Found argument " << *I);
            answer.push_back(*I);
          }
        }
        return answer;
      }

      pair<string, string> ParseOption(const string& arg) const {
        string opt_name;
        string opt_value = "";

        size_t i = 2;
        for (; i < arg.size() && arg[i] != '='; ++i) {
          opt_name = opt_name + arg[i];
        }
        for (; i < arg.size(); ++i) {
          opt_value = opt_value + arg[i];
        }

        TRACE("Name/Value " << opt_name << " " << opt_value);
        if (opt_value == "")
          opt_value = "true";


        return make_pair(opt_name, opt_value);
      }

      vector<string> ParseShortOption(const string& arg) const {
        vector<string> result;
        size_t i = 1;
        for (; i < arg.size(); ++i) {
          string s = "";
          s = s + arg[i];
          result.push_back(s);
        }
        return result;
      }

      void ParseArguments(const vector<string>& args) {
        for (size_t i = 0; i < args.size(); ++i) {
          TRACE("Parsing argument " << args[i]);
          if (args[i][0] == '-' && args[i][1] == '-') {
            //--smth=<smth>
            TRACE("it is an option");
            pair<string, string> opt_val = ParseOption(args[i]);
            options.insert(opt_val);
          }
          else if (args[i][0] == '-') {
            TRACE("it is a short option");
            const vector<string>& short_opt = ParseShortOption(args[i]);
            TRACE("short options in a vector " << short_opt);
            short_options.insert(short_opt.begin(), short_opt.end());
          }
          else {
            TRACE("it is a usual arg");
            arguments.push_back(args[i]);
          }
        }
      }
      map<string, string> options;
      set<string> short_options;
      vector<string> arguments;

      DECL_LOGGER("ArgumentList");
  };

}
