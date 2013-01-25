//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard_vis.hpp"

namespace online_vis_autocompletion {

  vector<string> COMMANDS;

  void initialize_readline();

  void Init(vector<string> commands) {
    COMMANDS = commands;
    initialize_readline(); 
  }

  inline char* dupstr(const char* s) {
    char* r;
    r = (char*) malloc(strlen(s) + 1);
    strcpy(r, s);

    return r;
  }

  char* command_generator(const char* text, int state) {
    static queue<string> list_possible_matches;

    if (state != 0) {
      if (!list_possible_matches.empty()) {
        char* answer = dupstr(list_possible_matches.front().c_str());
        list_possible_matches.pop();
        return answer;
      }
      else
        return NULL;
    }
    else {
      for (size_t i = 0; i < COMMANDS.size(); ++i) {
        string name = COMMANDS[i];
        if (!name.compare(0, strlen(text), text))
          list_possible_matches.push(name);
      }
      return command_generator(text, 1);
    }
    return NULL;
  }

  char** gaf_completion(const char* text, int start, int end) {
    if (start == 0) {
      typedef char* (*func_ptr)(const char*, int);
      func_ptr func = command_generator;
      return rl_completion_matches(text, func);
    }
    else
      return NULL;
  }

  void initialize_readline() {
    typedef char** (*func_ptr)(const char*, int, int);
    func_ptr func = gaf_completion;
    rl_attempted_completion_function = func;
  }
}
