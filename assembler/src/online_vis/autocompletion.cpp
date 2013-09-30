#include <vector>
#include <string>
#include "autocompletion.hpp"

namespace online_vis{

  std::vector<std::string> commands;

  char* command_generator(const char* text, int state) {
    static queue<string> list_possible_matches;

    if (state != 0) {
      if (!list_possible_matches.empty()) {
        char* answer = strdup(list_possible_matches.front().c_str());
        list_possible_matches.pop();
        return answer;
      }
      else
        return NULL;
    }
    else {
      for (size_t i = 0; i < commands.size(); ++i) {
        string name = commands[i];
        if (!name.compare(0, strlen(text), text))
          list_possible_matches.push(name);
      }
      return command_generator(text, 1);
    }
    return NULL;
  }

  char** gaf_completion(const char* text, int start, int /*end*/) {
    if (start == 0) {
      typedef char* (*func_ptr)(const char*, int);
      func_ptr func = command_generator;
      return rl_completion_matches(text, func);
    }
    else
      return NULL;
  }

  void Init(const vector<string>& available_commands) {
    commands = available_commands;
//    typedef char** (*func_ptr)(const char*, int, int);
//    func_ptr func = gaf_completion;
    rl_attempted_completion_function = gaf_completion/*func*/;
  }
}

