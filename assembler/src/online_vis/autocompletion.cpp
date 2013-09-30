#include <vector>
#include <string>
#include <queue>
#include <readline/readline.h>

namespace online_visualization {

std::vector<std::string> commands;

char* CommandGenerator(const char* text, int state) {
    static std::queue<std::string> list_possible_matches;

    if (state != 0) {
        if (!list_possible_matches.empty()) {
            char* answer = strdup(list_possible_matches.front().c_str());
            list_possible_matches.pop();
            return answer;
        } else
            return NULL;
    } else {
        for (size_t i = 0; i < commands.size(); ++i) {
            std::string name = commands[i];
            if (!name.compare(0, strlen(text), text))
                list_possible_matches.push(name);
        }
        return CommandGenerator(text, 1);
    }
    return NULL;
}

char** GafCompletion(const char* text, int start, int /*end*/) {
    if (start == 0) {
        return rl_completion_matches(text, CommandGenerator);
    } else
        return NULL;
}

void InitAutocompletion(const std::vector<std::string>& available_commands) {
    commands = available_commands;
    rl_attempted_completion_function = GafCompletion;
}

}

