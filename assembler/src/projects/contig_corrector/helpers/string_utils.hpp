#pragma once

#include <string>

inline std::string ToUpperStr(std::string s) {
    for (auto & c : s)
        c = toupper(c);
    return s;
}

inline std::string ToLowerStr(std::string s) {
    for (auto & c : s)
        c = tolower(c);
    return s;
}
