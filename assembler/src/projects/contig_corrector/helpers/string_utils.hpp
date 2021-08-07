#pragma once

#include <string>

namespace helpers {

inline std::string ToUpperStr(std::string s) {
    for (auto & c : s)
        c = static_cast<char>(toupper(c));
    return s;
}

inline std::string ToLowerStr(std::string s) {
    for (auto & c : s)
        c = static_cast<char>(tolower(c));
    return s;
}

} // namespace helpers
