#pragma once
#include <vector>
#include <map>
#include <string>

#include "path_helper.hpp"
namespace corrector {
std::vector<std::string> split(const std::string &s, char delim);

std::string ContigRenameWithLength(std::string name, size_t len);
}
;
