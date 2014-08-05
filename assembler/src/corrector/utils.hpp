#pragma once
#include <vector>
#include <map>
#include <string>
#include "logger/log_writers.hpp"
#include "logger/logger.hpp"

namespace corrector {
std::vector<std::string> split(const std::string &s, char delim);

std::string ContigRenameWithLength(std::string name, size_t len);
}
;
