#pragma once
#include <vector>
#include <map>
#include <string>
#include "logger/log_writers.hpp"
#include "logger/logger.hpp"
// FIXME: Get rid of this. There is nothing here which would require such junk

namespace corrector {
std::vector<std::string> split(const std::string &s, char delim);

std::string ContigRenameWithLength(std::string name, size_t len);
}
;
