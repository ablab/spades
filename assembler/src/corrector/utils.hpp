#pragma once
#include <vector>
#include <string>
namespace corrector {
std::vector<std::string> split(const std::string &s, char delim);

std::string ContigRenameWithLength( std::string name, size_t len);

//inline static bool IsValidVariant(char C);
};
