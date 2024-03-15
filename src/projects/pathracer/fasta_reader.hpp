#pragma once

#include <string>
#include <vector>
#include <utility>

std::string rev_comp(const std::string &s);
std::vector<std::string> read_fasta_edges(const std::string &filename, bool add_rc = false);
std::vector<std::pair<std::string, std::string>> read_fasta(const std::string &filename);
