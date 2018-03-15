#pragma once

#include <string>
#include <vector>

std::vector<std::string> read_fasta_edges(const std::string &filename, bool add_rc = false);
