#pragma once

#include "hmmer_fwd.h"

#include <string>
#include <iostream>
#include <vector>
#include <limits>
#include <cstdlib>

std::vector<std::string> read_fasta_edges(const std::string &filename, bool add_rc = false);
