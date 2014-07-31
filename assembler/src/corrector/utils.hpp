#pragma once
#include <vector>
#include <map>
#include <string>
#include "logger/log_writers.hpp"
#include "logger/logger.hpp"
#include "standard.hpp"

namespace corrector {
std::vector<std::string> split(const std::string &s, char delim);

std::string ContigRenameWithLength( std::string name, size_t len);
std::map<std::string, std::string> GetContigs(std::string filename);
void PutContig(std::string full_path, std::string contig_name, std::string contig_seq) ;
};
