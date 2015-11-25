#pragma once

#include <regex>
#include "io/single_read.hpp"

namespace debruijn_graph {

typedef std::string bin_id;
typedef std::string contig_id;
typedef std::pair<contig_id, std::vector<bin_id>> ContigAnnotation;

inline contig_id GetId(const io::SingleRead& contig) {
     std::smatch m;
     std::regex e ("ID_(\\d+)$");
     bool success = std::regex_search(contig.name(), m, e);
     VERIFY(success);
     return m[1];
}

}
