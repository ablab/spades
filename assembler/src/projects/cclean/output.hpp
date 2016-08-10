//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <string>
#include <vector>
#include <map>
#include <io/read.hpp>
#include <ostream>
#include "comparator.hpp"
#include "modules/ssw_cpp.h"

namespace cclean_output {

void print_n_times(std::ostream& output, char c, int n);

void print_alignment(std::ostream& output,
                     const StripedSmithWaterman::Alignment & data,
                     const std::string& ref,
                     const std::string& query, const std::string& name,
                     const std::string& database_name);

void print_match(std::ostream& output, std::ostream& bed, std::map<std::string*,
                 std::vector<int>, Compare>& res, const std::string& name,
                 const std::string& seq, const std::string &db_name);

void print_bad(std::ostream& output, const std::string & name,
               int start, int stop);

inline void print_read(std::ostream& output, const Read &read) {
    std::ofstream &stream =
    reinterpret_cast<std::ofstream&>(output);
    read.print(stream, Read::PHRED_OFFSET);
}

inline void print_bad(std::ostream& output, const std::string & name,
                      int start, int stop) {
         output << name << "\t" << start << "\t" << stop << std::endl;
}

// end of namespace
}
#endif /* OUTPUT_H_ */
