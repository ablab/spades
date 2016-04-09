//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <iostream>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "output.hpp"
#include "utils.hpp"

namespace cclean_output {

void print_n_times(std::ostream& output, char c, int n) {
  for (int i = 0; i < n; ++i) {
    output << c;
  }
}

void print_alignment(std::ostream& output, const StripedSmithWaterman::Alignment &data,
        const std::string& ref, const std::string& query,
        const std::string& name, const std::string& database_name) {

  output << "Alignment: input sequence (first line) " << name << " alignes "
         << std::endl
         << "sequence from database (last line) " << database_name << std::endl;

  std::string aligned_query, aligned_ref;
  cclean_utils::RestoreFromCigar(ref, query, aligned_ref, aligned_query, data);

  // case when pattern's start pos is less than text one
  int text_offset = data.ref_begin - data.query_begin < 0 ? data.query_begin
                                                            - data.ref_begin : 0;

  // ref = read
  print_n_times(output, ' ', text_offset);
  output << ref << std::endl;
  print_n_times(output, ' ', text_offset + data.ref_begin);
  output << aligned_ref << std::endl;

  // vertical dashes
  print_n_times(output, ' ', text_offset + data.ref_begin);
  for (int i = 0; i < (int)std::min(aligned_query.length(), aligned_ref.length()); ++i) {
   aligned_query.at(i) == aligned_ref.at(i) ? output << "|" : output << "*";
  }
  output << std::endl;

  // query = contamination
  print_n_times(output, ' ', text_offset + data.ref_begin);
  output << aligned_query << std::endl;
  print_n_times(output, ' ', data.ref_begin - data.query_begin);
  output << query << std::endl;
  output << std::endl;
 }

void print_match(std::ostream& output, std::ostream& bed, std::map<std::string*,
                  std::vector<int>, Compare>& res, const std::string& name,
                  const std::string& seq, const std::string &db_name) {
  for (auto it = res.begin(); it != res.end(); ++it) {
   for (auto it_pos = it->second.begin(); it_pos != it->second.end(); ++it_pos) {

    output << "Match: input sequence (first line) " << name << " matches "
           << std::endl
           << "sequence from database (2nd line) " << db_name << std::endl;

    output << seq << std::endl;
    print_n_times(output, ' ', *it_pos);
    print_n_times(output, '|', it->first->length());
    output << std::endl;
    print_n_times(output, ' ', *it_pos);
    output << *(it->first) << std::endl;
    output << std::endl;

    print_bad(bed, name, *it_pos, *it_pos + it->first->size());
   }
  }
}
//end of namespace
}
