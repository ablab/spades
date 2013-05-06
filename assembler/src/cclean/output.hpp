#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <string>
#include <vector>
#include <map>
#include <ostream>
#include "comparator.hpp"
#include "database.hpp"
#include "ssw/ssw_cpp.h"

void print_n_times(std::ostream& output, char c, int n);
void print_alignment(std::ostream& output, const StripedSmithWaterman::Alignment & data,
		const std::string& ref, const std::string& query, const std::string& name,
		const std::string& database_name, const std::string& database_comment);
void print_match(std::ostream& output, std::ostream& bed, std::map<std::string*, std::vector<int>, Compare>& res, const std::string& name, const std::string& seq, const Database * data);
void print_bed(std::ostream& output, const std::string & name, int start, int stop);

#endif /* OUTPUT_H_ */
