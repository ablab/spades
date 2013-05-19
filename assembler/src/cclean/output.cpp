#include <iostream>
#include <algorithm>
#include "output.hpp"

void print_n_times(std::ostream& output, char c, int n) {
	for(int i = 0; i < n; ++i) {
		output << c;
	}
}

void restoreFromCigar(const std::string& ref, const std::string& query, std::string& out_ref, std::string& out_query, const StripedSmithWaterman::Alignment& a) {
	std::vector<char> aligned_ref, aligned_query;
	int ref_pos = 0, query_pos = 0;
	for (std::vector<uint32_t>::const_iterator it = a.cigar.begin(); it != a.cigar.end(); ++it) {
		int num = (*it & 0xFFFFFFF0) >> 4;
		int op_code = *it & 0x0000000F;
		switch (op_code) {
			case 0: {//match
				for (int i = 0; i < num; ++i) {
					aligned_ref.push_back(ref[a.ref_begin + ref_pos++]);
					aligned_query.push_back(query[a.query_begin + query_pos++]);
				}
				break;
			}
			case 1: {//insert
				for (int i = 0; i < num; ++i) {
					aligned_ref.push_back('-');
					aligned_query.push_back(query[a.query_begin + query_pos++]);
				}
				break;
			}
			case 2: {//del
				for (int i = 0; i < num; ++i) {
					aligned_ref.push_back(ref[a.ref_begin + ref_pos++]);
					aligned_query.push_back('-');
				}
				break;
			}
			default:
				break;
		}
	}

	out_ref = std::string(aligned_ref.begin(), aligned_ref.end());
	out_query = std::string(aligned_query.begin(), aligned_query.end());
}

void print_alignment(std::ostream& output, const StripedSmithWaterman::Alignment & data,
		const std::string& ref, const std::string& query,
		const std::string& name, const std::string& database_name) {

	output << "Alignment: input sequence (first line) " << name << " alignes " << std::endl
			<< "sequence from database  (last line) " << database_name << std::endl;

	std::string aligned_query, aligned_ref;
	restoreFromCigar(ref, query, aligned_ref, aligned_query, data);

	// case when pattern's start pos is less than text one
	int text_offset = data.ref_begin - data.query_begin < 0 ? data.query_begin - data.ref_begin : 0;

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

void print_match(std::ostream& output, std::ostream& bed, std::map<std::string*, std::vector<int>, Compare>& res, const std::string& name, const std::string& seq) {
	for (std::map<std::string*, std::vector<int>, Compare>::const_iterator it = res.begin(); it != res.end(); ++it) {
		for (std::vector<int>::const_iterator it_pos = it->second.begin(); it_pos != it->second.end(); ++it_pos) {
			std::string database_name = "FIXME";
			output << "Match: input sequence (first line) " << name << " matches " << std::endl
					<< "sequence from database (2nd line) " << database_name << std::endl;

			output << seq << std::endl;
			print_n_times(output, ' ', *it_pos);
			print_n_times(output, '|', it->first->length());
			output << std::endl;
			print_n_times(output, ' ', *it_pos);
			output << *(it->first) << std::endl;
			output << std::endl;

			std::replace(database_name.begin(), database_name.end(), ' ', '_');

			print_bed(bed, name, *it_pos, *it_pos + it->first->size());
		}
	}
}

void print_bed(std::ostream& output, const std::string & name, int start, int stop) {
	output << name << "\t" << start << "\t" << stop << std::endl;
}
