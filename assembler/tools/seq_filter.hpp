/* 
 * seq_filter.hpp
 * 
 * Created on: 17.04.2011
 *     Author: Mariya Fomkina
 *   Modified: 15.05.2011 by author
 */

#ifndef _SEQ_FILTER_HPP_
#define _SEQ_FILTER_HPP_
 
#include "ireadstream.hpp"

// @description Class seq_filter is used for getting all 
// k-mers from all reads of a given file
// and selecting k-mers which amount is more than a given number.
// @paramater size the size of k-mer
// @parameter hm hash map class
// @parameter in input file name
// @parameter L given number for filtration

template<size_t size, class hm>
class seq_filter {		
public:
	static void filter(const std::string& in, const size_t L = 1) { 
    hm map;
    std::vector<Seq<size> > seqs;
		add_seqs_from_file_to_map(in, map);
    write_seqs_from_map_to_stdout(map, L);
	}

private:
	seq_filter();
	seq_filter(const seq_filter<size, hm>& sf);

private:
  static void add_seqs_from_file_to_map(const std::string& in, hm& map) {
    ireadstream irs(in);
    while (!irs.eof()) {
      Read r;
      irs >> r;
      add_seqs_from_read_to_map(r, map);
    }
  }

	static void add_seqs_from_read_to_map(const Read& read, hm& map) {
		Sequence s = read.getSequence();
		if (s.size() >= size) {
			Seq<size> seq = s.start<size>();
			add_seq_in_map(seq, map);
			size_t s_size = s.size();
			for (size_t i = size; i < s_size; ++i) {
				Seq<size> next = seq << s[i];
				seq = next;
				add_seq_in_map(seq, map);
			}
		}
	}

	static void add_seq_in_map(Seq<size>& seq, hm& map) {
		typename hm::iterator it = map.find(seq);
		if (it == map.end()) {
			map.insert(std::make_pair(seq, 1));
		} else {
			(*it).second++;
		}
	}

	static void write_seqs_from_map_to_stdout(hm& map, const size_t& L) {
		typename hm::iterator end = map.end();
		for (typename hm::iterator it = map.begin(); it != end; ++it) {
			if ((*it).second > L) std::cout << (*it).first << std::endl;		
		}
	}
};

#endif
