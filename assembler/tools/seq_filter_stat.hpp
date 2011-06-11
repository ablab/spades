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
#include "../test/memory.hpp"
#include <sys/time.h>

// @description Class seq_filter_stat is used for getting 
// all k-mers from all reads of a given file
// and selecting k-mers which amount is more than 1.
// @parameter size the size of k-mer
// @parameter hm hash map class
// @parameter in input file name
// @parameter name the name of hash map


template<size_t size, class hm>
class seq_filter_stat {		
public:
	static void filter(const std::string& in, 
                     const std::string& label = "") {
    double vm1 = 0;
    double rss1 = 0;
    process_mem_usage(vm1, rss1);
    timeval tim;
    gettimeofday(&tim, NULL);
    double t1 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    hm map;
		add_seqs_from_file_to_map(in, map);

    double vm2 = 0;
    double rss2 = 0;
    process_mem_usage(vm2, rss2);
    gettimeofday(&tim, NULL);
    double t2 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    size_t n = 0;
    Seq<size> seq;
    typename hm::iterator it;
    for (it = map.begin(); it != map.end(); ++it) {
      if (n % 2) {
        seq = (*it).first;
        std::string s = seq.str();
        s[size / 2] = s[size / 2 + 1];
        map.find(Seq<size>(s));
      }
      ++n;
    }

    gettimeofday(&tim, NULL);
    double t3 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    read_seqs_from_file(in);

    gettimeofday(&tim, NULL);
    double t4 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    std::cout << "Memory: " << (vm2 - vm1) << " " << label << std::endl;
    std::cout << "Insert: " << (t2 - t1) - (t4 - t3) << " " << label << std::endl;
    std::cout << "Find: " << (t3 - t2) << " " << label << std::endl;
	}

  // Only for cuckoo!!!
	static void filter(const std::string& in, 
                     const size_t d,
                     const double step,
                     const std::string& label) {
    timeval tim;
    gettimeofday(&tim, NULL);
    double t1 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    hm map(d, 100, 100, step);
    add_seqs_from_file_to_map(in, map);

    gettimeofday(&tim, NULL);
    double t2 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    size_t n = 0;
    Seq<size> seq;
    typename hm::iterator it;
    for (it = map.begin(); it != map.end(); ++it) {
      if (n % 2) {
        seq = (*it).first;
        std::string s = seq.str();
        s[size / 2] = s[size / 2 + 1];
        map.find(Seq<size>(s));
      }
      ++n;
    } 

    gettimeofday(&tim, NULL);
    double t3 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    read_seqs_from_file(in);

    gettimeofday(&tim, NULL);
    double t4 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    std::cout << "Memory: " << 1. * map.size() / map.length()  << " " << label << std::endl;
    std::cout << "Insert: " << (t2 - t1) - (t4 - t3) << " " << label << std::endl;
    std::cout << "Find: " << (t3 - t2) << " " << label << std::endl;
	}

private:
	seq_filter_stat();
	seq_filter_stat(const seq_filter_stat<size, hm>& sf);

private:
  static void add_seqs_from_file_to_map(const std::string& in, hm& map) {
    ireadstream irs(in);
    while (!irs.eof()) {
      Read r;
      irs >> r;
      add_seqs_from_read_to_map(r, map);
    }
  }

  static void read_seqs_from_file(const std::string& in) {
    ireadstream irs(in);
    while (!irs.eof()) {
      Read r;
      irs >> r;
    }
  }

	static void add_seqs_from_read_to_map(const Read& read, hm& map) {
		Sequence s = read.getSequence();
		if (s.size() >= size) {
			Seq<size> seq = s.start<size>();
			map.insert(std::make_pair(seq, 1));
			size_t s_size = s.size();
			for (size_t i = size; i < s_size; ++i) {
				Seq<size> next = seq << s[i];
				seq = next;
        map.insert(std::make_pair(seq, 1));
			}
		}
	}
};

#endif
