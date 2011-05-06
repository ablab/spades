/* 
 * seq_filter.hpp
 * 
 * Created on: 17.04.2011
 *     Author: Mariya Fomkina
 *   Modified: 25.04.2011 by author
 */

#ifndef _SEQ_FILTER_HPP_
#define _SEQ_FILTER_HPP_
 
#include "ireadstream.hpp"
#include "../test/memory.hpp"
#include <sys/time.h>

// Class seq_filter is used for getting all k-mer from all read of a given file
// and selecting k-mers which amount is more than a given number.
// Template parameters are:
// size_t size - the size of k-mer
// hm - hash map class
// The only public method filter takes 4 parameters:
// const std::string& in - input file name
// const size_t L - given number for filtration
// const bool stat - if false, k-mers are printed to output, 
// else stat info in printed
// const bool console - if true, k-mers are printed to output, 
// else k-mers are added to vector and stat option is ignored

template<size_t size, class hm>
class seq_filter {		
public:
	static std::vector<Seq<size> > filter(const std::string& in, 
                                        const size_t L = 1, 
                                        const bool stat = false, 
                                        const bool console = true,
                                        const bool find = false) {
    double vm1 = 0;
    double rss1 = 0;
    process_mem_usage(vm1, rss1);
    timeval tim;
    gettimeofday(&tim, NULL);
    double t1 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    hm map;
    std::vector<Seq<size> > seqs;
		add_seqs_from_file_to_map(in, map);
    if (console) {
      write_seqs_from_map_to_stdout(map, L, stat);
    } else {
      seqs = get_seqs_from_map(map, L);
    }

    double vm2 = 0;
    double rss2 = 0;
    process_mem_usage(vm2, rss2);
    gettimeofday(&tim, NULL);
    double t2 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    if (find) {
      Seq<size> seq;
      typename hm::iterator it;
      typename hm::iterator it2;
      for (it = map.begin(); it != map.end(); ++it) {
        seq = (*it).first;
        /*it2 = */map.find(seq);
      }
    }

    gettimeofday(&tim, NULL);
    double t3 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    if ((stat) && (console)) {
      std::cout << "Memory: " << (vm2 - vm1) << std::endl;
      std::cout << "Insert: " << (t2 - t1) << std::endl;
      std::cout << "Find: " << (t3 - t2) << std::endl;
    }
    return seqs;
	}

  // Only for cuckoo!!!
  // Code is dublicated, needs refactoring
	static std::vector<Seq<size> > filter(const std::string& in, 
                                        const size_t L, 
                                        const bool stat, 
                                        const bool console,
                                        const bool find,
                                        const size_t max_loop) {
    double vm1 = 0;
    double rss1 = 0;
    process_mem_usage(vm1, rss1);
    timeval tim;
    gettimeofday(&tim, NULL);
    double t1 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    hm map;
    map.set_up(4, 100, max_loop, 1.2);

    std::vector<Seq<size> > seqs;
		add_seqs_from_file_to_map(in, map);
    if (console) {
      write_seqs_from_map_to_stdout(map, L, stat);
    } else {
      seqs = get_seqs_from_map(map, L);
    }

    double vm2 = 0;
    double rss2 = 0;
    process_mem_usage(vm2, rss2);
    gettimeofday(&tim, NULL);
    double t2 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    if (find) {
      Seq<size> seq;
      typename hm::iterator it;
      for (it = map.begin(); it != map.end(); ++it) {
        seq = (*it).first;
        map.find(seq);
      }
    }

    gettimeofday(&tim, NULL);
    double t3 = tim.tv_sec + ((float)tim.tv_usec/1e6);

    if ((stat) && (console)) {
      std::cout << "Memory: " << (vm2 - vm1) << std::endl;
      std::cout << "Insert: " << (t2 - t1) << std::endl;
      std::cout << "Find: " << (t3 - t2) << std::endl;
    }
    return seqs;
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

	static std::vector<Seq<size> > get_seqs_from_map(hm& map, const size_t& L) {
		std::vector<Seq<size> > seqs;
		typename hm::iterator end = map.end();
		for (typename hm::iterator it = map.begin(); it != end; ++it) {
			if ((*it).second > L) seqs.push_back((*it).first);
		}
		return seqs;
	}

	static void write_seqs_from_map_to_stdout(hm& map, const size_t& L, bool stat) {
		typename hm::iterator end = map.end();
		size_t cnt = 0;
		for (typename hm::iterator it = map.begin(); it != end; ++it) {
			if ((*it).second > L) {
        if (stat) {
          ++cnt;
        } else {
          std::cout << (*it).first << std::endl;
				}
			}
		}
    if (stat) {
      std::cout << "Selected " << cnt << " k-mers from " << map.size() 
                << ", removed "  << ((float)map.size() - cnt)/map.size()*100 
                << "% of k-mers." << std::endl; 
    }
	}
};

#endif
