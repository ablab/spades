#ifndef POSITION_READ_HPP_
#define POSITION_READ_HPP_

#include <vector>
#include <string>
#include <map>
#include <fstream>

#include "kmer_stat.hpp"
#include "globals.hpp"
#include "kmerno.hpp"

class PositionRead {
	hint_t start_;
	uint32_t size_;
	hint_t readno_;
	bool done_;
	
  public:
	PositionRead(hint_t start, uint32_t size, hint_t readno) : start_(start), size_(size), readno_(readno), done_(false) { }
	PositionRead(hint_t start, uint32_t size, hint_t readno, bool bad) : start_(start), size_(size), readno_(readno), done_(bad) { }
	hint_t start() const { return start_; }
	uint32_t size() const { return size_; }
	char at(uint32_t pos) const;
	char operator [] (uint32_t pos) const;
	bool isDone() { return done_; }
	void done() { done_ = true; }

	pair<int, hint_t> nextKMerNo( int begin ) const;

	const std::string & getQualityString() const;

};

#endif

