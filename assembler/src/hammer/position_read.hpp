#ifndef POSITION_READ_HPP_
#define POSITION_READ_HPP_

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <bitset>

#include "kmer_stat.hpp"
#include "globals.hpp"
#include "kmerno.hpp"

class PositionRead {
	hint_t start_;
	uint32_t size_;
	hint_t readno_;
	bool done_;
	std::bitset<K> rc_;
	
  public:
	PositionRead(hint_t start, uint32_t size, hint_t readno) : start_(start), size_(size), readno_(readno), done_(false) { rc_.reset(); }
	PositionRead(hint_t start, uint32_t size, hint_t readno, bool bad) : start_(start), size_(size), readno_(readno), done_(bad) { rc_.reset(); }
	hint_t start() const { return start_; }
	uint32_t size() const { return size_; }
	char at(uint32_t pos) const;
	char operator [] (uint32_t pos) const;
	bool isDone() { return done_; }
	void done() { done_ = true; }

	bool getRCBit(size_t i) const { return rc_[i]; }
	void setRCBit(size_t i) { rc_.set(i); }
	string getRCBitString() const { return rc_.to_string(); }

	pair<int, hint_t> nextKMerNo( int begin ) const;

	const std::string & getQualityString() const;

	hint_t getRCPosition(hint_t pos) const { return Globals::revPos + start_ + size_ - pos - K + 1; }
};

#endif

