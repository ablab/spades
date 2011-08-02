#ifndef POSITION_READ_HPP_
#define POSITION_READ_HPP_

#include <vector>
#include <string>
#include <map>
#include <fstream>

#include "kmer_stat.hpp"

class PositionRead {
	hint_t start_;
	uint32_t size_;
	hint_t readno_;
	
  public:
	PositionRead(hint_t start, uint32_t size, hint_t readno) : start_(start), size_(size), readno_(readno) { }
	hint_t start() const { return start_; }
	uint32_t size() const { return size_; }
	char at(uint32_t pos) const;
	char operator [] (uint32_t pos) const;
	// std::map<uint32_t, hint_t> & kmers() { return kmers_; }
	// void clearKMers() { kmers_.clear(); }

	bool nextKMer( std::pair<uint32_t, hint_t> * it ) const;

	const std::string & getQualityString() const;
	const std::string & getName() const;
	const std::string & getSequenceString() const;
	std::string getPhredQualityString(int offset) const;

	void print(ofstream & ofs, int offset) const;

};

#endif

