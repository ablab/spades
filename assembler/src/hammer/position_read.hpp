#ifndef POSITION_READ_HPP_
#define POSITION_READ_HPP_

#include <vector>
#include <string>
#include <map>
#include <fstream>

#include "kmer_stat.hpp"
#include "globals.hpp"

class PositionRead {
	hint_t start_;
	uint32_t size_;
	hint_t readno_;
	bool bad_;
	
  public:
	PositionRead(hint_t start, uint32_t size, hint_t readno) : start_(start), size_(size), readno_(readno), bad_(false) { }
	PositionRead(hint_t start, uint32_t size, hint_t readno, bool bad) : start_(start), size_(size), readno_(readno), bad_(bad) { }
	hint_t start() const { return start_; }
	uint32_t size() const { return size_; }
	char at(uint32_t pos) const;
	char operator [] (uint32_t pos) const;
	bool bad() { return bad_; }
	void setBad(bool b) { bad_ = b; }

	bool nextKMer( std::pair<uint32_t, hint_t> * it ) const;

	const std::string & getQualityString() const;
	const std::string & getName() const;
	const std::string & getSequenceString() const;
	std::string getPhredQualityString(int offset) const;

	void print(ofstream & ofs, int offset = Globals::qvoffset) const;

};

#endif

