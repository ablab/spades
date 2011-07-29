#ifndef POSITION_READ_HPP_
#define POSITION_READ_HPP_

#include <vector>
#include <string>
#include <map>
#include <fstream>

class PositionRead {
	uint64_t start_;
	uint32_t size_;
	uint64_t readno_;
	
	// std::map<uint32_t, uint64_t> kmers_;

  public:
	PositionRead(uint64_t start, uint32_t size, uint64_t readno) : start_(start), size_(size), readno_(readno) { }
	uint64_t start() const { return start_; }
	uint32_t size() const { return size_; }
	char at(uint32_t pos) const;
	char operator [] (uint32_t pos) const;
	// std::map<uint32_t, uint64_t> & kmers() { return kmers_; }
	// void clearKMers() { kmers_.clear(); }

	bool nextKMer( std::pair<uint32_t, uint64_t> * it ) const;

	const std::string & getQualityString() const;
	const std::string & getName() const;
	const std::string & getSequenceString() const;
	std::string getPhredQualityString(int offset) const;

	void print(ofstream & ofs, int offset) const;

};

#endif

