#include <vector>
#include "position_kmer.hpp"
#include "position_read.hpp"

char PositionRead::at(uint32_t pos) const {
	return PositionKMer::blob[ start_ + pos ];
}

char PositionRead::operator [] (uint32_t pos) const {
	return PositionKMer::blob[ start_ + pos ];
}

void PositionRead::print(ofstream & outf, int offset) const {
	outf << "@" << getName().c_str() << "\n" 
	     << PositionKMer::rv->at(readno_).getSequenceString().c_str() << "\n"
	     << "+" << getName().c_str() << "\n"
	     << getPhredQualityString( offset ).c_str() << "\n";
}

std::string PositionRead::getPhredQualityString(int offset) const {
	return PositionKMer::rv->at(readno_).getPhredQualityString(offset);
}

const std::string & PositionRead::getName() const {
	return PositionKMer::rv->at(readno_).getName();
}

const std::string & PositionRead::getSequenceString() const {
	return PositionKMer::rv->at(readno_).getSequenceString();
}

const std::string & PositionRead::getQualityString() const {
	return PositionKMer::rv->at(readno_).getQualityString();
}

bool PositionRead::nextKMer( std::pair<uint32_t, hint_t> * it ) const {
	for ( uint32_t pos = it->first + 1; pos < size_; ++pos ) {
		if (PositionKMer::blobkmers[ start_ + pos ] < BLOBKMER_UNDEFINED) {
			it->first = pos;
			it->second = PositionKMer::blobkmers[ start_ + pos ];
			return true;
		}
	}
	return false;
}


