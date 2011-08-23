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

pair<int, KMerCount*> PositionRead::nextKMer( int begin ) const {
	for ( int pos = begin + 1; pos < (int)(size_-K+1); ++pos ) {
		KMerNoHashMap::const_iterator it_hash = PositionKMer::hm.find ( KMerNo(start_ + pos) );
		if ( it_hash != PositionKMer::hm.end() ) {
			return make_pair(pos, it_hash->second);
		}
	}
	return make_pair(-1, (KMerCount*)NULL);
}

