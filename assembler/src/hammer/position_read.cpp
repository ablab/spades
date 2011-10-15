#include <vector>
#include "position_kmer.hpp"
#include "position_read.hpp"

char PositionRead::at(uint32_t pos) const {
	return Globals::blob[ start_ + pos ];
}

char PositionRead::operator [] (uint32_t pos) const {
	return Globals::blob[ start_ + pos ];
}

void PositionRead::print(ofstream & outf, int offset) const {
	outf << "@" << Globals::rv->at(readno_).getName().c_str() << "\n";
	for (int i=0; i < Globals::rv->at(readno_).ltrim(); ++i) outf << "N";
	outf << Globals::rv->at(readno_).getSequenceString().c_str();
	for (int i=0; i < Globals::rv->at(readno_).initial_size() - Globals::rv->at(readno_).rtrim(); ++i) outf << "N";
	outf << "\n" << "+" << Globals::rv->at(readno_).getName().c_str();
	if (Globals::rv->at(readno_).ltrim() > 0) outf << " ltrim=" << Globals::rv->at(readno_).ltrim();
	if (Globals::rv->at(readno_).rtrim() < Globals::rv->at(readno_).initial_size())
		outf << " rtrim=" << (Globals::rv->at(readno_).initial_size() - Globals::rv->at(readno_).rtrim());
	outf << "\n";
	char badq = (char)( offset + 2 );
	for (int i=0; i < Globals::rv->at(readno_).ltrim(); ++i) outf << badq;
	outf << Globals::rv->at(readno_).getPhredQualityString( offset ).c_str();
	for (int i=0; i < Globals::rv->at(readno_).initial_size() - Globals::rv->at(readno_).rtrim(); ++i) outf << badq;
	outf << "\n";
}

std::string PositionRead::getPhredQualityString(int offset) const {
	return Globals::rv->at(readno_).getPhredQualityString(offset);
}

const std::string & PositionRead::getName() const {
	return Globals::rv->at(readno_).getName();
}

const std::string & PositionRead::getSequenceString() const {
	return Globals::rv->at(readno_).getSequenceString();
}

const std::string & PositionRead::getQualityString() const {
	return Globals::rv->at(readno_).getQualityString();
}

pair<int, KMerCount*> PositionRead::nextKMer( int begin ) const {
	for ( int pos = begin + 1; pos < (int)(size_-K+1); ++pos ) {
		// cout << "    looking for " << (start_ + pos) << ": " << PositionKMer(start_ + pos).str() << endl;
		KMerNoHashMap::const_iterator it_hash = Globals::hm.find ( KMerNo(start_ + pos) );
		if ( it_hash != Globals::hm.end() ) {
			return make_pair(pos, it_hash->second);
		}
	}
	return make_pair(-1, (KMerCount*)NULL);
}

pair<int, hint_t> PositionRead::nextKMerNo( int begin ) const {
	for ( int pos = begin + 1; pos < (int)(size_-K+1); ++pos ) {
		// cout << "    looking for " << (start_ + pos) << ": " << PositionKMer(start_ + pos).str() << endl;
		vector<hint_t>::const_iterator it_vec = lower_bound(Globals::kmernos->begin(), Globals::kmernos->end(), start_ + pos, PositionKMer::compareKMersDirect );
		if ( *it_vec == start_ + pos ) return make_pair(pos, *it_vec);
	}
	return make_pair(-1, BLOBKMER_UNDEFINED);
}
