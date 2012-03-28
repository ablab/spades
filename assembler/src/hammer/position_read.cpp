#include "standard.hpp"
#include "position_kmer.hpp"
#include "position_read.hpp"
#include "config_struct_hammer.hpp"

char PositionRead::at(uint32_t pos) const {
	return Globals::blob[ start_ + pos ];
}

char PositionRead::operator [] (uint32_t pos) const {
	return Globals::blob[ start_ + pos ];
}

pair<int, hint_t> PositionRead::nextKMerNo( int begin ) const {
	for ( int pos = begin + 1; pos < (int)(size_-K+1); ++pos ) {
		//cout << "    looking for " << (start_ + pos) << ": " << PositionKMer(start_ + pos).str();
		vector<hint_t>::const_iterator it_vec = lower_bound(Globals::kmernos->begin(), Globals::kmernos->end(), start_ + pos, PositionKMer::compareKMersDirect );
		//cout << "  lowerbound=" << (it_vec - Globals::kmernos->begin()) << "\tvalue=" << *it_vec;
		if ( it_vec != Globals::kmernos->end() && PositionKMer::equalKMersDirect(*it_vec, start_ + pos)  ) {
			//cout << "  ok" << endl;
			return make_pair(pos, (it_vec - Globals::kmernos->begin()));
		}
	}
	//cout << endl;
	return make_pair(-1, BLOBKMER_UNDEFINED);
}
