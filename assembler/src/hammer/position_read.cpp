#include <vector>
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
	for ( int cur_pos = begin + 1; cur_pos < (int)(size_-K+1); ++cur_pos ) {
		int pos = rc_[cur_pos] ? getRCPosition(cur_pos+1) : start_ + cur_pos;
		cout << "    looking for " << (pos) << " (from " << (cur_pos) << "): " << PositionKMer(pos).str() << endl;
		cout << "         direct:\t" <<  PositionKMer(start_ + cur_pos).str();
		vector<hint_t>::const_iterator it_vec = lower_bound(Globals::kmernos->begin(), Globals::kmernos->end(), pos, PositionKMer::compareKMersDirect );
		cout << "  lowerbound=" << (it_vec - Globals::kmernos->begin()) << "\tvalue=" << *it_vec;
		if ( PositionKMer::equalKMersDirect(*it_vec, pos)  ) {
			cout << "  ok" << endl;
			return make_pair(cur_pos, (it_vec - Globals::kmernos->begin()));
		}
		else cout << endl;
	}
	cout << endl;
	return make_pair(-1, BLOBKMER_UNDEFINED);
}
