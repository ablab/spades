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
		if (Globals::blob[start_ + cur_pos + K-1] == 'N') {
			cur_pos = cur_pos + K - 1;
			continue;
		}
		//cout << "    looking for " << (cur_pos) << ", direct:\t" <<  PositionKMer(start_ + cur_pos).str();
		int pos = rc_[cur_pos] ? getRCPosition(cur_pos+1) : start_ + cur_pos;
		//cout << "   rc=" << rc_[cur_pos] << "  pos = " << pos << "\t" << PositionKMer(pos).str() << endl;
		vector<hint_t>::const_iterator it_vec = lower_bound(Globals::kmernos->begin(), Globals::kmernos->end(), pos, PositionKMer::compareKMersDirect );
		//cout << "  lowerbound=" << (it_vec - Globals::kmernos->begin()) << "\tvalue=" << *it_vec;
		if ( PositionKMer::equalKMersDirect(*it_vec, pos)  ) {
			//cout << "  ok" << endl;
			return make_pair(cur_pos, (it_vec - Globals::kmernos->begin()));
		}
		//else cout << endl;
	}
	//cout << endl;
	return make_pair(-1, BLOBKMER_UNDEFINED);
}
