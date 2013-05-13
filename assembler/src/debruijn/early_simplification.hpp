#pragma once
#include "standard.hpp"
#include "debruijn_kmer_index.hpp"
#include "runtime_k.hpp"
#include "mph_index/kmer_index.hpp"

namespace debruijn_graph {

class LinkCleaner {
private:
	typedef DeBruijnExtensionIndex<runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq> > Index;
	Index &index_;

//	KmerWithHash CreateKmerWithHash(runtime_k::RtSeq kmer) const {
//		return KmerWithHash(kmer, index_);
//	}

	void CleanForwardLinks(Index::KmerWithHash &kh, char i) {
		if(index_.CheckOutgoing(kh.idx, i)) {
			Index::KmerWithHash next_kh = index_.CreateKmerWithHash(kh.kmer << i);
			if(index_.CheckIncoming(next_kh.idx, kh.kmer[0])) {
				index_.DeleteOutgoing(kh.idx, i);
			}
		}
	}

	void CleanBackwardLinks(Index::KmerWithHash &kh, char i) {
		if(index_.CheckIncoming(kh.idx, i)) {
			Index::KmerWithHash prev_kh = index_.CreateKmerWithHash(kh.kmer >> i);
			if(index_.CheckOutgoing(prev_kh.idx, kh.kmer[index_.K() - 1])) {
				index_.DeleteIncoming(kh.idx, i);
			}
		}
	}

public:
	LinkCleaner(DeBruijnExtensionIndex<runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq> > &index) : index_(index) {
	}

	//TODO make parallel
	void CleanLinks() {
		for (auto it  = index_.kmer_begin(), et = index_.kmer_end(); it != et; ++it) {
			Index::KmerWithHash kh(runtime_k::RtSeq(index_.K(), (*it).data()), index_);
			for(char i = 0; i < 4; i++) {
				CleanForwardLinks(kh, i);
				CleanBackwardLinks(kh, i);
			}
		}
	}
};


class EarlyTipClipper {
private:
	typedef DeBruijnExtensionIndex<runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq> > Index;
	Index &index_;
	size_t length_bound_;

//Not optimal with respect to the number of large array queries (the one that contains adjacency masks). Should be ok though in case cash works the way I think it does
	size_t RemoveForward(Index::KmerWithHash kh) {
		vector<Index::KmerWithHash> tip;
		do {
			tip.push_back(kh);
			kh = index_.CreateKmerWithHash(kh.kmer << index_.GetUniqueOutgoing(kh.idx));
		} while(tip.size() < length_bound_ && index_.CheckUniqueIncoming(kh.idx) && index_.CheckUniqueOutgoing(kh.idx));
		if(!index_.CheckUniqueOutgoing(kh.idx)) {
			for(size_t i = 0; i < tip.size(); i++) {
				index_.DeleteOutgoing(tip[i].idx, index_.GetUniqueOutgoing(kh.idx));
			}
			for(size_t i = 1; i < tip.size(); i++) {
				index_.DeleteIncoming(tip[i].idx, index_.GetUniqueIncoming(kh.idx));
			}
			return tip.size();
		}
		return 0;
	}

	size_t RemoveBackward(Index::KmerWithHash kh) {
		vector<Index::KmerWithHash> tip;
		do {
			tip.push_back(kh);
			kh = index_.CreateKmerWithHash(kh.kmer >> index_.GetUniqueIncoming(kh.idx));
		} while(tip.size() < length_bound_ && index_.CheckUniqueIncoming(kh.idx) && index_.CheckUniqueOutgoing(kh.idx));
		if(!index_.CheckUniqueIncoming(kh.idx)) {
			for(size_t i = 0; i < tip.size(); i++) {
				index_.DeleteIncoming(tip[i].idx, index_.GetUniqueIncoming(kh.idx));
			}
			for(size_t i = 1; i < tip.size(); i++) {
				index_.DeleteOutgoing(tip[i].idx, index_.GetUniqueOutgoing(kh.idx));
			}
			return tip.size();
		}
		return 0;
	}

	//TODO make parallel
	size_t RoughClipTips() {
		size_t result = 0;
		for (auto it  = index_.kmer_begin(), et = index_.kmer_end(); it != et; ++it) {
			Index::KmerWithHash kh = index_.CreateKmerWithHash(runtime_k::RtSeq(index_.K(), (*it).data()));
			if(index_.IsDeadEnd(kh.idx) && index_.CheckUniqueIncoming(kh.idx)) {
				result += RemoveBackward(kh);
			} else if(index_.IsDeadStart(kh.idx) && index_.CheckUniqueOutgoing(kh.idx)) {
				result += RemoveForward(kh);
			}
		}
		return result;
	}


public:
	EarlyTipClipper(DeBruijnExtensionIndex<runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq> > &index, size_t length_bound) : index_(index), length_bound_(length_bound) {
	}

	/*
	 * Method returns the number of removed edges
	 */
	size_t ClipTips() {
		size_t result = RoughClipTips();
		LinkCleaner(index_).CleanLinks();
		return result;
	}
};

}
