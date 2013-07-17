#pragma once
#include "standard.hpp"
#include "indices/debruijn_kmer_index.hpp"
#include "runtime_k.hpp"
#include "mph_index/kmer_index.hpp"

namespace debruijn_graph {

class LinkCleaner {
private:
	typedef DeBruijnExtensionIndex<> Index;
	typedef Index::KMer Kmer;
	Index &index_;

	void CleanForwardLinks(KmerWithHash<Kmer> &kh, char i) {
		if(index_.CheckOutgoing(kh.idx, i)) {
		    KmerWithHash<Kmer> next_kh = index_.CreateKmerWithHash(kh.kmer << i);
			if(!index_.CheckIncoming(next_kh.idx, kh.kmer[0])) {
				index_.DeleteOutgoing(kh.idx, i);
			}
		}
	}

	void CleanBackwardLinks(KmerWithHash<Kmer> &kh, char i) {
		if(index_.CheckIncoming(kh.idx, i)) {
		    KmerWithHash<Kmer> prev_kh = index_.CreateKmerWithHash(kh.kmer >> i);
			if(!index_.CheckOutgoing(prev_kh.idx, kh.kmer[index_.k() - 1])) {
				index_.DeleteIncoming(kh.idx, i);
			}
		}
	}

public:
	LinkCleaner(Index &index) : index_(index) {}

	//TODO make parallel
	void CleanLinks() {
		for (auto it  = index_.kmer_begin(); it.good(); ++it) {
		    KmerWithHash<Kmer> kh = index_.CreateKmerWithHash(runtime_k::RtSeq(index_.k(), *it));
			for(char i = 0; i < 4; i++) {
				CleanForwardLinks(kh, i);
				CleanBackwardLinks(kh, i);
			}
		}
	}
};


class EarlyTipClipper {
private:
	typedef DeBruijnExtensionIndex<> Index;
	typedef Index::KMer Kmer;
	Index &index_;
	size_t length_bound_;

//Not optimal with respect to the number of large array queries (the one that contains adjacency masks). Should be ok though in case cash works the way I think it does
	size_t RemoveForward(KmerWithHash<Kmer> kh) {
        std::vector<KmerWithHash<Kmer>> tip;
		do {
			tip.push_back(kh);
			kh = index_.CreateKmerWithHash(kh.kmer << index_.GetUniqueOutgoing(kh.idx));
		} while (tip.size() < length_bound_ && index_.CheckUniqueIncoming(kh.idx) && index_.CheckUniqueOutgoing(kh.idx));

        if (!index_.CheckUniqueIncoming(kh.idx)) {
			for (size_t i = 0; i < tip.size(); i++) {
				index_.IsolateVertex(tip[i].idx);
			}
			return tip.size();
		}

		return 0;
	}

	size_t RemoveBackward(KmerWithHash<Kmer> kh) {
        std::vector<KmerWithHash<Kmer>> tip;
		do {
			tip.push_back(kh);
			kh = index_.CreateKmerWithHash(kh.kmer >> index_.GetUniqueIncoming(kh.idx));
		} while(tip.size() < length_bound_ && index_.CheckUniqueIncoming(kh.idx) && index_.CheckUniqueOutgoing(kh.idx));

        if (!index_.CheckUniqueOutgoing(kh.idx)) {
			for (size_t i = 0; i < tip.size(); i++) {
				index_.IsolateVertex(tip[i].idx);
			}
			return tip.size();
		}
		return 0;
	}

	//TODO make parallel
	size_t RoughClipTips() {
		size_t result = 0;
		for (auto it  = index_.kmer_begin(); it.good(); ++it) {
		    KmerWithHash<Kmer> kh = index_.CreateKmerWithHash(runtime_k::RtSeq(index_.k(), *it));
			if (index_.IsDeadEnd(kh.idx) && index_.CheckUniqueIncoming(kh.idx)) {
				result += RemoveBackward(kh);
			} else if(index_.IsDeadStart(kh.idx) && index_.CheckUniqueOutgoing(kh.idx)) {
				result += RemoveForward(kh);
			}
		}
		return result;
	}


public:
	EarlyTipClipper(Index &index, size_t length_bound) :
            index_(index), length_bound_(length_bound) {}

	/*
	 * Method returns the number of removed edges
	 */
	size_t ClipTips() {
		INFO("Early tip clipping");
		size_t result = RoughClipTips();
		LinkCleaner(index_).CleanLinks();
		INFO(result << " " << (index_.k()+1) <<"-mers were removed by early tip clipper");
		return result;
	}
protected:
	DECL_LOGGER("Early tip clipping");
};


class AlternativeEarlyTipClipper {
private:
	typedef DeBruijnExtensionIndex<> Index;
    typedef Index::KMer Kmer;
	Index &index_;
	size_t length_bound_;

	/*
	 * This method starts from the kmer that is second in the tip counting from junction vertex. It records all kmers of a tip into tip vector.
	 * The method returns length of a tip.
	 * In case it did not end as a tip or if it was too long tip vector is cleared and infinite length is returned.
	 * Thus tip vector contains only kmers to be removed while returned length value gives reasonable information of what happend.
	 */
	size_t FindForward(KmerWithHash<Kmer> kh, vector<KmerWithHash<Kmer>> &tip) {
		while(tip.size() < length_bound_ && index_.CheckUniqueIncoming(kh.idx) && index_.CheckUniqueOutgoing(kh.idx)) {
			tip.push_back(kh);
			kh = index_.CreateKmerWithHash(kh.kmer << index_.GetUniqueOutgoing(kh.idx));
		}
		tip.push_back(kh);
		if(index_.CheckUniqueIncoming(kh.idx) && index_.IsDeadEnd(kh.idx)) {
			return tip.size();
		}
		tip.clear();
		return -1;
	}

	size_t FindBackward(KmerWithHash<Kmer> kh, vector<KmerWithHash<Kmer>> &tip) {
		while(tip.size() < length_bound_ && index_.CheckUniqueOutgoing(kh.idx) && index_.CheckUniqueIncoming(kh.idx)) {
			tip.push_back(kh);
			kh = index_.CreateKmerWithHash(kh.kmer >> index_.GetUniqueIncoming(kh.idx));
		}
		tip.push_back(kh);
		if(index_.CheckUniqueOutgoing(kh.idx) && index_.IsDeadStart(kh.idx)) {
			return tip.size();
		}
		tip.clear();
		return -1;
	}

    size_t RemoveTip(vector<KmerWithHash<Kmer> > &tip) {
		for(size_t i = 0; i < tip.size(); i++)
			index_.IsolateVertex(tip[i].idx);
		return tip.size();
	}

    size_t RemoveTips(vector<vector<KmerWithHash<Kmer> > > tips, size_t max) {
		size_t result = 0;
		for(char c = 0; c < 4; c++) {
			if(tips[c].size() < max) {
				result += RemoveTip(tips[c]);
			}
		}
		return result;
	}

	size_t RemoveForward(KmerWithHash<Kmer> kh) {
		vector<vector<KmerWithHash<Kmer> >> tips;
		tips.resize(4);
		size_t max = 0;
		for(char c = 0; c < 4; c++) {
			if(index_.CheckOutgoing(kh.idx, c)) {
				KmerWithHash<Kmer> khc = index_.CreateKmerWithHash(kh.kmer << c);
				size_t len = FindForward(khc, tips[c]);
				if(len > max)
					max = len;
			}
		}
		return RemoveTips(tips, max);
	}

	size_t RemoveBackward(KmerWithHash<Kmer> kh) {
		vector<vector<KmerWithHash<Kmer> >> tips;
		tips.resize(4);
		size_t max = 0;
		for(char c = 0; c < 4; c++) {
			if(index_.CheckIncoming(kh.idx, c)) {
				KmerWithHash<Kmer> khc = index_.CreateKmerWithHash(kh.kmer >> c);
				size_t len = FindBackward(khc, tips[c]);
				if(len > max)
					max = len;
			}
		}
		return RemoveTips(tips, max);
	}

	//TODO make parallel
	size_t RoughClipTips() {
		size_t result = 0;
		for (auto it  = index_.kmer_begin(); it.good(); ++it) {
			KmerWithHash<Kmer> kh = index_.CreateKmerWithHash(runtime_k::RtSeq(index_.k(), *it));
			if(index_.OutgoingEdgeCount(kh.idx)  >= 2) {
				result += RemoveForward(kh);
			}
			if(index_.IncomingEdgeCount(kh.idx)  >= 2) {
				result += RemoveBackward(kh);
			}
		}
		return result;
	}


public:
	AlternativeEarlyTipClipper(Index &index, size_t length_bound) : index_(index), length_bound_(length_bound) {
	}

	/*
	 * Method returns the number of removed edges
	 */
	size_t ClipTips() {
		INFO("Early tip clipping");
		size_t result = RoughClipTips();
		LinkCleaner(index_).CleanLinks();
		INFO(result << " " << (index_.k()+1) <<"-mers were removed by early tip clipper");
		return result;
	}
protected:
	DECL_LOGGER("Early tip clipping");
};

}
