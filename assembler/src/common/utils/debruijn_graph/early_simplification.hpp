//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "utils/standard_base.hpp"
#include "utils/indices/perfect_hash_map.hpp"
#include "utils/mph_index/kmer_index.hpp"

namespace debruijn_graph {

class LinkCleaner {
private:
    typedef DeBruijnExtensionIndex<> Index;
    typedef Index::KMer Kmer;
    typedef Index::KeyWithHash KeyWithHash;
    Index &index_;

    void CleanForwardLinks(KeyWithHash &kh, char i) {
        if(index_.CheckOutgoing(kh, i)) {
            KeyWithHash next_kh = index_.GetOutgoing(kh, i);
            if(!index_.CheckIncoming(next_kh, kh[0])) {
                index_.DeleteOutgoing(kh, i);
            }
        }
    }

    void CleanBackwardLinks(KeyWithHash &kh, char i) {
        if(index_.CheckIncoming(kh, i)) {
            KeyWithHash prev_kh = index_.GetIncoming(kh, i);
            if(!index_.CheckOutgoing(prev_kh, kh[index_.k() - 1])) {
                index_.DeleteIncoming(kh, i);
            }
        }
    }

public:
    LinkCleaner(Index &index) : index_(index) {}

    //TODO make parallel
    void CleanLinks() {
        vector<Index::kmer_iterator> iters = index_.kmer_begin(10 * omp_get_max_threads());
#   pragma omp parallel for schedule(guided)
        for(size_t i = 0; i < iters.size(); i++) {
            for (Index::kmer_iterator &it = iters[i]; it.good(); ++it) {
                KeyWithHash kh = index_.ConstructKWH(RtSeq(index_.k(), *it));
                if (kh.is_minimal()) {
                    KeyWithHash kh = index_.ConstructKWH(RtSeq(index_.k(), *it));
                    for (char i = 0; i < 4; i++) {
                        CleanForwardLinks(kh, i);
                        CleanBackwardLinks(kh, i);
                    }
                }
            }
        }
    }
};

class AlternativeEarlyTipClipper {
private:
    typedef DeBruijnExtensionIndex<> Index;
    typedef Index::KMer Kmer;
    typedef Index::KeyWithHash KeyWithHash;
    Index &index_;
    size_t length_bound_;

    /*
     * This method starts from the kmer that is second in the tip counting from junction vertex. It records all kmers of a tip into tip vector.
     * The method returns length of a tip.
     * In case it did not end as a tip or if it was too long tip vector is cleared and infinite length is returned.
     * Thus tip vector contains only kmers to be removed while returned length value gives reasonable information of what happend.
     */
    size_t FindForward(KeyWithHash kh, vector<KeyWithHash> &tip) {
        while(tip.size() < length_bound_ && index_.CheckUniqueIncoming(kh) && index_.CheckUniqueOutgoing(kh)) {
            tip.push_back(kh);
            kh = index_.GetUniqueOutgoing(kh);
        }
        tip.push_back(kh);
        if(index_.CheckUniqueIncoming(kh) && index_.IsDeadEnd(kh)) {
            return tip.size();
        }
        tip.clear();
        return -1;
    }

    size_t FindBackward(KeyWithHash kh, vector<KeyWithHash> &tip) {
        while(tip.size() < length_bound_ && index_.CheckUniqueOutgoing(kh) && index_.CheckUniqueIncoming(kh)) {
            tip.push_back(kh);
            kh = index_.GetUniqueIncoming(kh);
        }
        tip.push_back(kh);
        if(index_.CheckUniqueOutgoing(kh) && index_.IsDeadStart(kh)) {
            return tip.size();
        }
        tip.clear();
        return -1;
    }

    size_t RemoveTip(vector<KeyWithHash > &tip) {
        for(size_t i = 0; i < tip.size(); i++)
            index_.IsolateVertex(tip[i]);
        return tip.size();
    }

    size_t RemoveTips(vector<vector<KeyWithHash > > tips, size_t max) {
        size_t result = 0;
        for(char c = 0; c < 4; c++) {
            if(tips[c].size() < max) {
                result += RemoveTip(tips[c]);
            }
        }
        return result;
    }

    size_t RemoveForward(KeyWithHash kh) {
        vector<vector<KeyWithHash >> tips;
        tips.resize(4);
        size_t max = 0;
        for(char c = 0; c < 4; c++) {
            if(index_.CheckOutgoing(kh, c)) {
                KeyWithHash khc = index_.GetOutgoing(kh, c);
                size_t len = FindForward(khc, tips[c]);
                if(len > max)
                    max = len;
            }
        }
        return RemoveTips(tips, max);
    }

    size_t RemoveBackward(KeyWithHash kh) {
        vector<vector<KeyWithHash >> tips;
        tips.resize(4);
        size_t max = 0;
        for(char c = 0; c < 4; c++) {
            if(index_.CheckIncoming(kh, c)) {
                KeyWithHash khc = index_.GetIncoming(kh, c);
                size_t len = FindBackward(khc, tips[c]);
                if(len > max)
                    max = len;
            }
        }
        return RemoveTips(tips, max);
    }

    //TODO make parallel
    size_t RoughClipTips() {
        vector<Index::kmer_iterator> iters = index_.kmer_begin(10 * omp_get_max_threads());
        vector<size_t> result(iters.size());
#   pragma omp parallel for schedule(guided)
        for(size_t i = 0; i < iters.size(); i++) {
            for(Index::kmer_iterator &it = iters[i]; it.good(); ++it) {
                KeyWithHash kh = index_.ConstructKWH(RtSeq(index_.k(), *it));
                if(kh.is_minimal()) {
                    if (index_.OutgoingEdgeCount(kh) >= 2) {
                        result[i] += RemoveForward(kh);
                    }
                    if (index_.IncomingEdgeCount(kh) >= 2) {
                        result[i] += RemoveBackward(kh);
                    }
                }
            }
        }
        size_t sum = 0;
        for(size_t i = 0; i < result.size(); i++)
            sum += result[i];
        return sum;
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
