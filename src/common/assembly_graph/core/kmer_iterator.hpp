//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "graph.hpp"

#include "sequence/sequence.hpp"
#include "sequence/rtseq.hpp"

#include "adt/iterator_range.hpp"

#include <boost/iterator/iterator_facade.hpp>

namespace debruijn_graph {

class kmer_iterator : public boost::iterator_facade<kmer_iterator,
                                                    RtSeq,
                                                    boost::forward_traversal_tag,
                                                    const RtSeq&> {
  public:

    explicit kmer_iterator(size_t k = 0)
            : k_(k), pos_(-1ULL) {}

    kmer_iterator(const Sequence &seq, size_t k)
            : kmer_iterator(k) {
        reset(seq);
    }

    bool good() const { return pos_ != -1ULL; }

    void reset(const Sequence &seq) {
        seq_ = seq;
        if (seq_.size() < k_) {
            pos_ = -1ULL;
            return;
        }

        kmer_ = seq.start<RtSeq>(k_);
        pos_ = k_ - 1;
    }
    
  private:
    friend class boost::iterator_core_access;
    
    void increment() {
        pos_ += 1;
        if (pos_ >= seq_.size()) {
            pos_ = -1ULL;
            return;
        }
        kmer_ <<= seq_[pos_];
    }

    bool equal(const kmer_iterator &other) const {
        if (pos_ != other.pos_)
            return false;
        if (pos_ == -1ULL && other.pos_ == -1ULL) // end iterators always equal
            return true;
        
        return (k_ == other.k_ && seq_ == other.seq_);
    }
    
    const RtSeq &dereference() const {
        return kmer_;
    }

    Sequence seq_;
    size_t k_;
    size_t pos_;
    RtSeq kmer_;
};

// We expect that It dereferences to EdgeId
template<class It>
class kmer_graph_iterator : public boost::iterator_facade<kmer_graph_iterator<It>,
                                                          RtSeq,
                                                          boost::forward_traversal_tag,
                                                          const RtSeq&> {
  public:
    kmer_graph_iterator() {}
    
    kmer_graph_iterator(It begin, It end, size_t k,
                        const Graph &g)
            : begin_(std::move(begin)), end_(std::move(end)),
              inner_iterator_(k), g_(g) {
        if (begin_ == end_)
            return;
        
        inner_iterator_.reset(g_.get().EdgeNucls(*begin_));
    }
    
  private:
    friend class boost::iterator_core_access;

    bool good() const {
        return begin_ != end_;
    }
    
    void increment() {
        ++inner_iterator_;
        
        while (!inner_iterator_.good() &&
               ++begin_ != end_) {
            inner_iterator_.reset(g_.get().EdgeNucls(*begin_));
        }
    }

    bool equal(const kmer_graph_iterator &other) const {
        if (!good() && !other.good())
            return true;
        
        return begin_ == other.begin_ && inner_iterator_ == other.inner_iterator_;
    }
    
    const RtSeq &dereference() const {
        return *inner_iterator_;
    }

    It begin_, end_;
    kmer_iterator inner_iterator_;
    std::reference_wrapper<const Graph> g_;
};

template<class It>
auto kmer_begin(It begin, It end, const Graph &g, size_t k = 0) {
    return kmer_graph_iterator<It>(std::move(begin), std::move(end),
                                   k == 0 ? g.k() : k, g);
}

template<class It>
auto kmer_end(It end, const Graph &g) {
    return kmer_graph_iterator<It>(end, end, 0, g);
}

inline auto kmer_begin(const Graph &g, size_t k = 0) {
    return kmer_begin(g.e_begin(), g.e_end(), g, k);
}

inline auto kp1mer_begin(const Graph &g) {
    return kmer_begin(g.e_begin(), g.e_end(), g, g.k() + 1);
}

inline auto kmer_end(const Graph &g) {
    return kmer_end(g.e_end(), g);
}

inline auto kp1mer_end(const Graph &g) {
    return kmer_end(g.e_end(), g);
}

inline auto kmer_begin(EdgeId e, const Graph &g, size_t k = 0) {
    return kmer_iterator(g.EdgeNucls(e), k == 0 ? g.k() : k);
}

inline auto kmer_end(EdgeId, const Graph &) {
    return kmer_iterator();
}

inline auto kp1mer_begin(EdgeId e, const Graph &g) {
    return kmer_iterator(g.EdgeNucls(e), g.k() + 1);
}

inline auto kp1mer_end(EdgeId, const Graph &) {
    return kmer_iterator();
}

}
