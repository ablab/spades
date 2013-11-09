#ifndef __HAMMER_HKMER_DISTANCE_HPP__
#define __HAMMER_HKMER_DISTANCE_HPP__

#include "hkmer.hpp"
#include "err_helper_table.hpp"

namespace hammer {

enum IonEventType {
  kIonEventMismatch,
  kIonEventBaseInsertion,
  kIonEventRunInsertion,
  kIonEventBaseDeletion,
  kIonEventRunDeletion
};

template <typename It1, typename It2>
struct IonPairAlignEvent {
  IonEventType type;
  It1 x_iter;
  It2 y_iter;
  size_t length;
};

template<typename It1, typename It2>
class IonPairAligner {

  It1 x_it_;
  It1 x_end_;
  It2 y_it_;
  It2 y_end_;

  bool empty_;
  hammer::HomopolymerRun cx_, cy_;

  IonPairAlignEvent<It1, It2> front_;

  bool checkForZeroLengthRuns() {
    if (x_it_ == x_end_ || y_it_ == y_end_)
      return false;

    if (cx_.len > 0 && cy_.len > 0)
      return true;

    while (cx_.len == 0) {
      if (++x_it_ == x_end_)
        return false;
      cx_ = *x_it_;
    }

    while (cy_.len == 0) {
      if (++y_it_ == y_end_)
        return false;
      cy_ = *y_it_;
    }

    return true;
  }

  bool fetchNextX() {
    ++x_it_;
    if (x_it_ == x_end_)
      return false;
    cx_ = *x_it_;
    return true;
  }

  bool fetchNextY() {
    ++y_it_;
    if (y_it_ == y_end_)
      return false;
    cy_ = *y_it_;
    return true;
  }

  void updateEventOffset() {
    front_.x_iter = x_it_;
    front_.y_iter = y_it_;
  }

  void yieldBaseInsertion() {
    front_.type = kIonEventBaseInsertion;
    front_.length = cy_.len - cx_.len;
    updateEventOffset();
  }

  void yieldBaseDeletion() {
    front_.type = kIonEventBaseDeletion;
    front_.length = cx_.len - cy_.len;
    updateEventOffset();
  }

  void yieldRunEvent(IonEventType type) {
    front_.type = type;
    front_.length = 1;
    updateEventOffset();
  }

  void yieldMismatch() { yieldRunEvent(kIonEventMismatch); }
  void yieldRunInsertion() { yieldRunEvent(kIonEventRunInsertion); }
  void yieldRunDeletion() { yieldRunEvent(kIonEventRunDeletion); }

 public:
  IonPairAligner(const It1 &x_begin, const It1 &x_end,
                 const It2 &y_begin, const It2 &y_end)
    : x_it_(x_begin), x_end_(x_end), y_it_(y_begin), y_end_(y_end),
        empty_(false), cx_(*x_it_), cy_(*y_it_)
  {
    popFront();
  }

  bool empty() const { return empty_; }

  IonPairAlignEvent<It1, It2> front() const { return front_; }

  void popFront() {
    VERIFY(x_it_ <= x_end_);
    VERIFY(y_it_ <= y_end_);

    while (true) {
      if (x_it_ == x_end_ || y_it_ == y_end_ || !checkForZeroLengthRuns()) {
        if (x_it_ == x_end_ && y_it_ != y_end_ && cy_.len < y_it_->len) {
          cx_ = *(x_it_ - 1);
          cx_.len = 0;
          yieldBaseInsertion();
          fetchNextY();
          return;
        }

        if (y_it_ == y_end_ && x_it_ != x_end_ && cx_.len < x_it_->len) {
          cy_ = *(y_it_ - 1);
          cy_.len = 0;
          yieldBaseDeletion();
          fetchNextX();
          return;
        }

        empty_ = true;
        return;
      }

      bool end = false;
      while (cx_.raw == cy_.raw)
        if (!fetchNextX() || !fetchNextY()) {
          end = true;
          break;
        }

      if (!end)
        break;
    }

    if (cx_.nucl == cy_.nucl) {

      if (cx_.len >= 4 && cy_.len >= 4) {
        if (cx_.len < cy_.len)
          yieldBaseInsertion();
        else
          yieldBaseDeletion();
        fetchNextX();
        fetchNextY();
        return;
      } else {
        --cx_.len;
        --cy_.len;
        popFront();
      }

    } else {

      using namespace hammer::errHelper;
      auto hint = getHint(x_it_, x_end_, y_it_, y_end_, cx_.len, cy_.len);

      switch (hint) {
        case kMismatch:
          --cx_.len, --cy_.len;
          yieldMismatch();
          return;
        case kInsertion:
          --cy_.len;
          yieldRunInsertion();
          return;
        case kDeletion:
          --cx_.len;
          yieldRunDeletion();
          return;
      }
    }
  }
};

template <int kMismatchCost=1,
          int kBaseInsertionCost=1, int kRunInsertionCost=1,
          int kBaseDeletionCost=1, int kRunDeletionCost=1,
          typename It1, typename It2>
inline unsigned distanceHKMer(const It1 &x_begin, const It1 &x_end,
                              const It2 &y_begin, const It2 &y_end,
                              unsigned tau = -1) {
  unsigned dist = 0;

  IonPairAligner<It1, It2> aligner(x_begin, x_end, y_begin, y_end);

  while (!aligner.empty()) {
    auto event = aligner.front();
    switch (event.type) {
      case kIonEventMismatch:
        dist += kMismatchCost * event.length; break;
      case kIonEventBaseInsertion:
        dist += kBaseInsertionCost * event.length; break;
      case kIonEventBaseDeletion:
        dist += kBaseDeletionCost * event.length; break;
      case kIonEventRunInsertion:
        dist += kRunInsertionCost * event.length; break;
      case kIonEventRunDeletion:
        dist += kRunDeletionCost * event.length; break;
      default: break;
    }
    if (dist > tau)
      return dist;
    aligner.popFront();
  }

  return dist;
}

#include <cassert>
#include <iostream>
namespace unittest {

  namespace detail {

    typedef hammer::HKMer::StorageType::const_iterator It;

    inline unsigned distanceHKMer(It beg1, It end1, It beg2, It end2) {
      unsigned dist = 0;

      IonPairAligner<It, It> aligner(beg1, end1, beg2, end2);

      while (!aligner.empty()) {
        auto event = aligner.front();
        switch (event.type) {
        case kIonEventMismatch:
          std::cerr << event.length << 'X';
          dist += event.length; break;
        case kIonEventBaseInsertion:
          std::cerr << event.length << 'I';
          dist += event.length; break;
        case kIonEventBaseDeletion:
          std::cerr << event.length << 'D';
          dist += event.length; break;
        case kIonEventRunInsertion:
          std::cerr << event.length << 'I';
          dist += event.length; break;
        case kIonEventRunDeletion:
          std::cerr << event.length << 'D';
          dist += event.length; break;
        default: break;
        }
        aligner.popFront();
      }

      std::cerr << std::endl;
      return dist;
    }

    inline unsigned distance(const std::string &s, const std::string &t) {
      using namespace hammer;
      HKMer k1, k2;
      for (size_t i = 0; i < s.size(); ++i)
        k1 <<= s[i];
      for (size_t i = 0; i < t.size(); ++i)
        k2 <<= t[i];

      return distanceHKMer(k1.begin(), k1.end(), k2.begin(), k2.end());
    }
  }

  inline void hkmer_distance() {
    using namespace detail;
    assert(distance("GATAGCGATTTGTTCGGTTTAGGGGGGG",
                    "GATAGCGATTTGTTCGTTTAG") >= 7);

    assert(distance("ATAGCGATTTGTTCGGTTTAGGGGGGGT",
                    "ATAGCGATTTGTTCGTTTAGA") >= 7);

    assert(distance("GATTTGTTCGGTTTAGGGGGGGTAGGGGGATTA",
                    "GATTTGTTCGTTTAGGGGGGGTAGGGGGATTA") == 1);

    assert(distance("TTAAGGCTTACAAAGACTGCGTTT",
                    "TTAAGGCTTACAAAGACTGCGTTTT") == 1);

    assert(distance("AAGGCTTACAAAGACTGCGTTTAA",
                    "AAGGCTTACAAAGACTGCGTA") >= 2);

    assert(distance("ACAAAGACTGCGTTTAAGAGC",
                    "ACAAAGACTGCGTTTTAAGAGC") == 1);

    assert(distance("CTAGGAATGAAAAAGAGAACAAGAA",
                    "CTAGGAATGAAAAAGAGAAAAAAGAATG") == 2);

    assert(distance("ACACACAGGGTTTTTGAACTGGATT",
                    "ACACACAGGGTTTTGAACTGGATT") == 1);

    assert(distance("ACATAAGCCTTTGTACTTAGC",
                    "ACATAAGCCTTTGACTTAGCA") == 1);
  }
}


};
#endif // __HAMMER_HKMER_DISTANCE_HPP__
