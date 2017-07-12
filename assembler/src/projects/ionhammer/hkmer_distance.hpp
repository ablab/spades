//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_HKMER_DISTANCE_HPP__
#define __HAMMER_HKMER_DISTANCE_HPP__

#include "err_helper_table.hpp"
#include "hkmer.hpp"

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
  unsigned length;
};

template <typename It1, typename It2>
class IonPairAligner {
  It1 x_it_;
  It1 x_end_;
  It2 y_it_;
  It2 y_end_;

  bool empty_;
  hammer::HomopolymerRun cx_, cy_;
  int end_diff_;
  bool at_the_start_;  // turned off once we find a pair of runs with same
                       // nucleotide

  IonPairAlignEvent<It1, It2> front_;

  // true iff alignment process is not yet finished
  bool checkForZeroLengthRuns() {
    if (x_it_ == x_end_ || y_it_ == y_end_) return false;

    if (cx_.len > 0 && cy_.len > 0) return true;

    bool result = true;
    while (cx_.len == 0) {
      if (++x_it_ == x_end_) {
        result = false;
        break;
      }
      cx_ = *x_it_;
    }

    while (cy_.len == 0) {
      if (++y_it_ == y_end_) {
        result = false;
        break;
      }
      cy_ = *y_it_;
    }

    return result;
  }

  bool fetchNextX() {
    ++x_it_;
    if (x_it_ == x_end_) return false;
    cx_ = *x_it_;
    return true;
  }

  bool fetchNextY() {
    ++y_it_;
    if (y_it_ == y_end_) return false;
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

  void finishAlignmentProcess() {
    empty_ = true;
    if (x_it_ != x_end_) {
      end_diff_ += int(x_end_ - x_it_);
    }
    if (y_it_ != y_end_) {
      end_diff_ -= int(y_end_ - y_it_);
    }
  }

 public:
  IonPairAligner(const It1 &x_begin, const It1 &x_end, const It2 &y_begin,
                 const It2 &y_end)
      : x_it_(x_begin),
        x_end_(x_end),
        y_it_(y_begin),
        y_end_(y_end),
        empty_(false),
        cx_(*x_it_),
        cy_(*y_it_),
        end_diff_(0),
        at_the_start_(true) {
    popFront();
  }

  bool empty() const { return empty_; }

  IonPairAlignEvent<It1, It2> front() const { return front_; }

  // see comment to distanceHKMer function
  int endDiff() const { return end_diff_; }

  void popFront() {
    VERIFY(x_it_ <= x_end_);
    VERIFY(y_it_ <= y_end_);

    while (true) {
      if (x_it_ == x_end_ || y_it_ == y_end_ || !checkForZeroLengthRuns()) {
        if (x_it_ == x_end_ && y_it_ != y_end_ && cy_.len < y_it_->len) {
          cx_ = *(x_it_ - 1);
          cx_.len = 0;
          // end_diff_ -= 1;
          yieldBaseInsertion();
          fetchNextY();
          return;
        }

        if (y_it_ == y_end_ && x_it_ != x_end_ && cx_.len < x_it_->len) {
          cy_ = *(y_it_ - 1);
          cy_.len = 0;
          // end_diff_ += 1;
          yieldBaseDeletion();
          fetchNextX();
          return;
        }

        finishAlignmentProcess();
        return;
      }

      bool end = false;
      while (cx_.raw == cy_.raw) {
        // don't short-circuit for correct end_diff_ calculation!
        end = !fetchNextX();
        end |= !fetchNextY();
        if (end) break;
        at_the_start_ = false;
      }

      if (!end) break;
    }

    if (!checkForZeroLengthRuns()) {
      finishAlignmentProcess();
      return;
    }

    VERIFY(cx_.len > 0);
    VERIFY(cy_.len > 0);
    VERIFY(x_it_ < x_end_);
    VERIFY(y_it_ < y_end_);

    if (cx_.nucl == cy_.nucl) {
      if (cx_.len >= 4 && cy_.len >= 4) {
        if (cx_.len < cy_.len)
          yieldBaseInsertion();
        else
          yieldBaseDeletion();
        fetchNextX();
        fetchNextY();
        at_the_start_ = false;
        return;
      } else {
        --cx_.len;
        --cy_.len;
        at_the_start_ = false;
        popFront();
      }

    } else if (at_the_start_) {
      // alignment can't start with a deletion or insertion
      // unless it is in a homopolymer run
      --cx_.len, --cy_.len;
      yieldMismatch();
      return;
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

// returns distance between two homopolymer sequences;
// optionally, fills *end_diff:
//  [ --------- X ----------- ]
// [---------- Y -------]######
//                       \____/
//                       end_diff
template <int kMismatchCost = 1, int kBaseInsertionCost = 1,
          int kRunInsertionCost = 1, int kBaseDeletionCost = 1,
          int kRunDeletionCost = 1, typename It1, typename It2>
inline unsigned distanceHKMer(const It1 &x_begin, const It1 &x_end,
                              const It2 &y_begin, const It2 &y_end,
                              unsigned tau = -1, int *end_diff = NULL) {
  unsigned dist = 0;

  IonPairAligner<It1, It2> aligner(x_begin, x_end, y_begin, y_end);

  while (!aligner.empty()) {
    auto event = aligner.front();
    switch (event.type) {
      case kIonEventMismatch:
        dist += kMismatchCost * event.length;
        break;
      case kIonEventBaseInsertion:
        dist += kBaseInsertionCost * event.length;
        break;
      case kIonEventBaseDeletion:
        dist += kBaseDeletionCost * event.length;
        break;
      case kIonEventRunInsertion:
        dist += kRunInsertionCost * event.length;
        break;
      case kIonEventRunDeletion:
        dist += kRunDeletionCost * event.length;
        break;
      default:
        break;
    }
    if (dist > tau && end_diff == NULL) break;
    aligner.popFront();
  }

  if (end_diff != NULL) *end_diff = aligner.endDiff();

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
        dist += event.length;
        break;
      case kIonEventBaseInsertion:
        std::cerr << event.length << 'I';
        dist += event.length;
        break;
      case kIonEventBaseDeletion:
        std::cerr << event.length << 'D';
        dist += event.length;
        break;
      case kIonEventRunInsertion:
        std::cerr << event.length << 'I';
        dist += event.length;
        break;
      case kIonEventRunDeletion:
        std::cerr << event.length << 'D';
        dist += event.length;
        break;
      default:
        break;
    }
    aligner.popFront();
  }

  std::cerr << " (end. diff. = " << aligner.endDiff() << ")" << std::endl;
  return dist;
}

inline unsigned distance(const std::string &s, const std::string &t) {
  using namespace hammer;
  HKMer k1, k2;
  for (size_t i = 0; i < s.size(); ++i) k1 <<= s[i];
  for (size_t i = 0; i < t.size(); ++i) k2 <<= t[i];

  return distanceHKMer(k1.begin(), k1.end(), k2.begin(), k2.end());
}
}  // namespace detail

inline void hkmer_distance() {
  using namespace detail;

  assert(distance("ACGTACGTACGTACGT", "CGTACGTACGTACGTA") > 1);

  assert(distance("AACGTACGTACGTACGT", "CGTACGTACGTACGTA") > 1);

  assert(distance("GATAGCGATTTGTTCGGTTTAGGGGGGG", "GATAGCGATTTGTTCGTTTAG") >=
         7);

  assert(distance("ATAGCGATTTGTTCGGTTTAGGGGGGGT", "ATAGCGATTTGTTCGTTTAGA") >=
         7);

  assert(distance("GATTTGTTCGGTTTAGGGGGGGTAGGGGGATTA",
                  "GATTTGTTCGTTTAGGGGGGGTAGGGGGATTA") == 1);

  assert(distance("TTAAGGCTTACAAAGACTGCGTTT", "TTAAGGCTTACAAAGACTGCGTTTT") ==
         1);

  assert(distance("AAGGCTTACAAAGACTGCGTTTAA", "AAGGCTTACAAAGACTGCGTA") >= 2);

  assert(distance("ACAAAGACTGCGTTTAAGAGC", "ACAAAGACTGCGTTTTAAGAGC") == 1);

  assert(distance("CTAGGAATGAAAAAGAGAACAAGAA",
                  "CTAGGAATGAAAAAGAGAAAAAAGAATG") == 2);

  assert(distance("ACACACAGGGTTTTTGAACTGGATT", "ACACACAGGGTTTTGAACTGGATT") ==
         1);

  assert(distance("ACATAAGCCTTTGTACTTAGC", "ACATAAGCCTTTGACTTAGCA") == 1);
}
}  // namespace unittest

};      // namespace hammer
#endif  // __HAMMER_HKMER_DISTANCE_HPP__
