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

  bool fetchNextPairOfRuns() {
    ++x_it_;
    ++y_it_;
    if (x_it_ == x_end_ || y_it_ == y_end_)
      return false;
    cx_ = *x_it_;
    cy_ = *y_it_;
    return checkForZeroLengthRuns();
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
    if (checkForZeroLengthRuns())
      popFront();
    else
      empty_ = true;
  }

  bool empty() const { return empty_; }

  IonPairAlignEvent<It1, It2> front() const { return front_; }

  void popFront() {
    while (true) {
      VERIFY(cx_.len > 0 && cy_.len > 0);
      if (cx_.raw == cy_.raw) {
        if (!fetchNextPairOfRuns())
          break;
        else
          continue;
      }

      if (cx_.nucl == cy_.nucl) {

        if (cx_.len >= 4 && cy_.len >= 4) {
          if (cx_.len < cy_.len)
            yieldBaseInsertion();
          else 
            yieldBaseDeletion();
          if (!fetchNextPairOfRuns())
            break;
          return;
        } else {
          --cx_.len;
          --cy_.len;
          if (!checkForZeroLengthRuns())
            break;
        }

      } else {

        using namespace hammer::errHelper;
        auto hint = getHint(x_it_, x_end_, y_it_, y_end_, cx_.len, cy_.len);
        
        switch (hint) {
          case kMismatch:
            --cx_.len, --cy_.len;
            yieldMismatch();
            break;
          case kInsertion:
            --cy_.len;
            yieldRunInsertion();
            break;
          case kDeletion:
            --cx_.len;
            yieldRunDeletion();
            break;
        }

        if (!checkForZeroLengthRuns())
          break;
        return;
      }
    }

    empty_ = true;
  }
};

template <int kMismatchCost=1, 
          int kBaseInsertionCost=1, int kRunInsertionCost=1,
          int kBaseDeletionCost=1, int kRunDeletionCost=1,
          typename It1, typename It2>
inline size_t distanceHKMer(const It1 &x_begin, const It1 &x_end,
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

};
#endif // __HAMMER_HKMER_DISTANCE_HPP__
