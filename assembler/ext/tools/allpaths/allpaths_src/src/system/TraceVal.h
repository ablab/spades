///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   File: TraceVal.h

   Debugging utils for tracing the origin of bad values in a program.
   It lets you store for each integer, "I'm the result of the 12345th integer
   assignment!", so for a given incorrect value you can re-run the program
   and stop it at the point where the incorrect value was assigned.
   
   See class TraceVal.
*/

#ifndef __INCLUDE_TraceVal_h
#define __INCLUDE_TraceVal_h

#include <iostream>
#include "system/Types.h"
#include "system/Assert.h"

// Macros: Tracing macros
//   TRACEVAL_ON - undefine if you're not doing tracing
//   TRACEVAL_STOP_TRACING_COPIES -- call before an operation that shuffles
//       a vector of TraceVals.
//   TRACEVAL_START_TRACING_COPIES -- call after an operation that shuffles
//       a vector of TraceVals.

//#define TRACEVAL_ON


#ifdef TRACEVAL_ON
#define TRACEVAL_STOP_TRACING_COPIES TraceValCommon::StopTracingCopies()
#define TRACEVAL_START_TRACING_COPIES TraceValCommon::StartTracingCopies()
#else
#define TRACEVAL_STOP_TRACING_COPIES 
#define TRACEVAL_START_TRACING_COPIES
#endif

/**
   Class: TraceValCommon

   Defines things common to all instantiations of TraceVal.
*/
class TraceValCommon {
 public:
  typedef longlong timestamp_t;

  static void SetStopVal( timestamp_t stopVal ) {
    stopAtStamp_ = stopVal;
    cout << "\n*** TraceVal: will stop at time " << stopAtStamp_ << endl << endl;
  }
  
  static void StopTracingCopies() { tracingCopiesStopped_ = True; }
  static void StartTracingCopies() { tracingCopiesStopped_ = False; }

 private:
  template <class T> friend class TraceVal;
  
  static timestamp_t nextTimeStamp_;
  static timestamp_t stopAtStamp_;
  static Bool tracingCopiesStopped_;
};

/**
   Class: TraceVal

   Debugging class to help us trace the origin of a value that appears at some
   point in the program, by letting you stop the program at the precise moment
   the given value is created. Keeps a global counter (clock) of all values created,
   and saves with each value the time of its creation.  So, a TraceVal<int> holds
   its integer value AND the clock value ("I was the 193234th integer created!").
   Then, you can re-run the program with the special argument TV=193234, and it will stop and print
   a stacktrace at the moment the 193234th integer is created.  If the traced value
   was being constructed from another traced value, the program will also print
   the timestamp where _that_ value was created, so you can now re-run the program
   and stop it when _that_ value was created.   In this way you can trace the origin
   of suspect values.
*/
template <class T>
class TraceVal {
 public:
  
  TraceVal(): val_(0) { makeStamp(); }
  TraceVal(const TraceVal<T>& tv): val_(tv.val_) {
    makeStamp(tv);
  }
  TraceVal(const T& val): val_(val) {  makeStamp(); }

  operator T() const { return val_; }

  TraceVal<T>& operator=( const TraceVal<T>& tv ) { val_ = tv.val_; makeStamp(tv); return *this; }
  TraceVal<T>& operator=( const T& val ) { val_ = val; makeStamp(); return *this; }

  TraceVal<T>& operator++() { val_++; return *this; }
  TraceVal<T> operator++(int) { TraceVal<T> oldVal = *this; this->operator++(); return oldVal; }
  TraceVal<T>& operator--() { val_--; return *this; }
  TraceVal<T> operator--(int) { TraceVal<T> oldVal = *this; this->operator--(); return oldVal; }
  TraceVal<T>& operator+=( const TraceVal<T>& tv ) {
    val_ += tv.val_;
    makeStamp(tv);
    return *this;
  }
  TraceVal<T>& operator-=( const TraceVal<T>& tv ) {
    val_ -= tv.val_;
    makeStamp(tv);
    return *this;
  }

  template <class T2> friend std::istream& operator>> ( std::istream& in, TraceVal<T2>& tv );

  TraceValCommon::timestamp_t GetTimeStamp() const { return timeStamp_; }


protected:
  T val_;
  TraceValCommon::timestamp_t timeStamp_;


  void makeStamp() {
    if ( ( timeStamp_ = TraceValCommon::nextTimeStamp_++ ) == TraceValCommon::stopAtStamp_ ) {
      cout << "TraceVal: STOPPING BY REQUEST AT " << timeStamp_ << "th val." << endl;
      TracebackThisProcess();
    }
  }

  void makeStamp( const TraceVal<T>& src ) {
    if ( TraceValCommon::tracingCopiesStopped_ ) {
      timeStamp_ = src.timeStamp_;
     } else {
      if ( ( timeStamp_ = TraceValCommon::nextTimeStamp_++ ) == TraceValCommon::stopAtStamp_ ) {
	cout << "TraceVal: STOPPING BY REQUEST AT " << timeStamp_ << "th val." << endl;
	cout << "This val is being created from val created at time " << src.GetTimeStamp() << endl;
	TracebackThisProcess();
      }
    }
  }
  
  friend class TraceValInitializer;
  
};  // class TraceVal

template <class T>
inline std::istream& operator>> ( std::istream& in, TraceVal<T>& tv ) {
  in >> tv.val_;
  return in;
}

typedef TraceVal<int> TraceInt;

#endif
// #ifndef __INCLUDE_TraceVal_h
