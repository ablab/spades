// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef READ_FILL_RECORD
#define READ_FILL_RECORD

#include "system/System.h"

/// The file reads.filled_fillrecords.k48 in the run_dir is a 
/// vec<ReadFillRecord> (one entry per read) created by CloseAllReadGaps.
/// Hand either the filename or the vector to a ReadFillDatabase
/// to get at this information in an encapsulated way.

class ReadFillRecord {
public:
  ReadFillRecord() {}
  ReadFillRecord(int first, int num, bool expl, bool merg)
    : first_filling(first), 
      num_fillings(num), 
      exploded(expl),
      merger_limit(merg)
  {}

  int FirstFilling() const { return first_filling; }
  int LastFilling() const { return first_filling + num_fillings - 1; }
  int NumFillings() const { return num_fillings; }
  bool Exploded() const { return exploded; }
  bool HitMergerLimit() const { return merger_limit; }

  void SetFirstFilling(int first) { first_filling = first; }
  void SetNumFillings(int num) { num_fillings = num; }
  void SetExploded(bool expl) { exploded = expl; }
  void SetHitMergerLimit(bool merg) { merger_limit = merg; }

  void AddResults( int num, bool expl, bool merg )
  { num_fillings += num; exploded |= expl; merger_limit |= merg; }

  // So we can save and load a vec of these:
  friend ostream& operator<<(ostream& out, const ReadFillRecord& my) {
    out << my.first_filling << " " 
	<< my.num_fillings << " "
	<< my.exploded << " "
	<< my.merger_limit << "\n";
    return out;
  }
  friend istream& operator>>(istream& in, ReadFillRecord& my) {
    in >> my.first_filling 
       >> my.num_fillings 
       >> my.exploded
       >> my.merger_limit;
    in.ignore(1);
    return in;
  }


private:
  int first_filling;
  int num_fillings;
  Bool exploded;
  Bool merger_limit;
};


#endif
