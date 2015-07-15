//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_EXPANDER_HPP__
#define __HAMMER_EXPANDER_HPP__

class KMerData;
namespace io {
class SingleRead;
}

#include <cstring>

class Expander {
  KMerData &data_;
  size_t changed_;
  
 public:
  Expander(KMerData &data) 
      : data_(data), changed_(0) {}

  size_t changed() const { return changed_; }

  bool operator()(const io::SingleRead &r);
};

#endif
