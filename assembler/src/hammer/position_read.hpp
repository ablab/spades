//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef POSITION_READ_HPP_
#define POSITION_READ_HPP_

#include <vector>
#include <string>
#include <map>
#include <fstream>

#include "kmer_stat.hpp"
#include "globals.hpp"
#include "kmerno.hpp"

class PositionRead {
  hint_t start_  : 48;
  uint16_t size_ : 16;
  hint_t readno_ : 48;
  bool done_     : 1;
  unsigned res   : 15;

  public:
  PositionRead(hint_t start, unsigned size, hint_t readno, bool bad = false)
      : start_(start), size_(size), readno_(readno), done_(bad)  {
    VERIFY(size_ < 65536);
  }
  hint_t start() const { return start_; }
  uint32_t size() const { return size_; }
  char at(uint32_t pos) const;
  char operator [] (uint32_t pos) const;
  bool isDone() const { return done_; }
  void set_done(bool val = true) { done_ = val; }
  bool valid() const { return size_ >= K; }

  std::pair<int, hint_t> nextKMerNo(int begin) const;
};

#endif
