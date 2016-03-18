//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_IT_FLOW_SPACE_READ_HPP__
#define __HAMMER_IT_FLOW_SPACE_READ_HPP__

#include "io/reads/single_read.hpp"
#include "HSeq.hpp"

#include <deque>
#include <cstddef>
#include <string>

namespace hammer {

/// Read interpreted as series of homopolymer runs
class FlowSpaceRead {
  std::string name_;
  std::deque<HomopolymerRun> runs_;
 public:
  FlowSpaceRead(const io::SingleRead& read) : name_(read.name()) {
    const auto& seq = read.GetSequenceString();
    hammer::iontorrent::toHomopolymerRuns(seq, runs_);
  }

  template <typename It>
  FlowSpaceRead(It runs_beg, It runs_end) :
    runs_(runs_beg, runs_end) {}

  size_t size() const {
    return runs_.size();
  }

  const std::string& name() const {
    return name_;
  }

  HomopolymerRun operator[](size_t index) const {
    return runs_[index];
  }

  HomopolymerRun& operator[](size_t index) {
    return runs_[index];
  }

  void TrimLeft(size_t n_runs) {
    if (n_runs >= runs_.size())
      runs_.clear();
    else
      runs_.erase(runs_.begin(), runs_.begin() + n_runs);
  }

  void TrimRight(size_t n_runs) {
    if (n_runs >= runs_.size())
      runs_.clear();
    else
      runs_.erase(runs_.end() - n_runs, runs_.end());
  }

  std::string GetSequenceString() const {
    std::string seq;
    for (size_t i = 0; i < runs_.size(); ++i)
      seq += runs_[i].str();
    return seq;
  }

  const std::deque<hammer::HomopolymerRun>& data() const {
    return runs_;
  }
};

} // namespace hammer
#endif
