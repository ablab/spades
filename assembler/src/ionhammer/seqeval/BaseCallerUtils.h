/* Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     BaseCallerUtils.h
//! @ingroup  BaseCaller
//! @brief    Small helper classes and functions for BaseCaller

#ifndef BASECALLERUTILS_H
#define BASECALLERUTILS_H

#include <string>
#include <algorithm>
#include <vector>

namespace ion {

class FlowOrder {
public:
  FlowOrder() : cycle_nucs_("TACG"), num_flows_(0) {}
  FlowOrder(const std::string& cycle_nucs, int num_flows) : cycle_nucs_(cycle_nucs), num_flows_(num_flows) {
    BuildFullArrays();
  }

  void SetFlowOrder(const std::string& cycle_nucs, int num_flows) {
    cycle_nucs_ = cycle_nucs;
    num_flows_ = num_flows;
    BuildFullArrays();
  }
  void SetFlowOrder(const std::string& cycle_nucs) {
    cycle_nucs_ = cycle_nucs;
    BuildFullArrays();
  }
  void SetNumFlows(int num_flows) {
    num_flows_ = num_flows;
    BuildFullArrays();
  }

  char operator[](int flow) const { return full_nucs_.at(flow); }

  const std::string& str() const { return cycle_nucs_; }
  const char *c_str() const { return cycle_nucs_.c_str(); }
  bool is_ok() const { return num_flows_ > 0 and not cycle_nucs_.empty(); }

  int int_at(int flow) const { return full_ints_.at(flow); }
  char nuc_at(int flow) const { return full_nucs_.at(flow); }

  int num_flows() const { return num_flows_; }
  const char *full_nucs() const { return full_nucs_.c_str(); }

//  void BasesToFlows (const string& basespace, vector<int> &flowspace) const
//  void BasesToFlows (const string& basespace, vector<int> &flowspace, int num_flows) const;

  // Potential use case: In BaseCaller, BaseCallerLite:
  //void FlowsToBases, automatically resizes bases, returns number of bases
  //void FlowsToFlowIndex, also generates flow_index...
  //


  int BasesToFlows (const std::string& seq, int *ionogram, int ionogram_space) const
  {
    int flows = 0;
    unsigned int bases = 0;
    while (flows < ionogram_space and flows < num_flows_ and bases < seq.length()) {
      ionogram[flows] = 0;
      while (bases < seq.length() and full_nucs_[flows] == seq[bases]) {
        ionogram[flows]++;
        bases++;
      }
      flows++;
    }
    return flows;
  }

  int FlowsToBases(const std::vector<char>& ionogram, std::string &bases) const
  {
    int num_bases = 0;
    for(int flow = 0; flow < (int)ionogram.size() and flow < num_flows_; ++flow)
      num_bases += ionogram[flow];
    bases.clear();
    bases.reserve(num_bases+1);
    for(int flow = 0; flow < (int)ionogram.size() and flow < num_flows_; ++flow)
      for (int hp = 0; hp < ionogram[flow]; ++hp)
        bases.push_back(nuc_at(flow));
    return num_bases;
  }

  static int NucToInt (char nuc) {
      if (nuc=='a' or nuc=='A') return 0;
      if (nuc=='c' or nuc=='C') return 1;
      if (nuc=='g' or nuc=='G') return 2;
      if (nuc=='t' or nuc=='T') return 3;
      return -1;
  }
  static char IntToNuc (int idx) {
      if (idx == 0) return 'A';
      if (idx == 1) return 'C';
      if (idx == 2) return 'G';
      if (idx == 3) return 'T';
      return 'N';
  }

private:
  void BuildFullArrays() {
    if (cycle_nucs_.empty())
      return;
    full_nucs_.resize(num_flows_);
    full_ints_.resize(num_flows_);
    for (int flow = 0; flow < num_flows_; ++flow) {
      full_nucs_[flow] = cycle_nucs_[flow % cycle_nucs_.size()];
      full_ints_[flow] = NucToInt(full_nucs_[flow]);
    }
  }

  std::string       cycle_nucs_;
  std::string       full_nucs_;
  std::vector<int>  full_ints_;
  int               num_flows_;
};

}

/*
// Idea: class that manages splitting the chip into equal-size regions.
// Application: BaseCaller's wells reading, PhasingEstimator wells reading and cafie regions.
class ChipPartition {

};
*/

#define MAX_KEY_FLOWS     64

class KeySequence {
public:
  KeySequence() : bases_length_(0), flows_length_(0) { flows_[0]=0; }
  KeySequence(const ion::FlowOrder& flow_order, const std::string& key_string, const std::string& key_name) {
    Set(flow_order, key_string, key_name);
  }

/*  void Set(const std::string& flow_order, const std::string& key_string, const std::string& key_name) {
    name_ = key_name;
    bases_ = key_string;
    transform(bases_.begin(), bases_.end(), bases_.begin(), ::toupper);
    bases_length_ = key_string.length();
    flows_length_ = seqToFlow(bases_.c_str(), bases_length_, &flows_[0], MAX_KEY_FLOWS,
          (char *)flow_order.c_str(), flow_order.length());
  }*/
  void Set(const ion::FlowOrder& flow_order, const std::string& key_string, const std::string& key_name) {
    name_ = key_name;
    bases_ = key_string;
    transform(bases_.begin(), bases_.end(), bases_.begin(), ::toupper);
    bases_length_ = key_string.length();
    flows_length_ = flow_order.BasesToFlows(bases_, flows_, MAX_KEY_FLOWS);
  }

  const std::string &   name() const { return name_; }
  const std::string &   bases() const { return bases_; }
  size_t                bases_length() const { return bases_length_; }
  const int *           flows() const { return flows_; }
  int                   flows_length() const { return flows_length_; }
  int                   operator[](int i) const { return flows_[i]; }
  const char *          c_str() const { return bases_.c_str(); }

  // Feature request: keypass a basespace read
  // Feature request: keypass a flowspace read (int or float)

private:
  std::string           name_;
  std::string           bases_;
  size_t                bases_length_;
  int                   flows_length_;
  int                   flows_[MAX_KEY_FLOWS];
};




#endif // BASECALLERUTILS_H
