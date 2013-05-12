#ifndef JOBWRAPERS_H_
#define JOBWRAPERS_H_

#include "running_modes.hpp"
#include "QcException.hpp"
#include "aho_corasick.hpp"
#include "output.hpp"
#include "config_struct_cclean.hpp"
#include "io/read_processor.hpp"
#include "ssw/ssw_cpp.h"

class AlignmentJobWrapper {
public:
  AlignmentJobWrapper(const Database * data, std::ostream& output, std::ostream& bed)
      :data(data), output(output), bed(bed),
       mismatch_threshold(cfg::get().mismatch_threshold), aligned_part_fraction(cfg::get().aligned_part_fraction),
       aligned_(0) {};

  bool operator()(const Read &r);

  size_t aligned() const { return aligned_; }

private:
  const Database * data;
  std::ostream& output;
  std::ostream& bed;
  unsigned mismatch_threshold;
  double aligned_part_fraction;
  size_t aligned_;
};

class ExactMatchJobWrapper {
public:
  ExactMatchJobWrapper(const Database * data, std::ostream& output, std::ostream& bed, const AhoCorasick &a)
      :data(data), output(output), bed(bed), ahoCorasick(a), aligned_(0) {};

  bool operator()(const Read &r);

  size_t aligned() const { return aligned_; }

private:
  const Database * data;
  std::ostream& output;
  std::ostream& bed;
  AhoCorasick ahoCorasick;
  size_t aligned_;
};

class ExactAndAlignJobWrapper {
public:
  ExactAndAlignJobWrapper(const Database * data, std::ostream& output, std::ostream& bed, const AhoCorasick &a, const AhoCorasick &b)
      :data(data), output(output), bed(bed), dbAhoCorasick(a), kmersAhoCorasick(b),
       mismatch_threshold(cfg::get().mismatch_threshold), aligned_part_fraction(cfg::get().aligned_part_fraction),
       aligned_(0) {};

  bool operator()(const Read &r);

  size_t aligned() const { return aligned_; }

private:
  const Database * data;
  std::ostream& output;
  std::ostream& bed;
  AhoCorasick dbAhoCorasick;
  AhoCorasick kmersAhoCorasick;
  unsigned mismatch_threshold;
  double aligned_part_fraction;
  size_t aligned_;
};

#endif /* JOBWRAPERS_H_ */
