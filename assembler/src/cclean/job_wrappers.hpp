#ifndef JOBWRAPERS_H_
#define JOBWRAPERS_H_

#include "running_modes.hpp"
#include "output.hpp"
#include "config_struct_cclean.hpp"
#include "io/read_processor.hpp"

namespace cclean {
class AdapterIndex;
}

class SimpleClean {
public:
  SimpleClean(std::ostream& output, std::ostream& bed, const std::string &db,
              const cclean::AdapterIndex &index)
      : output(output), bed(bed),
        mismatch_threshold(cfg::get().mismatch_threshold),
        aligned_part_fraction(cfg::get().aligned_part_fraction),
        aligned_(0), index_(index), db_(db) {}

  bool operator()(const Read &r);
  inline size_t aligned() const { return aligned_; }

private:
  std::ostream& output;
  std::ostream& bed;
  const std::string& db_;
  unsigned mismatch_threshold;
  double aligned_part_fraction;
  size_t aligned_;
  const cclean::AdapterIndex &index_;
};

#endif /* JOBWRAPERS_H_ */
