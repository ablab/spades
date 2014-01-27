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
  SimpleClean(std::ostream& aligned_output, std::ostream& output,
              std::ostream& bed, const std::string &db,
              const cclean::AdapterIndex &index)
    : aligned_output_(aligned_output), bed_(bed), output_stream_(output),
        mismatch_threshold_(cfg::get().mismatch_threshold),
        aligned_part_fraction_(cfg::get().aligned_part_fraction),
        aligned_(0), index_(index), db_(db) {}

  bool operator()(const Read &r);

  inline size_t aligned() const { return aligned_; }

private:
  std::ostream& output_stream_;
  std::ostream& aligned_output_;
  std::ostream& bed_;
  const std::string& db_;
  unsigned mismatch_threshold_;
  double aligned_part_fraction_;
  size_t aligned_;
  const cclean::AdapterIndex &index_;
};

#endif /* JOBWRAPERS_H_ */
