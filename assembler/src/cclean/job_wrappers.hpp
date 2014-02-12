#ifndef JOBWRAPERS_H_
#define JOBWRAPERS_H_

#include "running_modes.hpp"
#include "output.hpp"
#include "config_struct_cclean.hpp"
#include "io/read_processor.hpp"
#include "additional.cpp"

namespace cclean {
class AdapterIndex;
}

class SimpleClean {
public:
  SimpleClean(std::ostream& aligned_output, std::ostream& output,
              std::ostream& bed, const std::string &db,
              const cclean::AdapterIndex &index,
              const additional::WorkModeType &mode,
              const std::unordered_map<std::string, std::string> &options)
    : aligned_output_(aligned_output), bed_(bed), output_stream_(output),
        mismatch_threshold_(cfg::get().mismatch_threshold),
        aligned_part_fraction_(cfg::get().aligned_part_fraction),
        aligned_(0), index_(index), db_(db), options_(options), mode_(mode)
        { read_mlen_ = atoi(options_["mlen"].c_str()); }

  bool operator()(const Read &r);
  inline void SingleEndClean(const Read &r);
  inline void SingleEndWQualityClean(const Read &r);

  inline size_t aligned() const { return aligned_; }

private:
  std::unordered_map<std::string, std::string> options_;
  const std::string& db_;
  std::ostream& output_stream_;
  std::ostream& aligned_output_;
  std::ostream& bed_;
  const additional::WorkModeType &mode_;
  double aligned_part_fraction_;
  size_t aligned_;
  const cclean::AdapterIndex &index_;
  unsigned mismatch_threshold_;
  unsigned read_mlen_;
};

#endif /* JOBWRAPERS_H_ */
