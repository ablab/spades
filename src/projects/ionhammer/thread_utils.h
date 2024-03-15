//
// Created by Vasiliy Ershov on 02/12/2016.
//

#ifndef PROJECT_THREAD_UTILS_H
#define PROJECT_THREAD_UTILS_H

namespace n_computation_utils {

template <class AdditiveStat>
class ParallelStatisticsCalcer {
 private:
  size_t num_threads_;

 public:
  ParallelStatisticsCalcer(size_t num_threads) : num_threads_(num_threads) {}

  template <class TFunction>
  AdditiveStat Calculate(size_t n, std::function<AdditiveStat()>&& factory,
                          TFunction&& func) const {
    std::vector<AdditiveStat> aggregated_stats;
    for (uint i = 0; i < num_threads_; ++i) {
      aggregated_stats.push_back(factory());
    }

#pragma omp parallel for num_threads(num_threads_)
    for (size_t i = 0; i < n; ++i) {
      const auto tid = omp_get_thread_num();
      func(aggregated_stats[tid], i);
    }

    for (size_t i = 1; i < aggregated_stats.size(); ++i) {
      aggregated_stats[0] += aggregated_stats[i];
    }
    return aggregated_stats[0];
  }
};

template <class TStat, class TAdditiveStat>
class TAdditiveStatisticsCalcer {
 private:
  const std::vector<TStat>& stats_;
  size_t num_threads_;

 public:
  TAdditiveStatisticsCalcer(const std::vector<TStat>& stats, size_t num_threads)
      : stats_(stats), num_threads_(num_threads) {}

  TAdditiveStat Calculate(std::function<TAdditiveStat()>&& factory) const {
    ParallelStatisticsCalcer<TAdditiveStat> parallel_calcer(num_threads_);
    return parallel_calcer.Calculate(
        stats_.size(), std::move(factory),
        [&](TAdditiveStat& stat, size_t i) { stat.Add(stats_[i]); });
  }
};
}  // namespace n_computation_utils
#endif  // PROJECT_THREAD_UTILS_H
