//
// Created by Vasiliy Ershov on 16/07/16.
//

#ifndef PROJECT_QUALITY_THRESHOLDS_ESTIMATOR_H
#define PROJECT_QUALITY_THRESHOLDS_ESTIMATOR_H

#include <algorithm>
#include <cmath>
#include <vector>

class SimpleTwoClassClustering {
 private:
  std::vector<double> samples_;

  static inline double Score(double sum, double weight) {
    return weight > 3 ? -sum * sum * (1 + 2 * log(weight + 1)) / weight : 0;
  }

  struct BestSplit {
    double score_;
    double split_;
    double left_sum_;
    double left_weight_;
    double right_sum_;
    double right_weight_;
  };

 public:
  SimpleTwoClassClustering(size_t maxSize = 100) { samples_.reserve(maxSize); }

  void Add(double sample) { samples_.push_back(sample); }

  double FindBestSpit() {
    std::sort(samples_.begin(), samples_.end());
    auto bestSplit = SimpleThresholdEstimation(samples_.begin(), samples_.end());
    return bestSplit.split_;
  }

  double EstimateAlpha() {
    double minSample = 0;
    std::sort(samples_.begin(), samples_.end());

    for (auto sampl : samples_) {
      minSample = std::min(sampl, minSample);
    }
    auto bestSplit = SimpleThresholdEstimation(samples_.begin(), samples_.end());
    const double p = 0.5;
    double alpha = log(1.0 - p) / bestSplit.split_;
    return alpha;
  }

  // it's simple decision tree with quality target
  static BestSplit SimpleThresholdEstimation(std::vector<double>::const_iterator from,
                                              std::vector<double>::const_iterator to) {
    const double total_sum = [&]() -> double {
      double sum = 0;
      for (auto sorted_samples_iterator = from; sorted_samples_iterator != to; ++sorted_samples_iterator) {
        const auto sample = *sorted_samples_iterator;
        sum += sample;
      }
      return sum;
    }();
    const double total_weight = (double)(to - from);

    double best_score = 0;
    double best_left_sum = 0;
    double best_left_weight = 0;
    double best_split = 0;

    double sum = 0;
    double weight = 0;

    for (auto sorted_samples_iterator = from; sorted_samples_iterator != to; ++sorted_samples_iterator) {
      const auto sample = *sorted_samples_iterator;
      sum += sample;
      ++weight;

      const double right_leaf_sum = total_sum - sum;
      const double right_leaf_weight = total_weight - weight;
      const double split_score =
          Score(sum, weight) + Score(right_leaf_sum, right_leaf_weight);

      if (split_score <= best_score) {
        best_score = split_score;
        best_left_weight = weight;
        best_left_sum = sum;
        best_split = sample;
      }
    }

    return {best_score,
            best_split,
            best_left_sum,
            best_left_weight,
            total_sum - best_left_sum,
            total_weight - best_left_weight};
  }
};

#endif  // PROJECT_QUALITY_THRESHOLDS_ESTIMATOR_H
