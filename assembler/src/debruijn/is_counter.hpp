/*
 * is_counter.hpp
 *
 *  Created on: May 25, 2014
 *      Author: andrey
 */

#ifndef IS_COUNTER_HPP_
#define IS_COUNTER_HPP_


#include "de/insert_size_refiner.hpp"
#include "sequence_mapper_notifier.hpp"

namespace debruijn_graph {

class InsertSizeCounter: public SequenceMapperListener {
    typedef std::map<int, size_t> HistType;

public:

    InsertSizeCounter(const conj_graph_pack& gp,
            size_t edge_length_threshold,
            bool ignore_negative = false)
        : gp_(gp), hist_(), tmp_hists_(),
          total_(), counted_(), negative_(),
          k_(gp.k_value), edge_length_threshold_(edge_length_threshold),
          ignore_negative_(ignore_negative) {
    }

    HistType hist() { return hist_; }
    size_t total() const { return total_.total_; }
    size_t mapped() const { return counted_.total_; }
    size_t negative() const { return negative_.total_; }


    virtual void StartProcessLibrary(size_t threads_count) {
        hist_.clear();
        tmp_hists_ = vector<HistType*>(threads_count);
        tmp_hists_[0] = &hist_;
        for (size_t i = 1; i < threads_count; ++i)
            tmp_hists_[i] = new HistType();

        total_ = count_data(threads_count);
        counted_ = count_data(threads_count);
        negative_ = count_data(threads_count);
    }

    virtual void StopProcessLibrary() {
        for (size_t i = 1; i < tmp_hists_.size(); ++i) {
            MergeBuffer(i);
            delete tmp_hists_[i];
        }
        total_.merge();
        counted_.merge();
        negative_.merge();
    }

    virtual void ProcessPairedRead(size_t thread_index,
                                   const MappingPath<EdgeId>& read1,
                                   const MappingPath<EdgeId>& read2,
                                   size_t dist) {

        ++total_.arr_[thread_index];

        if (read1.size() == 1 && read2.size() == 1 &&
                read2.simple_path().front() == read1.simple_path().front() &&
                gp_.g.length(read1.simple_path().front()) >= edge_length_threshold_) {

            auto mapping_edge_1 = read1.front().second;
            auto mapping_edge_2 = read2.front().second;

            int read1_start = (int) mapping_edge_1.mapped_range.start_pos - (int) mapping_edge_1.initial_range.start_pos ;
            TRACE("Read 1: " << (int) mapping_edge_1.mapped_range.start_pos << " - " << (int) mapping_edge_1.initial_range.start_pos << " = " << read1_start);
            int read2_start = (int) mapping_edge_2.mapped_range.start_pos - (int) mapping_edge_2.initial_range.start_pos;
            TRACE("Read 2: " << (int) mapping_edge_2.mapped_range.start_pos << " - " << (int) mapping_edge_2.initial_range.start_pos << " = " << read2_start);
            int is = read2_start - read1_start + (int) dist;
            TRACE("IS: " << read2_start << " - " <<  read1_start << " + " << (int) dist << "(" << dist << ") = " << is);

            if (is > 0 || !ignore_negative_) {
                (*tmp_hists_[thread_index])[is] += 1;
                ++counted_.arr_[thread_index];
            } else {
                ++negative_.arr_[thread_index];
            }

        }

    }

    virtual void ProcessSingleRead(size_t /*thread_index*/, const MappingPath<EdgeId>& /*read*/) {

    }

    virtual void MergeBuffer(size_t thread_index) {
        if (thread_index != 0) {
            for (auto it = tmp_hists_[thread_index]->begin(); it != tmp_hists_[thread_index]->end(); ++it) {
                (*tmp_hists_[0])[it->first] += it->second;
            }
            tmp_hists_[thread_index]->clear();
        }
    }

    void FindMean(double& mean, double& delta, std::map<size_t, size_t>& percentiles) const {
        double median = get_median(hist_);
        double mad = get_mad(hist_, median);
        double low = median - 5. * 1.4826 * mad;
        double high = median + 5. * 1.4826 * mad;

        DEBUG("Median IS: " << median);
        DEBUG("MAD: " << mad);
        DEBUG("Thresholds set to: [" << low << ", " << high << "]");

        size_t n = 0;
        double sum = 0.;
        double sum2 = 0.;
        DEBUG("Counting average");
        for (auto iter = hist_.begin(); iter != hist_.end(); ++iter) {
            if (iter->first < low || iter->first > high) {
                continue;
            }
            n += iter->second;
            sum += (double) iter->second * 1. * (double) iter->first;
            sum2 += (double)iter->second * 1. * (double)iter->first * (double)iter->first;
        }
        mean = sum / (double) n;
        delta = sqrt(sum2 / (double) n - mean * mean);

        low = mean - 5 * delta;
        high = mean + 5 * delta;

        DEBUG("Mean IS: " << mean);
        DEBUG("sd: " << delta);
        DEBUG("Thresholds set to: [" << low << ", " << high << "]");

        n = 0;
        sum = 0.;
        sum2 = 0.;
        for (auto iter = hist_.begin(); iter != hist_.end(); ++iter) {
            if (iter->first < low || iter->first > high) {
                continue;
            }
            n += iter->second;
            sum += (double) iter->second * 1. * (double) iter->first;
            sum2 += (double) iter->second * 1. * (double) iter->first * (double) iter->first;
        }
        mean = sum / (double) n;
        delta = sqrt(sum2 / (double) n - mean * mean);

        DEBUG("Mean IS: " << mean);
        DEBUG("sd: " << delta);

        size_t m = 0;

        DEBUG("Counting percentiles");
        //todo optimize
        size_t q[19];
        for (size_t i = 1; i < 20; ++i) {
            q[i - 1] = 5 * i;
        }
        for (auto iter = hist_.begin(); iter != hist_.end(); ++iter) {
            if (iter->first < low || iter->first > high) {
                continue;
            }
            size_t mm = m + iter->second;
            for (size_t i = 0; i < utils::array_size(q); i++) {
                size_t scaled_q_i((size_t) ((double) q[i] / 100. * (double) n));
                if (m < scaled_q_i && mm >= scaled_q_i) {
                    percentiles[q[i]] = iter->first;
                }
            }
            m = mm;
        }
    }

    void FindMedian(double& median, double& mad, HistType& histogram) const {
        DEBUG("Counting median and MAD");
        median = get_median(hist_);
        mad = get_mad(hist_, median);
        double low = median - 5. * 1.4826 * mad;
        double high = median + 5. * 1.4826 * mad;
        omnigraph::hist_crop(hist_, low, high, histogram);
        median = get_median(histogram);
        mad = get_mad(histogram, median);
    }

private:

    struct count_data {
      size_t total_;
      vector<size_t> arr_;
      count_data(): total_(0) {
      }
      count_data(size_t nthreads): total_(0), arr_(nthreads, 0) {
      }
      void inc(size_t i) {
        ++arr_[i];
      }
      void merge() {
        for (size_t i = 0; i < arr_.size(); ++i) {
          total_ += arr_[i];
        }
      }
    };

private:
    const conj_graph_pack& gp_;

    HistType hist_;
    vector<HistType*> tmp_hists_;

    count_data total_;
    count_data counted_;
    count_data negative_;

    size_t k_;
    size_t edge_length_threshold_;
    bool ignore_negative_;
};

}


#endif /* IS_COUNTER_HPP_ */
