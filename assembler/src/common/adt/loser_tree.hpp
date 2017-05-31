#pragma once

#include "iterator_range.hpp"
#include <vector>

namespace adt {

template<typename IntegerType>
IntegerType ilog2(IntegerType x) {
    IntegerType lg = 0;
    while (x >= 256) {
        x >>= 8;
        lg += 8;
    }
    while (x >>= 1) lg += 1;

    return lg;
}

template<typename IntegerType>
IntegerType ilog2ceil(IntegerType x) {
    return ilog2(x - 1) + 1;
}

template<class It, class Cmp>
class loser_tree {
    typedef typename std::iterator_traits<It>::value_type value_type;

    size_t log_k_;
    size_t k_;
    std::vector<size_t> entry_;
    Cmp inner_cmp_;

    bool cmp(const adt::iterator_range<It> &a, const adt::iterator_range<It> &b) const {
        // Emulate sentinels
        if (b.end() == b.begin())
            return true;
        if (a.end() == a.begin())
            return false;

        return inner_cmp_(*a.begin(), *b.begin());
    }

    size_t init_winner(size_t root) {
        if (root >= k_)
            return root - k_;

        size_t left = init_winner(2 * root);
        size_t right = init_winner(2 * root + 1);
        if (cmp(runs_[left], runs_[right])) {
            entry_[root] = right;
            return left;
        } else {
            entry_[root] = left;
            return right;
        }
    }

public:
    loser_tree(const std::vector<adt::iterator_range<It>> &runs,
               Cmp inner_cmp = Cmp())
            : inner_cmp_(inner_cmp), runs_(runs) {
        log_k_ = ilog2ceil(runs.size());
        k_ = (size_t(1) << log_k_);

        // fprintf(stderr, "k: %zu, logK: %zu, nruns: %zu\n", k_, log_k_, runs.size());

        entry_.resize(2 * k_);
        for (size_t i = 0; i < k_; ++i)
            entry_[k_ + i] = i;

        // Insert sentinels
        for (size_t i = runs.size(); i < k_; ++i)
            runs_.emplace_back(adt::make_range(runs_[0].end(), runs_[0].end()));

        // Populate tree
        entry_[0] = init_winner(1);

        // for (const auto &entry : entry_)
        //    fprintf(stderr, "%zu, ", entry);
        // fprintf(stderr, "\n");
    }

    size_t replay(size_t winner_index) {
        auto &winner = runs_[winner_index];
        if (winner.begin() == winner.end())
            return winner_index;

        winner = adt::make_range(std::next(winner.begin()), winner.end());
        for (size_t i = (winner_index + k_) >> 1; i > 0; i >>= 1)
            if (cmp(runs_[entry_[i]], runs_[winner_index]))
                std::swap(entry_[i], winner_index);

        return winner_index;
    }

    bool empty() const {
        size_t winner_index = entry_[0];
        const auto &winner = runs_[winner_index];
        return (winner.begin() == winner.end());
    }


    template<class It2>
    size_t multi_merge(It2 out, size_t amount = -1ULL) {
        size_t cnt = 0;
        size_t winner_index = entry_[0];

        for (cnt = 0; cnt < amount; ++cnt) {
            auto &winner = runs_[winner_index];
            if (winner.begin() == winner.end())
                break;

            *out++ = *winner.begin();

            winner_index = replay(winner_index);
        }

        entry_[0] = winner_index;

        return cnt;
    }

    value_type pop() {
        size_t winner_index = entry_[0];
        value_type res = *runs_[winner_index].begin();
        entry_[0] = replay(winner_index);

        return res;
    }


private:
    std::vector<adt::iterator_range<It>> runs_;
};

} //adt