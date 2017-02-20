#pragma once

namespace debruijn_graph {
using omnigraph::Path;
using omnigraph::MappingPath;
using omnigraph::MappingRange;

inline double PairedReadCountWeight(const std::pair<EdgeId, EdgeId>&,
                                    const MappingRange&, const MappingRange&) {
    return 1.;
}

inline double KmerCountProductWeight(const std::pair<EdgeId, EdgeId>&,
                                     const MappingRange& mr1, const MappingRange& mr2) {
    return (double)(mr1.initial_range.size() * mr2.initial_range.size());
}

class WeightDEWrapper {
private:

    vector<double> new_hist;
    int left_x;
    int insert_size;

    void ExtendLinear(const std::map<int, size_t> & hist) {
        size_t sum_weight = 0;

        for (auto iter = hist.begin(); iter != hist.end(); ++iter)
            sum_weight += iter->second;
        DEBUG(sum_weight);

        VERIFY(hist.size() > 0);
        auto iter = hist.begin();

        left_x = iter->first;

        int prev = iter->first;
        size_t prev_val = iter->second;

        new_hist.push_back((double)prev_val / (double)sum_weight);
        ++iter;

        for (; iter != hist.end(); ++iter) {
            int x = iter->first;
            size_t y = iter->second;
            double tan = ((double)y - (double)prev_val) / (x - prev);

            VERIFY(prev < x);
            for (int i = prev + 1; i <= x; ++i) {
                new_hist.push_back(((double)prev_val + tan * (i - prev)) / (double)sum_weight);
            }
            prev = x;
            prev_val = y;
            DEBUG("hist " << x << " " << y);
        }
    }

public:
    WeightDEWrapper(const map<int, size_t>& hist, double IS) {
        DEBUG("WeightDEWrapper " << IS);
        insert_size = (int) IS;
        DEBUG("Extending linear");
        ExtendLinear(hist);
    }

    ~WeightDEWrapper() {
    }


    double CountWeight(int x) const {
        int xx = insert_size - left_x + x - 1;

        if (!(xx >= 0 && xx < (int) new_hist.size())) return 0.;
        VERIFY(math::le(new_hist[xx], 1.));
        return 1000. * new_hist[xx];
    }
};

inline double UnityFunction(int /*x*/) {
    return 1.;
}
}
