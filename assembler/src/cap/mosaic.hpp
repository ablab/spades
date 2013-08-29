#include "standard_base.hpp"
#include "utils.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/sequence_reader.hpp"
#include "diff_masking.hpp"
#include "adt/bag.hpp"

namespace cap {

template<class gp_t>
class MosaicStructureAnalyzer {
    typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;

    Sequence genome_;
    const gp_t& gp_;
    const Graph& g_;
    vector<EdgeId> genome_path_;
    //coords are in k+1-mers
    multimap<EdgeId, Range> occurences_;

    size_t min_support_block_length_;
    size_t max_support_block_multiplicity_;
    size_t max_inter_block_length_;

    size_t curr_block_;

    void FillOccurences() {
        size_t cumm_coord = 0;
        for (EdgeId e : genome_path_) {
            occurences_.insert(
                    make_pair(e, Range(cumm_coord, cumm_coord + g_.length(e))));
            cumm_coord += g_.length(e);
        }
    }

//    void Construct() {
//        io::ReadStreamVector streams(
//                new io::CleanRCReaderWrapper<io::SingleRead>(
//                        new io::SequenceReader<io::SingleRead>(genome_,
//                                                               "genome")));
//        CapConstructGraph(k_, streams, gp_.g, gp_.index);
//        RefineGP(gp_, k_ / 2);
//    }

    size_t Multiplicity(EdgeId e) const {
        return get_all(occurences_, e).size();
    }

    bool CheckSupporting(EdgeId e) {
        size_t mult = Multiplicity(e);
        return mult > 1 && g_.length(e) > min_support_block_length_
                && mult < max_support_block_multiplicity_;
    }

    bool CheckSupporting() {
        return CheckSupporting(genome_path_[curr_block_]);
    }

    /*
     * returns number of k+1-mers to the next support block or -1u if it wasn't found
     * curr_block is updated as a side effect!
     */
    size_t MoveToNextSupportBlock() {
        size_t distance = 0;
        curr_block_++;
        while (curr_block_ < genome_path_.size() && !CheckSupporting()) {
            distance += g_.length(genome_path_[curr_block_]);
            curr_block_++;
        }
        return curr_block_ == genome_path_.size() ? -1u : distance;
    }

    struct MosaicInterval {
        vector<EdgeId> support_blocks;
        vector<EdgeId> all_blocks;
        Range path_position;
        Range genome_position;

        MosaicInterval(EdgeId e, size_t path_pos, size_t genome_pos_start,
                       size_t genome_pos_end)
                : path_position(path_pos, path_pos + 1),
                  genome_position(genome_pos_start, genome_pos_end) {
            support_blocks.push_back(e);
            all_blocks.push_back(e);
        }
    };

    class MosaicIntervalSet {
        const Graph& g_;
        vector<MosaicInterval> raw_intervals_;
        vector<MosaicInterval> final_intervals_;
        bag<string> struct_cnt_;

        string SupportBlockKey(const vector<EdgeId>& blocks) const {
            std::stringstream ss;
            string delim = "";
            for (EdgeId e : blocks) {
                ss << delim;
                ss << g_.int_id(e);
                delim = "$";
            }
            return ss.str();
        }

        template<class It1, class It2>
        It2 Find(It1 pattern_begin, It1 pattern_end, It2 target_begin,
                 It2 target_end) const {
            for (It2 it = target_begin;; ++it) {
                size_t i = 0;
                bool flag = true;
                for (It1 it2 = pattern_begin; it2 != pattern_end; ++it2) {
                    if (it + i == target_end)
                        return target_end;
                    if (*(it + i) != *it2) {
                        flag = false;
                        break;
                    }
                    ++i;
                }
                if (flag) {
                    return it;
                }
            }
            return target_end;
        }

        template<class Container1, class Container2>
        size_t Find(const Container1& pattern, const Container2& target) const {
            auto it = Find(pattern.begin(), pattern.end(), target.begin(),
                           target.end());
            if (it == target.end())
                return -1u;
            else
                return it - target.begin();
        }

        bool IsSubinterval(const MosaicInterval& a,
                           const MosaicInterval& b) const {
            return Find(a.support_blocks, b.support_blocks) != -1u;
        }

        bool IsSubintervalOfAnswer(const MosaicInterval& to_check) const {
            for (const MosaicInterval& interval : final_intervals_) {
                if (IsSubinterval(to_check, interval)) {
                    return true;
                }
            }
            return false;
        }

        void CountSubIntervals(const vector<EdgeId>& interval) {
            for (size_t i = 0; i < interval.size(); ++i) {
                for (size_t j = i; j < interval.size(); ++j) {
                    vector<EdgeId> sub_interval;
                    std::copy(interval.begin() + i, interval.begin() + j + 1,
                              std::back_inserter<vector<EdgeId>>(sub_interval));
                    struct_cnt_.put(SupportBlockKey(sub_interval));
                }
            }
        }

    public:
        MosaicIntervalSet(const Graph& g) : g_(g) {

        }

        const vector<MosaicInterval>& final_intervals() const {
            return final_intervals_;
        }

//        const bag<string>& struct_cnt() const {
//            return struct_cnt_;
//        }

        void ProcessInterval(const MosaicInterval& interval) {
            raw_intervals_.push_back(interval);
        }

        void FinalAnalysis() {
            std::sort(
                    raw_intervals_.begin(),
                    raw_intervals_.end(),
                    [](const MosaicInterval& a, const MosaicInterval& b) {
                        return a.support_blocks.size() > b.support_blocks.size();
                    });
            for (const MosaicInterval& interval : raw_intervals_) {
                CountSubIntervals(interval.support_blocks);
                if (!IsSubintervalOfAnswer(interval)) {
                    final_intervals_.push_back(interval);
                }
            }
        }

        string IdsConcat(const vector<EdgeId>& interval) const {
            std::stringstream ss;
            string delim = "";
            for (EdgeId e : interval) {
                ss << delim;
                ss << g_.int_id(e);
                delim = " ";
            }
            return ss.str();
        }

        void ReportSubIntervalCount(const vector<EdgeId>& interval) const {
            for (size_t i = 0; i < interval.size(); ++i) {
                for (size_t j = i; j < interval.size(); ++j) {
                    vector<EdgeId> sub_interval;
                    std::copy(interval.begin() + i, interval.begin() + j + 1,
                              std::back_inserter<vector<EdgeId>>(sub_interval));
                    std::cout << "Support multiplicity of " << IdsConcat(sub_interval) << " is "
                    << struct_cnt_.mult(SupportBlockKey(sub_interval)) << std::endl;
                }
            }
        }

    };

    void Report(const vector<EdgeId>& interval) const {
        using std::cout;
        using std::endl;
        for (EdgeId e : interval) {
            cout << g_.int_id(e);
            cout << " (length: " << g_.length(e);
            cout << ", mult: " << Multiplicity(e);
            cout << "); ";
        }
        cout << endl;
    }

    void Report(const MosaicInterval& interval) const {
        using std::cout;
        using std::endl;
        cout << "--------------------" << endl;
        cout << "Mosaic:" << endl;
        cout << "Genome pos: " << interval.genome_position << endl;
        cout << "Path pos: " << interval.path_position << endl;
        cout << "Support blocks:" << endl;
        Report(interval.support_blocks);
        cout << "All blocks:" << endl;
        Report(interval.all_blocks);
    }

    void Report(const MosaicIntervalSet& intervals) const {
        for (const MosaicInterval& interval : intervals.final_intervals()) {
            Report(interval);
            intervals.ReportSubIntervalCount(interval.support_blocks);
        }
    }

public:
    MosaicStructureAnalyzer(const gp_t& gp, const Sequence& genome,
                            size_t min_support_length, size_t max_support_mult,
                            size_t max_inter_length)
            : genome_(genome),
              gp_(gp),
//              gp_(k_, "tmp", genome_),
              g_(gp_.g),
              min_support_block_length_(min_support_length),
              max_support_block_multiplicity_(max_support_mult),
              max_inter_block_length_(max_inter_length),
              curr_block_(0) {
    }

    void Analyze() {
//        Construct();
        genome_path_ = MapperInstance(gp_)->MapSequence(genome_).simple_path()
                .sequence();
        FillOccurences();
        MosaicIntervalSet interval_set(g_);

        size_t genome_pos = MoveToNextSupportBlock();
        while (curr_block_ < genome_path_.size()) {
            EdgeId e = genome_path_[curr_block_];
            MosaicInterval interval(e, curr_block_, genome_pos,
                                    genome_pos + g_.length(e));
            genome_pos += g_.length(e);
            size_t dist;
            while ((dist = MoveToNextSupportBlock()) < max_inter_block_length_
                    && curr_block_ < genome_path_.size()) {
                e = genome_path_[curr_block_];
                genome_pos += (dist + g_.length(e));
                interval.genome_position.end_pos = genome_pos;
                std::copy(
                        genome_path_.begin() + interval.path_position.end_pos,
                        genome_path_.begin() + curr_block_ + 1,
                        std::back_insert_iterator<std::vector<EdgeId>>(
                                interval.all_blocks));
                interval.path_position.end_pos = (curr_block_ + 1);
                interval.support_blocks.push_back(e);
            }
            if (interval.support_blocks.size() > 1) {
                interval_set.ProcessInterval(interval);
            }
        }
        interval_set.FinalAnalysis();
        Report(interval_set);
    }

};

template<class gp_t>
void PerformMosaicAnalysis(const gp_t& gp, const Sequence& genome,
                        size_t min_support_length, size_t max_support_mult,
                        size_t max_inter_length) {
    MosaicStructureAnalyzer<gp_t> analyzer(gp, genome, min_support_length, max_support_mult, max_inter_length);
    analyzer.Analyze();
}

}
