#include "standard_base.hpp"
#include "utils.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/sequence_reader.hpp"
#include "diff_masking.hpp"
#include "adt/bag.hpp"

namespace cap {

namespace mosaic {

//todo temporary for correct eclipse highlight
using std::make_pair;

typedef ConjugateDeBruijnGraph Graph;
typedef Graph::EdgeId EdgeId;
typedef size_t Block;
typedef size_t Pos;

class BlockInfoProvider {
    const Graph& g_;

public:
    BlockInfoProvider(const Graph& g) : g_(g) {
    }

    size_t length(const Block& block) const {
        return g_.length(edge_id(block));
    }

    Sequence seq(const Block& block) const {
        return g_.EdgeNucls(edge_id(block));
    }

    EdgeId edge_id(const Block& block) const {
        return g_.int_ids().ReturnEdgeId(block);
    }
};

class GenomeBlockComposition {
    vector<Block> blocks_;
    multimap<Block, Pos> occurences_;
    map<Pos, Range> genome_coordinates_;

    bool CheckPath(const MappingPath<EdgeId>& mapping_path) const {
        for (Pos i = 1; i < mapping_path.size(); ++i)
            if (mapping_path[i].second.initial_range.start_pos != mapping_path[i - 1].second.initial_range.end_pos)
                return false;
        return true;
    }

public:
    GenomeBlockComposition(const Graph& g, const MappingPath<EdgeId>& mapping_path) {
        VERIFY(CheckPath(mapping_path));
        for (EdgeId e : mapping_path.simple_path().sequence()) {
            blocks_.push_back(g.int_id(e));
        }
        for (Pos i = 0; i < blocks_.size(); ++i) {
            genome_coordinates_.insert(make_pair(i, mapping_path[i].second.initial_range));
            occurences_.insert(make_pair(blocks_[i], i));
        }
    }

    Pos size() const {
        return blocks_.size();
    }

    const vector<Block>& blocks() const {
        return blocks_;
    }

    Block block(Pos pos) const {
        return blocks_[pos];
    }

    Range genome_coords(Pos pos) const {
        return get(genome_coordinates_, pos);
    }

    size_t multiplicity(Block block) const {
        return occurences_.count(block);
    }

    vector<Pos> occurences(Block block) const {
        return get_all(occurences_, block);
    }

    vector<Range> all_genome_coords(Block block) const {
        vector<Range> answer;
        for (size_t pos : occurences(block)) {
            answer.push_back(genome_coords(pos));
        }
        return answer;
    }

    Range genome_coords(Range pos_range) const {
        VERIFY(pos_range.end_pos > 0);
        return Range(genome_coords(pos_range.start_pos).start_pos, genome_coords(pos_range.end_pos - 1).end_pos);
    }

    size_t genome_span(Range pos_range) const {
        return genome_coords(pos_range).size();
    }
};

//todo maybe use ranges for latter parts of analysis
struct MosaicInterval {
    Range pos_range;
    vector<Pos> support_blocks;

    MosaicInterval(Pos pos)
            : pos_range(pos, pos + 1) {
        support_blocks.push_back(pos);
    }

    MosaicInterval(Range pos_range_, const vector<Pos>& support_blocks_)
            : pos_range(pos_range_), support_blocks(support_blocks_) {
    }

    MosaicInterval SubInterval(Range pos_range_) const {
        vector<Pos> sub_support_blocks;
        for (Pos pos : support_blocks) {
            if (pos >= pos_range_.start_pos && pos < pos_range_.end_pos) {
                sub_support_blocks.push_back(pos);
            }
        }
        return MosaicInterval(pos_range_, sub_support_blocks);
    }

    size_t support_size() const {
        return support_blocks.size();
    }

    bool operator<(const MosaicInterval &other) const {
      return pos_range < other.pos_range;
    }
};

template<class It1, class It2>
It2 Find(It1 pattern_begin, It1 pattern_end, It2 target_begin,
         It2 target_end) {
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
size_t Find(const Container1& pattern, const Container2& target) {
    auto it = Find(pattern.begin(), pattern.end(), target.begin(),
                   target.end());
    if (it == target.end())
        return -1u;
    else
        return it - target.begin();
}

class MosaicStructure {
    vector<Block> blocks_;
    vector<MosaicInterval> occurences_;

    vector<Block> PosToBlocks(const GenomeBlockComposition& block_composition, const vector<Pos>& poss) {
        vector<Block> answer;
        for (Pos pos : poss) {
            answer.push_back(block_composition.block(pos));
        }
        return answer;
    }

    //todo simplify after switching to Ranges
    vector<MosaicInterval> SubIntervals(size_t block_start, size_t block_end) const {
        vector<MosaicInterval> answer;
        for (const MosaicInterval& interval : occurences_) {
            answer.push_back(interval.SubInterval(Range(interval.support_blocks[block_start],
                                                                interval.support_blocks[block_end] + 1)));
        }
        return answer;
    }

public:
    explicit MosaicStructure(const vector<Block>& blocks)
    : blocks_(blocks) {
    }

    MosaicStructure(const vector<Block>& blocks, const MosaicInterval& interval) :
                        blocks_(blocks) {
        occurences_.push_back(interval);
    }

    MosaicStructure(const vector<Block>& blocks, const vector<MosaicInterval>& occurences) :
                        blocks_(blocks), occurences_(occurences) {
    }

    MosaicStructure(const GenomeBlockComposition& block_composition, const MosaicInterval& interval) :
                        blocks_(PosToBlocks(block_composition, interval.support_blocks)) {
        occurences_.push_back(interval);
    }

    void AddOccurence(const MosaicInterval& interval) {
//        VERIFY(blocks_ == PosToBlocks(interval.support_blocks));
        occurences_.push_back(interval);
    }

    size_t mult() const {
        return occurences_.size();
    }

    const vector<MosaicInterval>& occurences() const {
        return occurences_;
    }

    const vector<Range> occurence_ranges() const {
        vector<Range> answer;
        for (const auto& interval : occurences_) {
            answer.push_back(interval.pos_range);
        }
        return answer;
    }

    const vector<Block>& blocks() const {
        return blocks_;
    }

    size_t block_size() const {
        return blocks_.size();
    }

    //block end incl
    MosaicStructure SubMosaic(size_t block_start, size_t block_end) const {
        VERIFY(block_start < blocks_.size());
        VERIFY(block_end < blocks_.size());
        VERIFY(block_start <= block_end);

        const MosaicInterval& interval = occurences_.front();
        return MosaicStructure(vector<Block>(blocks_.begin() + block_start, blocks_.begin() + block_end + 1),
            SubIntervals(block_start, block_end));
    }

    string Fingerprint() const {
        std::stringstream ss;
        string delim = "";
        for (Block block : blocks_) {
            ss << delim;
            ss << block;
            delim = " ";
        }
        return ss.str();
    }

    bool SameBlocks(const MosaicStructure& that) const {
        return blocks_ == that.blocks_;
    }

    void Merge(const MosaicStructure& that) {
        VERIFY(SameBlocks(that));
        push_back_all(occurences_, that.occurences_);
    }

    bool IsContainedIn(const MosaicStructure& that) const {
        return Find(blocks_, that.blocks_) != -1u;
    }

};

ostream& operator << (ostream& out, const MosaicStructure& structure) {
    return out << "Mosaic. size=" << structure.block_size() << " blocks=" << structure.blocks();
}

class MosaicPrintHelper {
    const GenomeBlockComposition& block_composition_;
    const BlockInfoProvider& block_info_;
    const multimap<string, size_t>& different_irred_presence_;
    const multimap<string, Range>& all_substruct_pos_;
    ostream& out_;

public:
    MosaicPrintHelper(const GenomeBlockComposition& block_composition,
                      const BlockInfoProvider& block_info,
                      const multimap<string, size_t>& different_irred_presence,
                      const multimap<string, Range>& all_substruct_pos,
                      ostream& out)
    : block_composition_(block_composition), block_info_(block_info), different_irred_presence_(different_irred_presence),
      all_substruct_pos_(all_substruct_pos), out_(out) {
    }

    void BlockInfo(Block b) {
        out_ << b << " (" << block_info_.length(b) << ")";
    }

    double AvgSpan(const MosaicStructure& mosaic) {
        double avg = 0.;
        for (Range r : mosaic.occurence_ranges()) {
            avg += r.size();
        }
        return avg / mosaic.mult();
    }

    double AvgGenomicSpan(const MosaicStructure& mosaic) {
        double avg = 0.;
        for (Range r : mosaic.occurence_ranges()) {
            Range genomic_range = block_composition_.genome_coords(r);
            avg += genomic_range.size();
        }
        return avg / mosaic.mult();
    }

    vector<double> AvgGenomicInterLengths(const MosaicStructure& mosaic) {
        vector<double> answer(mosaic.block_size() - 1, 0.);
        for (const MosaicInterval& interval : mosaic.occurences()) {
            for (size_t i = 1; i < interval.support_size(); ++i) {
                answer[i - 1] += block_composition_.genome_coords(Range(interval.support_blocks[i - 1] + 1,
                                                                        interval.support_blocks[i])).size();
            }
        }
        for (size_t i = 0; i < answer.size(); ++i) {
            answer[i] /= mosaic.mult();
        }
        return answer;
    }

    size_t TotalBlockLength(const MosaicStructure& mosaic) {
        size_t block_length = 0;
        for (Block b : mosaic.blocks()) {
            block_length += block_info_.length(b);
        }
        return block_length;
    }

    void ReportBasicInfo(const MosaicStructure& mosaic) {
        out_ << "Support block cnt = " << mosaic.block_size();
        out_ << "; Total support block length = " << TotalBlockLength(mosaic);
        out_ << "; Full mosaic multiplicity = " << mosaic.mult();
        out_ << "; Avg block span = " << AvgSpan(mosaic);
        out_ << "; Avg genome span = " << AvgGenomicSpan(mosaic);
        out_ << endl;
        if (mosaic.occurences().size() > 1) {
            out_ << "WARN! Full mosaic multiplicity = " << mosaic.mult();
            out_ << endl;
        }
        out_ << "Structure: ";
        vector<double> inter_lengths = AvgGenomicInterLengths(mosaic);
        BlockInfo(mosaic.blocks().front());
        for (size_t i = 1; i < mosaic.block_size(); ++i) {
            out_ << " |...";
            out_ << inter_lengths[i - 1];
            out_ << "...| ";
            BlockInfo(mosaic.blocks()[i]);
        }
        out_ << endl;
    }

    bool Reported(Range range, const set<Range>& reported) const {
        for (auto r : reported)
            if (r.contains(range))
                return true;
        return false;
    }

    vector<Range> FilterReported(const vector<Range>& ranges, const set<Range>& reported) const {
        vector<Range> answer;
        for (auto range : ranges)
            if (!Reported(range, reported))
                answer.push_back(range);
        return answer;
    }

    void ReportSubMosaic(const MosaicStructure& mosaic, const vector<Range>& ranges) {
        string finger = mosaic.Fingerprint();
        out_ << "------" << endl;
        out_ << "Sub_mosaic. Block cnt = " << mosaic.block_size() << endl;
        out_ << "Blocks " << finger;
        out_ << " ; Found in " << get_all(different_irred_presence_, finger).size() << " different irreducible mosaics";
        string delim = " (";
        for (size_t idx : get_all(different_irred_presence_, finger)) {
            out_ << delim;
            out_ << idx;
            delim = ", ";
        }
        out_ << ")" << endl;

        delim = "";
        out_ << "Ranges: ";
        for (Range r : ranges) {
            out_ << delim;
            out_ << block_composition_.genome_coords(r);
            out_ << " (Pos: ";
            out_ << r;
            out_ << ")";
            delim = "; ";
        }
        out_ << endl;
    }

    void ReportSubMosaics(const MosaicStructure& mosaic) {
        out_ << "Sub_mosaics" << endl;
        set<Range> reported_ranges;
        for (size_t d = mosaic.block_size(); d > 0; --d) {
            for (size_t i = 0; i + d < mosaic.block_size(); ++i) {
                MosaicStructure sub_mosaic = mosaic.SubMosaic(i, i + d);
                vector<Range> ranges = FilterReported(get_all(all_substruct_pos_, sub_mosaic.Fingerprint()), reported_ranges);
                if (!ranges.empty()) {
                    ReportSubMosaic(sub_mosaic, ranges);
                    insert_all(reported_ranges, ranges);
                }
            }
        }
    }

    void ReportIrredMosaic(const MosaicStructure& mosaic, size_t idx) {
        out_ << "Irreducible Mosaic " << idx << endl;
        ReportBasicInfo(mosaic);
        ReportSubMosaics(mosaic);
        out_ << "..................................." << endl;
    }

};

class MosaicStructureSet {
    const GenomeBlockComposition& block_composition_;
    vector<MosaicInterval> raw_intervals_;

    vector<MosaicStructure> irreducible_structures_;
    //fixme currently not used in any way!!!
    map<string, MosaicStructure> nested_structures_;
    multimap<string, size_t> different_irred_presence_;
    multimap<string, Range> all_substruct_pos_;

    void IndexSubIntervals(const MosaicStructure& mosaic) {
        for (size_t i = 0; i < mosaic.block_size(); ++i) {
            for (size_t j = i; j < mosaic.block_size(); ++j) {
                MosaicStructure sub_mosaic = mosaic.SubMosaic(i, j);
                all_substruct_pos_.insert(make_pair(sub_mosaic.Fingerprint(), sub_mosaic.occurences().front().pos_range));
            }
        }
    }

    void CountDifferentIrred(const MosaicStructure& mosaic, size_t idx) {
        for (size_t i = 0; i < mosaic.block_size(); ++i) {
            for (size_t j = i; j < mosaic.block_size(); ++j) {
                MosaicStructure sub_mosaic = mosaic.SubMosaic(i, j);
                different_irred_presence_.insert(make_pair(sub_mosaic.Fingerprint(), idx));
            }
        }
    }

    void CountDifferentIrred() {
        for (size_t i = 0; i < irreducible_structures_.size(); ++i) {
            CountDifferentIrred(irreducible_structures_[i], i);
        }
    }

    bool CountNested(const MosaicStructure& nested, const MosaicStructure& /*outer*/) {
        string finger = nested.Fingerprint();
        if (nested_structures_.count(finger) == 0) {
            nested_structures_.insert(make_pair(finger, nested));
            return true;
        } else {
            get(nested_structures_, finger).Merge(nested);
            return false;
        }
    }

    bool AnalyzeStructure(const MosaicStructure& mosaic) {
//        cout << "Analyzing " << mosaic << endl;
//        cout << irreducible_structures_.size() << endl;
        for (auto& irred_struct : irreducible_structures_) {
//            cout << "Irred " << irred_struct << endl;
            if (irred_struct.SameBlocks(mosaic)) {
//                cout << "same" << endl;
                irred_struct.Merge(mosaic);
                return false;
            }
            if (mosaic.IsContainedIn(irred_struct)) {
//                cout << "contained" << endl;
                CountNested(mosaic, irred_struct);
                return false;
            }
        }
//        cout << "new_irred" << endl;
        irreducible_structures_.push_back(mosaic);
        return true;
    }

public:
    MosaicStructureSet(const GenomeBlockComposition& block_composition) : block_composition_(block_composition) {

    }

    void ProcessInterval(const MosaicInterval& interval) {
        //todo fix two conversions
        IndexSubIntervals(MosaicStructure(block_composition_, interval));
        if (interval.support_size() > 1) {
            raw_intervals_.push_back(interval);
        }
    }

    void Analysis() {
        INFO("Sorting raw intervals");
        std::sort(
                raw_intervals_.begin(),
                raw_intervals_.end(),
                [](const MosaicInterval& a, const MosaicInterval& b) {
                    return a.support_blocks.size() > b.support_blocks.size();
                });
        INFO("Analyzing sorted intervals");
        for (const MosaicInterval& interval : raw_intervals_) {
            AnalyzeStructure(MosaicStructure(block_composition_, interval));
        }
        INFO("Counting distinct irreducible");
        CountDifferentIrred();
    }

    void Report(BlockInfoProvider& block_info, ostream& out) {
        MosaicPrintHelper printer(block_composition_, block_info,
                                  different_irred_presence_, all_substruct_pos_, out);
        for (size_t i = 0; i < irreducible_structures_.size(); ++i) {
            printer.ReportIrredMosaic(irreducible_structures_[i], i);
        }
    }

};

template<class gp_t>
class MosaicStructureAnalyzer {
    typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;

    GenomeBlockComposition block_composition_;
    BlockInfoProvider block_info_;

    size_t min_support_block_length_;
    size_t max_support_block_multiplicity_;
    size_t max_inter_block_length_;
    ostream& out_;

    Pos curr_pos_;

    bool CheckSupporting(Block b) {
        size_t mult = block_composition_.multiplicity(b);
        return mult > 1 && block_info_.length(b) > min_support_block_length_
                && mult < max_support_block_multiplicity_;
    }

    bool CheckSupporting() {
        return CheckSupporting(block_composition_.block(curr_pos_));
    }

    /*
     * returns distance to the next support block or -1u if it wasn't found
     * curr_block is updated as a side effect!
     */
    size_t MoveToNextSupportBlock() {
        Pos init_pos = curr_pos_;
        curr_pos_++;
        while (curr_pos_ < block_composition_.size() && !CheckSupporting()) {
            curr_pos_++;
        }
        return curr_pos_ != block_composition_.size()
                ? block_composition_.genome_coords(curr_pos_).start_pos
                - block_composition_.genome_coords(init_pos).end_pos
                : -1u;
    }

public:
    MosaicStructureAnalyzer(const gp_t& gp, const Sequence& genome,
                            size_t min_support_length, size_t max_support_mult,
                            size_t max_inter_length,
                            ostream& out)
            : block_composition_(gp.g, MapperInstance(gp)->MapSequence(genome)),
              block_info_(gp.g),
              min_support_block_length_(min_support_length),
              max_support_block_multiplicity_(max_support_mult),
              max_inter_block_length_(max_inter_length),
              out_(out),
              curr_pos_(0) {
    }

    void Analyze() {
        MosaicStructureSet interval_set(block_composition_);

        INFO("Collecting mosaic intervals");
        MoveToNextSupportBlock();
        while (curr_pos_ < block_composition_.size()) {
            MosaicInterval interval(curr_pos_);
            while (MoveToNextSupportBlock() < max_inter_block_length_) {
                interval.pos_range.end_pos = (curr_pos_ + 1);
                interval.support_blocks.push_back(curr_pos_);
            }
            interval_set.ProcessInterval(interval);
        }
        INFO("Analyzing intervals and forming mosaic structures");
        interval_set.Analysis();
        INFO("Reporting mosaic structures");
        interval_set.Report(block_info_, out_);
    }

};

template<class gp_t>
void PerformMosaicAnalysis(const gp_t& gp, const Sequence& genome,
                        size_t min_support_length, size_t max_support_mult,
                        size_t max_inter_length, ostream& out) {
    MosaicStructureAnalyzer<gp_t> analyzer(gp, genome, min_support_length, max_support_mult, max_inter_length, out);
    analyzer.Analyze();
}

}
}
