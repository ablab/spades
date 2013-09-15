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
public:
    GenomeBlockComposition(const Graph& g, const MappingPath<EdgeId>& mapping_path) {
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
        VERIFY(pos_range.end_pos > pos_range.start_pos);
        return Range(genome_coords(pos_range.start_pos).start_pos, genome_coords(pos_range.end_pos - 1).end_pos);
    }

    size_t genome_span(Range pos_range) const {
        return genome_coords(pos_range).size();
    }
};

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

public:
    explicit MosaicStructure(const vector<Block>& blocks)
    : blocks_(blocks) {
    }

    MosaicStructure(const vector<Block>& blocks, const MosaicInterval& interval) :
                        blocks_(blocks) {
        occurences_.push_back(interval);
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

//    const vector<Range>& occurence_ranges() const {
//        vector<Range> answer;
//        for (const auto& interval : occurences_) {
//            answer.push_back(interval.pos_range);//todo
//        }
//        return occurences_;
//    }

    size_t block_size() const {
        return blocks_.size();
    }

    //block end incl
    MosaicStructure SubMosaic(size_t block_start, size_t block_end) const {
        VERIFY(occurences_.size() == 1);
        VERIFY(block_start < blocks_.size());
        VERIFY(block_end < blocks_.size());
        VERIFY(block_start <= block_end);

        const MosaicInterval& interval = occurences_.front();
        return MosaicStructure(vector<Block>(blocks_.begin() + block_start, blocks_.begin() + block_end + 1),
                               interval.SubInterval(Range(interval.support_blocks[block_start],
                                                    interval.support_blocks[block_end] + 1)));
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

class MosaicStructureSet {
    const GenomeBlockComposition& block_composition_;
    vector<MosaicInterval> raw_intervals_;

    vector<MosaicStructure> irreducible_structures_;
//        bag<string> struct_cnt_;
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
        if (nested_structures_.count(finger) > 0) {
            nested_structures_.insert(make_pair(finger, nested));
            return true;
        } else {
            get(nested_structures_, finger).Merge(nested);
            return false;
        }
    }

    bool AnalyzeStructure(const MosaicStructure& mosaic) {
        for (auto& irred_struct : irreducible_structures_) {
            if (irred_struct.SameBlocks(mosaic)) {
                irred_struct.Merge(mosaic);
                return false;
            }
            if (mosaic.IsContainedIn(irred_struct)) {
                CountNested(mosaic, irred_struct);
                return false;
            }
        }
        irreducible_structures_.push_back(mosaic);
        return true;
    }

    void ReportBasicInfo(const MosaicStructure& mosaic, BlockInfoProvider& block_info, ostream& out) const {
    }

    void ReportSubMosaics(const MosaicStructure& mosaic, BlockInfoProvider& block_info, ostream& out) const {
        for (size_t d = mosaic.block_size() - 1; d > 0; --d) {
            for (size_t i = 0; i + d < mosaic.block_size(); ++i) {
                MosaicStructure sub_mosaic = mosaic.SubMosaic(i, i + d);
                get_all(all_substruct_pos_, sub_mosaic.Fingerprint());
            }
        }
    }

    void ReportIrredMosaic(const MosaicStructure& mosaic, BlockInfoProvider& block_info, ostream& out) const {
        set<Range> reported_ranges;
        insert_all(reported_ranges, mosaic.occurence_ranges());
        ReportBasicInfo(mosaic, block_info, out);
        ReportSubMosaics(mosaic, block_info, out);
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
        std::sort(
                raw_intervals_.begin(),
                raw_intervals_.end(),
                [](const MosaicInterval& a, const MosaicInterval& b) {
                    return a.support_blocks.size() > b.support_blocks.size();
                });
        for (const MosaicInterval& interval : raw_intervals_) {
            AnalyzeStructure(MosaicStructure(block_composition_, interval));
        }
        CountDifferentIrred();
    }

    void Report(BlockInfoProvider& block_info, ostream& out) {
        for (size_t i = 0; i < irreducible_structures_.size(); ++i) {
            out << "Irreducible Mosaic " << i << endl;
            ReportIrredMosaic(irreducible_structures_[i], block_info, out);
        }
    }

//    void ReportSubIntervalCount(const vector<EdgeId>& interval, ostream& out) const {
//        for (size_t i = 0; i < interval.size(); ++i) {
//            for (size_t j = i; j < interval.size(); ++j) {
//                vector<EdgeId> sub_interval;
//                std::copy(interval.begin() + i, interval.begin() + j + 1,
//                          std::back_inserter<vector<EdgeId>>(sub_interval));
//                out << "Support multiplicity of " << IdsConcat(sub_interval) << " is "
//                << struct_cnt_.mult(SupportBlockKey(sub_interval)) << std::endl;
//            }
//        }
//    }

};

template<class gp_t>
class MosaicStructureAnalyzer {
    typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;

//    Sequence genome_;
//    const gp_t& gp_;
//    const Graph& g_;
//    vector<EdgeId> genome_path_;
    //coords are in k+1-mers
//    multimap<EdgeId, Range> occurences_;

    GenomeBlockComposition block_composition_;
    BlockInfoProvider block_info_;

    size_t min_support_block_length_;
    size_t max_support_block_multiplicity_;
    size_t max_inter_block_length_;
    ostream& out_;

    Pos curr_pos_;

//    void FillOccurences() {
//        size_t cumm_coord = 0;
//        for (EdgeId e : genome_path_) {
//            occurences_.insert(
//                    make_pair(e, Range(cumm_coord, cumm_coord + g_.length(e))));
//            cumm_coord += g_.length(e);
//        }
//    }

//    void Construct() {
//        io::ReadStreamVector streams(
//                new io::CleanRCReaderWrapper<io::SingleRead>(
//                        new io::SequenceReader<io::SingleRead>(genome_,
//                                                               "genome")));
//        CapConstructGraph(k_, streams, gp_.g, gp_.index);
//        RefineGP(gp_, k_ / 2);
//    }

//    size_t Multiplicity(EdgeId e) const {
//        return get_all(occurences_, e).size();
//    }

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

//    void Report(const vector<EdgeId>& interval) const {
//        using std::endl;
//        for (EdgeId e : interval) {
//            out_ << g_.int_id(e);
//            out_ << " (length: " << g_.length(e);
//            out_ << ", mult: " << Multiplicity(e);
//            out_ << "); ";
//        }
//        out_ << endl;
//    }
//
//    void Report(const MosaicInterval& interval) const {
//        using std::endl;
//        out_ << "--------------------" << endl;
//        out_ << "Mosaic:" << endl;
//        out_ << "Genome pos: " << interval.genome_position << endl;
//        out_ << "Path pos: " << interval.path_position << endl;
//        out_ << "Support blocks:" << endl;
//        Report(interval.support_blocks);
//        out_ << "All blocks:" << endl;
//        Report(interval.all_blocks);
//    }
//
//    void Report(const MosaicIntervalSet& intervals) const {
//        for (const MosaicInterval& interval : intervals.final_intervals()) {
//            Report(interval);
//            intervals.ReportSubIntervalCount(interval.support_blocks, out_);
//        }
//    }

public:
    MosaicStructureAnalyzer(const gp_t& gp, const Sequence& genome,
                            size_t min_support_length, size_t max_support_mult,
                            size_t max_inter_length,
                            ostream& out)
            : /*genome_(genome),
              gp_(gp),
//              gp_(k_, "tmp", genome_),
              g_(gp_.g),*/
              block_composition_(gp.g, MapperInstance(gp)->MapSequence(genome)),
              block_info_(gp.g),
              min_support_block_length_(min_support_length),
              max_support_block_multiplicity_(max_support_mult),
              max_inter_block_length_(max_inter_length),
              out_(out),
              curr_pos_(0) {
    }

    void Analyze() {
//        Construct();
//        genome_path_ = MapperInstance(gp_)->MapSequence(genome_).simple_path()
//                .sequence();
//        FillOccurences();
        MosaicStructureSet interval_set(block_composition_);

        MoveToNextSupportBlock();
        while (curr_pos_ < block_composition_.size()) {
            MosaicInterval interval(curr_pos_);
            while (MoveToNextSupportBlock() < max_inter_block_length_) {
                interval.pos_range.end_pos = (curr_pos_ + 1);
                interval.support_blocks.push_back(curr_pos_);
            }
            interval_set.ProcessInterval(interval);
        }
        interval_set.Analysis();
        interval_set.Report(out_);
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
