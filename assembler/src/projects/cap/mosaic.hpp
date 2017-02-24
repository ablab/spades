//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/standard_base.hpp"
#include "io/reads/rc_reader_wrapper.hpp"
#include "io/reads/sequence_reader.hpp"
#include "diff_masking.hpp"
#include "common/adt/bag.hpp"
#include "io/reads/vector_reader.hpp"
#include "visualization/graph_colorer.hpp"

namespace cap {

namespace mosaic {

/*
 * todo
 * 1. Somehow report information on single edges
 * 2. Show multiplicity information
 */

//todo temporary for correct eclipse highlight
using std::make_pair;

typedef ConjugateDeBruijnGraph Graph;
typedef Graph::EdgeId EdgeId;

//block and its conjugate are represented by different natural
//which are not related in any trivial way
typedef size_t Block;

typedef size_t Pos;
typedef pair<Pos, bool> StrandPos;
typedef pair<Range, bool> StrandRange;

class BlockInfoProvider {
    const Graph& g_;
    omnigraph::GraphElementFinder<Graph> finder_;

public:
    BlockInfoProvider(const Graph& g) : g_(g), finder_(g) {
        finder_.Init();
    }

    Block conjugate(const Block& block) const {
        return g_.int_id(g_.conjugate(edge_id(block)));
    }

    size_t length(const Block& block) const {
        return g_.length(edge_id(block));
    }

    Sequence seq(const Block& block) const {
        return g_.EdgeNucls(edge_id(block));
    }

    EdgeId edge_id(const Block& block) const {
        return finder_.ReturnEdgeId(block);
    }

    const Graph& g() const {
        return g_;
    }
};

class GenomeBlockComposition {
    vector<Block> blocks_;
    multimap<Block, StrandPos> occurences_;
    map<Pos, Range> genome_coordinates_;
    //in k-mers!!!
    const size_t genome_length_;
    const BlockInfoProvider& block_info_;

    bool CheckPath(const MappingPath<EdgeId>& mapping_path) const {
        for (Pos i = 1; i < mapping_path.size(); ++i)
            if (mapping_path[i].second.initial_range.start_pos != mapping_path[i - 1].second.initial_range.end_pos)
                return false;
        return true;
    }

    vector<StrandPos> AddStrand(vector<Pos> poss, bool strand) const {
        vector<StrandPos> answer;
        for (auto pos : poss) {
            answer.push_back(make_pair(pos, strand));
        }
        return answer;
    }

public:
    GenomeBlockComposition(const Graph& g, const MappingPath<EdgeId>& mapping_path, size_t genome_length, const BlockInfoProvider& block_info)
    : genome_length_(genome_length), block_info_(block_info) {
        VERIFY(CheckPath(mapping_path));
        for (EdgeId e : mapping_path.simple_path()) {
            blocks_.push_back(g.int_id(e));
        }
        for (Pos i = 0; i < blocks_.size(); ++i) {
            genome_coordinates_.insert(make_pair(i, mapping_path[i].second.initial_range));
            StrandPos strand_pos = make_pair(i, true);
            occurences_.insert(make_pair(blocks_[i], strand_pos));
            occurences_.insert(make_pair(ConjBlock(blocks_[i]), ConjStrandPos(strand_pos)));
        }
    }

    Block ConjBlock(Block b) const {
        return block_info_.conjugate(b);
    }

    const BlockInfoProvider& block_info() const {
        return block_info_;
    }

    Pos ConjPos(Pos pos) const {
        return size() - pos - 1;
    }

    StrandPos ConjStrandPos(StrandPos pos) const {
        return make_pair(ConjPos(pos.first), !pos.second);
    }

    Range ConjCoordsRange(Range r) const {
        return Range(genome_length_ - r.end_pos, genome_length_ - r.start_pos);
    }

    StrandRange ConjStrandRange(StrandRange r) const {
        return make_pair(Range(size() - r.first.end_pos, size() - r.first.start_pos), !r.second);
    }

    Pos size() const {
        return blocks_.size();
    }

    const vector<Block>& blocks() const {
        return blocks_;
    }

    Block block(StrandPos pos) const {
        return block(pos.first, pos.second);
    }

    Block block(Pos pos, bool strand) const {
        if (strand)
            return blocks_[pos];
        else
            return block_info_.conjugate(blocks_[ConjPos(pos)]);
    }

    Range genome_coords(StrandPos pos) const {
        return genome_coords(pos.first, pos.second);
    }

    Range genome_coords(Pos pos, bool strand) const {
        if (strand)
            return get(genome_coordinates_, pos);
        else
            return ConjCoordsRange(get(genome_coordinates_, ConjPos(pos)));
    }

    size_t multiplicity(Block block) const {
        return occurences_.count(block);
    }

    vector<StrandPos> occurences(Block block) const {
        return get_all(occurences_, block);
    }

//    vector<Range> all_genome_coords(Block block) const {
//        vector<Range> answer;
//        for (size_t pos : occurences(block)) {
//            answer.push_back(genome_coords(pos));
//        }
//        return answer;
//    }

    Range genome_coords(StrandRange pos_range) const {
        return genome_coords(pos_range.first, pos_range.second);
    }

    Range genome_coords(Range pos_range, bool strand) const {
        VERIFY(pos_range.end_pos > 0);
        return Range(genome_coords(pos_range.start_pos, strand).start_pos, genome_coords(pos_range.end_pos - 1, strand).end_pos);
    }

    size_t genome_span(Range pos_range, bool strand) const {
        return genome_coords(pos_range, strand).size();
    }
};

//todo maybe use ranges for latter parts of analysis
struct MosaicInterval {
    bool strand;
    Range pos_range;
    vector<Pos> support_blocks;

    MosaicInterval(bool strand_, Pos pos)
            : strand(strand_), pos_range(pos, pos + 1) {
        support_blocks.push_back(pos);
    }

    MosaicInterval(bool strand_, Range pos_range_, const vector<Pos>& support_blocks_)
            : strand(strand_), pos_range(pos_range_), support_blocks(support_blocks_) {
    }

    MosaicInterval SubInterval(Range pos_range_) const {
        vector<Pos> sub_support_blocks;
        for (Pos pos : support_blocks) {
            if (pos >= pos_range_.start_pos && pos < pos_range_.end_pos) {
                sub_support_blocks.push_back(pos);
            }
        }
        return MosaicInterval(strand, pos_range_, sub_support_blocks);
    }

    StrandRange strand_range() const {
        return make_pair(pos_range, strand);
    }

    size_t support_size() const {
        return support_blocks.size();
    }

//    bool operator<(const MosaicInterval &other) const {
//      return pos_range < other.pos_range;
//    }
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
//    vector<MosaicInterval> occurences_;

    vector<Block> ToBlocks(const GenomeBlockComposition& block_composition, const MosaicInterval& interval) {
        vector<Block> answer;
        for (Pos pos : interval.support_blocks) {
            answer.push_back(block_composition.block(pos, interval.strand));
        }
        return answer;
    }

    //todo simplify after switching to Ranges
//    vector<MosaicInterval> SubIntervals(size_t block_start, size_t block_end) const {
        vector<MosaicInterval> answer;
//        for (const MosaicInterval& interval : occurences_) {
//            answer.push_back(interval.SubInterval(Range(interval.support_blocks[block_start],
//                                                                interval.support_blocks[block_end] + 1)));
//        }
//        return answer;
//    }

public:
    explicit MosaicStructure(const vector<Block>& blocks)
    : blocks_(blocks) {
    }

    MosaicStructure(const GenomeBlockComposition& block_composition, const MosaicInterval& interval) :
                        blocks_(ToBlocks(block_composition, interval)) {
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

        return MosaicStructure(vector<Block>(blocks_.begin() + block_start, blocks_.begin() + block_end + 1));
    }

    //sorted by length in decreasing order
    vector<MosaicStructure> SubMosaics() const {
        vector<MosaicStructure> answer;
        for (size_t d = block_size(); d > 0; --d) {
            for (size_t i = 0; i + d < block_size(); ++i) {
                answer.push_back(SubMosaic(i, i + d));
            }
        }
        return answer;
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

    bool IsContainedIn(const MosaicStructure& that) const {
        return Find(blocks_, that.blocks_) != -1u;
    }

    MosaicStructure conjugate(const BlockInfoProvider& block_info) const {
        vector<Block> answer;
        for (auto it = blocks_.rbegin(); it != blocks_.rend(); ++it) {
            answer.push_back(block_info.conjugate(*it));
        }
        return MosaicStructure(answer);
    }

};

ostream& operator << (ostream& out, const MosaicStructure& structure) {
    return out << "Mosaic. size=" << structure.block_size() << " blocks=" << structure.blocks();
}

class IntervalIndex {
    multimap<string, StrandRange> all_substruct_pos_;
    const GenomeBlockComposition& block_composition_;

    bool Reported(StrandRange range, const set<StrandRange>& reported) const {
        for (auto r : reported)
            if (r.second == range.second && r.first.contains(range.first))
                return true;
        return false;
    }

public:
    IntervalIndex(const GenomeBlockComposition& block_composition)
    : block_composition_(block_composition) {

    }

    void IndexSubIntervals(const MosaicInterval& interval) {
        VERIFY(interval.support_blocks.front() == interval.pos_range.start_pos);
        VERIFY(interval.support_blocks.back() + 1 == interval.pos_range.end_pos);
        vector<Pos> support = interval.support_blocks;
        for (size_t i = 0; i < support.size(); ++i) {
            for (size_t j = i; j < support.size(); ++j) {
                MosaicInterval sub_interval = interval.SubInterval(Range(support[i], support[j] + 1));
                MosaicStructure sub_mosaic(block_composition_, sub_interval);
                all_substruct_pos_.insert(make_pair(sub_mosaic.Fingerprint(), sub_interval.strand_range()));
                all_substruct_pos_.insert(make_pair(sub_mosaic.conjugate(block_composition_.block_info()).Fingerprint(),
                                                    block_composition_.ConjStrandRange(sub_interval.strand_range())));
            }
        }
    }

    vector<StrandRange> Occurences(const MosaicStructure& mosaic) const {
        return get_all(all_substruct_pos_, mosaic.Fingerprint());
    }

    size_t Mult(const MosaicStructure& mosaic) const {
        return all_substruct_pos_.count(mosaic.Fingerprint());
    }

    vector<StrandRange> UnReportedOccurences(const MosaicStructure& mosaic, const set<StrandRange>& reported) const {
        vector<StrandRange> answer;
        for (StrandRange range : Occurences(mosaic))
            if (!Reported(range, reported) && !Reported(block_composition_.ConjStrandRange(range), reported))
                answer.push_back(range);
        return answer;
    }

    const multimap<string, StrandRange>& all_substruct_pos() const {
        return all_substruct_pos_;
    }
};

class MosaicHelper {
    const BlockInfoProvider& block_info_;
    const GenomeBlockComposition& block_composition_;
    const IntervalIndex& interval_index_;

    vector<Pos> FindSupportBlocks(Range r, bool strand, const MosaicStructure& mosaic) const {
        vector<Pos> answer;
        size_t i = r.start_pos;
        for (Block b : mosaic.blocks()) {
            while (i < r.end_pos && block_composition_.block(i, strand) != b) {
                ++i;
            }
            VERIFY(i < r.end_pos);
            answer.push_back(i);
            ++i;
        }
        VERIFY(answer.front() == r.start_pos && answer.back() == r.end_pos - 1);
        return answer;
    }

public:
    MosaicHelper(const GenomeBlockComposition& block_composition,
                 const IntervalIndex& interval_index)
    : block_info_(block_composition.block_info()),
      block_composition_(block_composition),
      interval_index_(interval_index) {
    }

    size_t Mult(const MosaicStructure& mosaic) const {
        return interval_index_.Mult(mosaic);
    }

    double AvgSpan(const MosaicStructure& mosaic) const {
        double avg = 0.;
        for (StrandRange r : interval_index_.Occurences(mosaic)) {
            avg += double(r.first.size());
        }
        return avg / double(Mult(mosaic));
    }

    double AvgGenomicSpan(const MosaicStructure& mosaic) const {
        double avg = 0.;
        for (StrandRange r : interval_index_.Occurences(mosaic)) {
            Range genomic_range = block_composition_.genome_coords(r);
            avg += double(genomic_range.size());
        }
        return avg / double(interval_index_.Mult(mosaic));
    }

    vector<double> AvgGenomicInterLengths(const MosaicStructure& mosaic) const {
        vector<double> answer(mosaic.block_size() - 1, 0.);
        for (StrandRange r : interval_index_.Occurences(mosaic)) {
            bool strand = r.second;
            vector<Pos> support = FindSupportBlocks(r.first, strand, mosaic);
            for (size_t i = 1; i < support.size(); ++i) {
                answer[i - 1] += double(block_composition_.genome_coords(Range(
                                          support[i - 1] + 1,
                                          support[i]), strand).size());
            }
        }
        for (size_t i = 0; i < answer.size(); ++i) {
            answer[i] /= double(interval_index_.Mult(mosaic));
        }
        return answer;
    }

    size_t length(Block b) const {
        return block_info_.length(b);
    }

    size_t TotalBlockLength(const MosaicStructure& mosaic) const {
        size_t block_length = 0;
        for (Block b : mosaic.blocks()) {
            block_length += length(b);
        }
        return block_length;
    }

    MosaicStructure GetStructure(const MosaicInterval& interval) const {
        return MosaicStructure(block_composition_, interval);
    }

    const BlockInfoProvider& block_info() const {
        return block_info_;
    }

    const GenomeBlockComposition& genome_composition() const {
        return block_composition_;
    }

};

class MosaicPrinter {
public:
    virtual void StartRecord(const MosaicStructure& /*mosaic*/) {

    }

    virtual void ReportSubMosaic(const MosaicStructure& /*mosaic*/, const vector<StrandRange>& /*ranges*/) {

    }

    virtual void EndRecord() {

    }

    virtual ~MosaicPrinter() {}
};

class TxtFileMosaicPrinter : public MosaicPrinter {
    size_t cnt_;
    const MosaicHelper& helper_;
//    const multimap<string, size_t>& different_irred_presence_;
    ofstream out_;

    void BlockInfo(Block b) {
        out_ << b << " (" << helper_.length(b) << ")";
    }

public:

    TxtFileMosaicPrinter(const MosaicHelper& helper,
//                         const multimap<string, size_t>& different_irred_presence,
                         const string& filename) :
                             cnt_(0),
                             helper_(helper),
//                             different_irred_presence_(different_irred_presence),
                             out_(filename) {

    }

    virtual void StartRecord(const MosaicStructure& mosaic) {
        out_ << "Irreducible Mosaic " << cnt_++ << endl;
        out_ << "Support block cnt = " << mosaic.block_size();
        out_ << "; Total support block length = " << helper_.TotalBlockLength(mosaic);
        out_ << "; Full mosaic multiplicity = " << helper_.Mult(mosaic);
        out_ << "; Avg block span = " << helper_.AvgSpan(mosaic);
        out_ << "; Avg genome span = " << helper_.AvgGenomicSpan(mosaic);
        out_ << endl;
        if (helper_.Mult(mosaic) > 1) {
            out_ << "WARN! Full mosaic multiplicity = " << helper_.Mult(mosaic);
            out_ << endl;
        }
        out_ << "Structure: ";
        vector<double> inter_lengths = helper_.AvgGenomicInterLengths(mosaic);
        BlockInfo(mosaic.blocks().front());
        for (size_t i = 1; i < mosaic.block_size(); ++i) {
            out_ << " |...";
            out_ << inter_lengths[i - 1];
            out_ << "...| ";
            BlockInfo(mosaic.blocks()[i]);
        }
        out_ << endl;
        out_ << "Sub_mosaics" << endl;
    }

    virtual void EndRecord() {
        out_ << "............................" << endl;
    }

    virtual void ReportSubMosaic(const MosaicStructure& mosaic, const vector<StrandRange>& ranges) {
        string finger = mosaic.Fingerprint();
        out_ << "------" << endl;
        out_ << "Sub_mosaic. Block cnt = " << mosaic.block_size() << endl;
        out_ << "Blocks " << finger;
//        set<size_t> different_irred;
//        insert_all(different_irred, get_all(different_irred_presence_, finger));
//        out_ << " ; Found in " << different_irred.size() << " different irreducible mosaics";
//        string delim = " (";
//        for (size_t idx : different_irred) {
//            out_ << delim;
//            out_ << idx;
//            delim = ", ";
//        }
//        out_ << ")" << endl;

        string delim = "";
        out_ << "Ranges: ";
        for (StrandRange r : ranges) {
            out_ << delim;
            out_ << "strand: " << (r.second ? "+" : "-") << " ";
            out_ << helper_.genome_composition().genome_coords(r);
            out_ << " (Pos: ";
            out_ << (r.second ? r.first : helper_.genome_composition().ConjStrandRange(r).first);
            out_ << ")";
            delim = "; ";
        }
        out_ << endl;
    }

};

class ParsableFormatPrinter : public MosaicPrinter {
    size_t cnt_;
    const MosaicHelper& helper_;
    vector<pair<MosaicStructure, vector<StrandRange>>> submosaics_;
    ofstream out_;

    void BlockInfo(Block b) {
        out_ << b << " " << helper_.length(b) << endl;
    }

public:

    ParsableFormatPrinter(const MosaicHelper& helper,
                         const string& filename) :
                             cnt_(0),
                             helper_(helper),
                             out_(filename) {

    }

    virtual void StartRecord(const MosaicStructure& mosaic) {
        out_ << cnt_++ << endl; // (the index of the irreducible mosaic)
        out_ << mosaic.block_size() << endl; // (number of the the support blocks)
        out_ << helper_.TotalBlockLength(mosaic) << endl; // (total length)
//        out_ << "; Full mosaic multiplicity = " << helper_.Mult(mosaic);
//        out_ << "; Avg block span = " << helper_.AvgSpan(mosaic);
//        out_ << "; Avg genome span = " << helper_.AvgGenomicSpan(mosaic);
//        out_ << endl;
//        if (helper_.Mult(mosaic) > 1) {
//            out_ << "WARN! Full mosaic multiplicity = " << helper_.Mult(mosaic);
//            out_ << endl;
//        }
//        out_ << "Structure: ";
        vector<double> inter_lengths = helper_.AvgGenomicInterLengths(mosaic);
        BlockInfo(mosaic.blocks().front());
        for (size_t i = 1; i < mosaic.block_size(); ++i) {
            out_ << inter_lengths[i - 1] << endl;
            BlockInfo(mosaic.blocks()[i]);
        }
        submosaics_.clear();
    }

    virtual void ReportSubMosaic(const MosaicStructure& mosaic, const vector<StrandRange>& ranges) {
        submosaics_.push_back(make_pair(mosaic, ranges));
    }

    virtual void EndRecord() {
        out_ << submosaics_.size() << endl;
        size_t cnt = 1;
        for (auto pair : submosaics_) {
            auto mosaic = pair.first;
            auto ranges = pair.second;
            //1 399 733 735 + 1630584 1634815// (the first sub_mosaic structure---1, the blocks---399 733 735, genomic position info--- + 1630584 1634815)
//            string finger = mosaic.Fingerprint();
    //        out_ << "Sub_mosaic. Block cnt = " << mosaic.block_size() << endl;
    //        out_ << "Blocks " << finger;
    //        out_ << " ; Found in " << get_all(different_irred_presence_, finger).size() << " different irreducible mosaics";
    //        string delim = " (";
    //        for (size_t idx : get_all(different_irred_presence_, finger)) {
    //            out_ << delim;
    //            out_ << idx;
    //            delim = ", ";
    //        }
    //        out_ << ")" << endl;

            out_ << cnt++;
            string delim = " ";
    //        out_ << "Ranges: ";
            for (Block b : mosaic.blocks()) {
                out_ << delim;
                out_ << b;
            }
            for (StrandRange r : ranges) {
                out_ << delim;
                out_ << (r.second ? "+" : "-") << " ";
                out_ << helper_.genome_composition().genome_coords(r);
    //            out_ << " (Pos: ";
    //            out_ << (r.second ? r.first : helper_.genome_composition().ConjStrandRange(r).first);
    //            out_ << ")";
    //            delim = "; ";
            }
            out_ << endl;
        }
    }

};

class NotTandemFilter : public func::Predicate<MosaicStructure> {
    const BlockInfoProvider& block_info_;
public:
    NotTandemFilter(const BlockInfoProvider& block_info) :
        block_info_(block_info) {

    }

    //todo use const references!!!
    bool Check (MosaicStructure mosaic) const {
        bag<Block> block_cnts;
        for (Block b : mosaic.blocks()) {
            block_cnts.put(b);
            block_cnts.put(block_info_.conjugate(b));
        }
        size_t tandem_block_cnt = 0;
        for (Block b : mosaic.blocks()) {
            if (block_cnts.mult(b) > 1)
                tandem_block_cnt += 1;
        }
        return tandem_block_cnt * 2 <= mosaic.block_size();
    }

};

class LengthFilter : public func::Predicate<MosaicStructure> {
    const MosaicHelper& helper_;
    size_t min_span_length_;
public:
    LengthFilter(const MosaicHelper& helper, size_t min_span_length) :
        helper_(helper),
        min_span_length_(min_span_length) {

    }

    //todo use const references!!!
    bool Check (MosaicStructure mosaic) const {
        return size_t(math::round(helper_.AvgGenomicSpan(mosaic))) >= min_span_length_;
    }

};

class FullMosaicTracker : public MosaicPrinter {
//    vector<vector<StrandRange>> full_mosaic_ranges_;
    vector<StrandRange> full_mosaic_ranges_;

    size_t curr_length_;
public:
    FullMosaicTracker() : curr_length_(0) {
    }

    virtual void StartRecord(const MosaicStructure& mosaic) {
        curr_length_ = mosaic.block_size();
    }

    virtual void ReportSubMosaic(const MosaicStructure& mosaic, const vector<StrandRange>& ranges) {
        if (mosaic.block_size() == curr_length_)
            full_mosaic_ranges_.push_back(ranges.front());
    }

    const vector<StrandRange>& full_mosaic_ranges() const {
        return full_mosaic_ranges_;
    }

};

class AllRangesTracker : public MosaicPrinter {
    vector<StrandRange> all_ranges_;

public:
    AllRangesTracker() {
    }

    virtual void ReportSubMosaic(const MosaicStructure& /*mosaic*/, const vector<StrandRange>& ranges) {
        all_ranges_.push_back(ranges.front());
    }

    const vector<StrandRange>& all_ranges() const {
        return all_ranges_;
    }

};

class MosaicStructureSet {
    const GenomeBlockComposition& block_composition_;
    IntervalIndex& interval_index_;

    vector<MosaicInterval> raw_intervals_;
    vector<MosaicStructure> irreducible_structures_;
//    multimap<string, size_t> different_irred_presence_;

    shared_ptr<func::Predicate<MosaicStructure>> filter_;
    shared_ptr<func::Predicate<MosaicStructure>> sub_filter_;


    MosaicStructure ConjugateMosaic(const MosaicStructure mosaic) const {
        return mosaic.conjugate(block_composition_.block_info());
    }

//    void CountDifferentIrred(const MosaicStructure& mosaic, size_t idx) {
//        for (size_t i = 0; i < mosaic.block_size(); ++i) {
//            for (size_t j = i; j < mosaic.block_size(); ++j) {
//                MosaicStructure sub_mosaic = mosaic.SubMosaic(i, j);
//                different_irred_presence_.insert(make_pair(sub_mosaic.Fingerprint(), idx));
//                different_irred_presence_.insert(
//                        make_pair(ConjugateMosaic(sub_mosaic).Fingerprint(), idx));
//            }
//        }
//    }
//
//    void CountDifferentIrred() {
//        for (size_t i = 0; i < irreducible_structures_.size(); ++i) {
//            CountDifferentIrred(irreducible_structures_[i], i);
//        }
//    }

    bool AnalyzeStructure(const MosaicStructure& mosaic) {
        for (auto& irred_struct : irreducible_structures_) {
            if (mosaic.IsContainedIn(irred_struct)) {
                DEBUG("Contained in some irred structure");
                return false;
            }
            if (mosaic.IsContainedIn(ConjugateMosaic(irred_struct))) {
                DEBUG("Contained in conjugate of some irred structure");
                return false;
            }
        }
        irreducible_structures_.push_back(mosaic);
        return true;
    }

public:
    MosaicStructureSet(const GenomeBlockComposition& block_composition,
                       IntervalIndex& interval_index) :
                           block_composition_(block_composition),
                           interval_index_(interval_index) {

    }

    void ProcessInterval(const MosaicInterval& interval) {
        interval_index_.IndexSubIntervals(interval);
        if (interval.support_size() > 1) {
            raw_intervals_.push_back(interval);
        }
    }

//    const multimap<string, size_t>& different_irred_presence() const {
//        return different_irred_presence_;
//    }

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
        //CountDifferentIrred();
    }

    void set_structure_filter(shared_ptr<func::Predicate<MosaicStructure>> filter) {
        filter_ = filter;
    }

    void set_substructure_filter(shared_ptr<func::Predicate<MosaicStructure>> sub_filter) {
        sub_filter_ = sub_filter;
    }

    void Report(MosaicPrinter& printer) {
        for (size_t i = 0; i < irreducible_structures_.size(); ++i) {
            MosaicStructure mosaic = irreducible_structures_[i];
            if (!filter_ || filter_->Check(mosaic)) {
                printer.StartRecord(mosaic);
                set<StrandRange> reported_ranges;
                for (const MosaicStructure& submosaic: mosaic.SubMosaics()) {
                    vector<StrandRange> ranges = interval_index_
                            .UnReportedOccurences(submosaic, reported_ranges);
                    if (ranges.empty())
                        continue;
                    if (!sub_filter_ || sub_filter_->Check(mosaic))
                        printer.ReportSubMosaic(submosaic, ranges);
                    insert_all(reported_ranges, ranges);
                }
                printer.EndRecord();
            }
        }
    }

    void Report(const vector<shared_ptr<MosaicPrinter>>& printers) {
        for (auto printer_ptr : printers) {
            Report(*printer_ptr);
        }
    }
};

//todo reduce duplication
const io::SingleRead MakeRead(const string& read, const string& name = "") {
    //todo fill with good quality
    std::string qual;
    qual.resize(read.size());
    return io::SingleRead(name, read, qual);
}

const vector<io::SingleRead> MakeReads(const vector<string>& reads) {
    vector<io::SingleRead> ans;
    for (size_t i = 0; i < reads.size(); ++i) {
        ans.push_back(MakeRead(reads[i]));
    }
    return ans;
}

const vector<io::SingleRead> MakeReads(const vector<string>& reads, const vector<string>& names) {
    vector<io::SingleRead> ans;
    for (size_t i = 0; i < reads.size(); ++i) {
        ans.push_back(MakeRead(reads[i], names[i]));
    }
    return ans;
}

vector<string> mosaic_names(size_t n) {
    vector<string> ans;
    for (size_t i = 0; i < n; ++i) {
        ans.push_back("mosaic_" + ToString(i));
    }
    return ans;
}

shared_ptr<io::ReadStream<io::SingleRead>> StreamInstance(const vector<string>& sequences) {
    return make_shared<io::VectorReadStream<io::SingleRead>>(MakeReads(sequences));
}

shared_ptr<io::ReadStream<io::SingleRead>> StreamInstance(const vector<string>& sequences, const vector<string>& names) {
    return make_shared<io::VectorReadStream<io::SingleRead>>(MakeReads(sequences, names));
}
//multimap<string, StrandRange> all_substruct_pos_;
//const GenomeBlockComposition& block_composition_;

string ExtractSequence(const StrandRange& strand_range, const GenomeBlockComposition& block_composition) {
    const BlockInfoProvider& block_info = block_composition.block_info();
    const Graph& g = block_info.g();
    vector<Graph::EdgeId> edges;
    Range range = strand_range.first;
    bool strand = strand_range.second;
    for (size_t i = range.start_pos; i < range.end_pos; ++i) {
        edges.push_back(block_info.edge_id(block_composition.block(i, strand)));
    }
    return MergeSequences(g, edges).str();
}

vector<string> ExtractSequences(const vector<StrandRange>& strand_ranges, const GenomeBlockComposition& block_composition) {
    vector<string> answer;
    for (StrandRange range : strand_ranges) {
        answer.push_back(ExtractSequence(range, block_composition));
    }
    return answer;
}

void DrawGraph(const vector<StrandRange>& all_ranges,
               const vector<StrandRange>& full_mosaic_ranges,
               const GenomeBlockComposition& block_composition) {
    make_dir("tmp");
    graph_pack<Graph, RtSeq> gp(block_composition.block_info().g().k(), "tmp", 0);

    auto stream = io::RCWrap(StreamInstance(ExtractSequences(all_ranges, block_composition)));
    auto streams = io::ReadStreamList<io::SingleRead>(stream);
//    ConstructGraphUsingOldIndex(streams, gp.g, gp.index);
    ConstructGraph(config::debruijn_config::construction(), omp_get_max_threads(), destreams, gp.g, gp.index);

    auto full_mosaic_pos_stream = io::RCWrap(StreamInstance(ExtractSequences(full_mosaic_ranges, block_composition), mosaic_names(full_mosaic_ranges.size())));
    INFO("Threading " << full_mosaic_ranges.size() << " full mosaics");
    visualization::position_filler::FillPos(gp, *full_mosaic_pos_stream);

    visualization::graph_labeler::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);

    shared_ptr<GraphSplitter<Graph>> splitter = omnigraph::ReliableSplitter(gp.g,
            numeric_limits<size_t>::max(),
            50
            /*numeric_limits<size_t>::max()*/);

    path::remove_if_exists("mosaic_pics");
    path::make_dir("mosaic_pics");
    INFO("Writing components");
    visualization::visualization_utils::WriteComponents(gp.g, "mosaic_pics/", splitter,
            visualization::graph_colorer::DefaultColorer(gp.g), labeler);
    INFO("Components written");
}

class MosaicStructureAnalyzer {
    const BlockInfoProvider& block_info_;
    const GenomeBlockComposition& block_composition_;

    size_t min_support_block_length_;
    size_t max_support_block_multiplicity_;
    size_t max_inter_block_length_;
    size_t min_reportable_mosaic_length_;
    size_t min_reportable_submosaic_length_;
    string folder_;

    Pos curr_pos_;

    bool CheckSupporting(Block b) {
        size_t mult = block_composition_.multiplicity(b);
        return mult > 1 && block_info_.length(b) > min_support_block_length_
                && mult < max_support_block_multiplicity_;
    }

    bool CheckSupporting() {
        return CheckSupporting(block_composition_.block(curr_pos_, true));
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
                ? block_composition_.genome_coords(curr_pos_, true).start_pos
                - block_composition_.genome_coords(init_pos, true).end_pos
                : -1u;
    }

public:
    MosaicStructureAnalyzer(const GenomeBlockComposition& block_composition,
                            size_t min_support_length, size_t max_support_mult,
                            size_t max_inter_length, size_t min_reportable_mosaic_length,
                            size_t min_reportable_submosaic_length,
                            const string& folder)
             : block_info_(block_composition.block_info()),
               block_composition_(block_composition),
               min_support_block_length_(min_support_length),
               max_support_block_multiplicity_(max_support_mult),
               max_inter_block_length_(max_inter_length),
               min_reportable_mosaic_length_(min_reportable_mosaic_length),
               min_reportable_submosaic_length_(min_reportable_submosaic_length),
               folder_(folder),
               curr_pos_(0) {
    }

    void Analyze() {
        IntervalIndex interval_index(block_composition_);
        MosaicHelper helper(block_composition_,
                            interval_index);

        MosaicStructureSet interval_set(block_composition_, interval_index);
        INFO("Collecting mosaic intervals");
        MoveToNextSupportBlock();
        while (curr_pos_ < block_composition_.size()) {
            MosaicInterval interval(true, curr_pos_);
            while (MoveToNextSupportBlock() < max_inter_block_length_) {
                interval.pos_range.end_pos = (curr_pos_ + 1);
                interval.support_blocks.push_back(curr_pos_);
            }
            interval_set.ProcessInterval(interval);
        }
        INFO("Analyzing intervals and forming mosaic structures");
        interval_set.Analysis();
        INFO("Reporting mosaic structures");

        //might have problems if only largest occurence is tandem making whole mosaic invalid
        //todo magic constant!
        auto filter = func::And<MosaicStructure>(make_shared<NotTandemFilter>(block_info_),
                                                 make_shared<LengthFilter>(helper, min_reportable_mosaic_length_));
//        auto filter = make_shared<func::AlwaysTrue<MosaicStructure>>();

        interval_set.set_structure_filter(filter);
        //todo move filter to set constructor

        //pics
        FullMosaicTracker tracker;
        interval_set.Report(tracker);
        vector<StrandRange> full_mosaic_ranges = tracker.full_mosaic_ranges();
        AllRangesTracker all_tracker;
        interval_set.Report(all_tracker);
        DrawGraph(all_tracker.all_ranges(), full_mosaic_ranges, block_composition_);
        //end pics

        interval_set.set_substructure_filter(make_shared<LengthFilter>(helper, min_reportable_submosaic_length_));

        ParsableFormatPrinter parsable_printer(helper, folder_ + "mosaic_to_parse.txt");
        interval_set.Report(parsable_printer);

        TxtFileMosaicPrinter readable_printer(helper, /*interval_set.different_irred_presence(),*/
                                              folder_ + "mosaic_to_read.txt");
        interval_set.Report(readable_printer);

    }

};

template<class gp_t>
void PerformMosaicAnalysis(const gp_t& gp, const MappingPath<EdgeId> mapping_path, const Sequence& genome,
                           size_t min_support_length, size_t max_support_mult,
                           size_t max_inter_length, size_t min_reportable_mosaic_length,
                           size_t min_reportable_submosaic_length, const string& folder) {
    BlockInfoProvider block_info(gp.g);
    GenomeBlockComposition block_composition(gp.g, mapping_path, genome.size() - gp.g.k(), block_info);

    MosaicStructureAnalyzer analyzer(block_composition, min_support_length, max_support_mult,
                                     max_inter_length, min_reportable_mosaic_length,
                                     min_reportable_submosaic_length, folder);
    analyzer.Analyze();
}

}
}
