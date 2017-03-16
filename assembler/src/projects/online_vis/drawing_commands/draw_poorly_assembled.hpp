//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "io/reads/wrapper_collection.hpp"
#include <boost/algorithm/string.hpp>
#include "assembly_graph/core/basic_graph_stats.hpp"

#include <boost/algorithm/string/predicate.hpp>

namespace online_visualization {

//class RepeatProcessor {
//    const Graph& g_;
//
//protected:
//    const Graph& g() const {
//        return g_;
//    }
//
//public:
//    virtual void ProcessResolved(EdgeId e1, EdgeId e2, size_t gap_length, const vector<EdgeId> repeat_edges) {
//    }
//
//    virtual void ProcessUnresolved(EdgeId e1, EdgeId e2, size_t gap_length, const vector<EdgeId> repeat_edges) {
//    }
//
//    virtual ~RepeatProcessor() {
//    }
//};
//
//class CompositeRepeatProcessor : public RepeatProcessor {
//    vector<shared_ptr<RepeatProcessor>> processors_;
//
//public:
//    virtual void ProcessResolved(EdgeId e1, EdgeId e2, size_t gap_length, const vector<EdgeId> repeat_edges) {
//        for (auto p : processors_) {
//            p->ProcessResolved(e1, e2, gap_length, repeat_edges);
//        }
//    }
//
//    virtual void ProcessUnresolved(EdgeId e1, EdgeId e2, size_t gap_length, const vector<EdgeId> repeat_edges) {
//        for (auto p : processors_) {
//            p->ProcessUnresolved(e1, e2, gap_length, repeat_edges);
//        }
//    }
//
//};

struct RepeatInfo {
    EdgeId e1;
    EdgeId e2;
    size_t genomic_gap;
    vector<EdgeId> ref_path;
    string seq_name;
    //number of repeat in seq
    size_t local_cnt;

    RepeatInfo(EdgeId e1_,
    EdgeId e2_,
    size_t genomic_gap_,
    const vector<EdgeId>& ref_path_,
    string seq_name_,
    size_t local_cnt_) :
        e1(e1_), e2(e2_),
        genomic_gap(genomic_gap_), ref_path(ref_path_),
        seq_name(seq_name_), local_cnt(local_cnt_) {

    }
};

class RepeatProcessor {
public:
    virtual ~RepeatProcessor() {
    }

    virtual void ProcessUnresolved(DebruijnEnvironment&, const RepeatInfo&) const {}

    virtual void ProcessResolved(DebruijnEnvironment&, const RepeatInfo&) const {}
};

class StructuredFileLogger : public RepeatProcessor {

    string Log(const GraphPack& gp, const RepeatInfo& repeat_info) const {
        return fmt::format("{:d} {:d} {:d} {:d} {:d}", repeat_info.genomic_gap, repeat_info.ref_path.size(),
                                CumulativeLength(gp.g, repeat_info.ref_path),
                                gp.g.int_id(repeat_info.e1), gp.g.int_id(repeat_info.e2));
    }

public:
    virtual void ProcessUnresolved(DebruijnEnvironment& curr_env, const RepeatInfo& repeat_info) const {
        cerr << "RI: 0 " << Log(curr_env.graph_pack(), repeat_info) << endl;

    }

    virtual void ProcessResolved(DebruijnEnvironment& curr_env, const RepeatInfo& repeat_info) const {
        cerr << "RI: 1 " << Log(curr_env.graph_pack(), repeat_info) << endl;
    }
};

class ReadableUnresolvedLogger : public RepeatProcessor {

public:
    virtual void ProcessUnresolved(DebruijnEnvironment& curr_env, const RepeatInfo& repeat_info) const {
        LOG(fmt::format("Genomic gap: {:d}, number of edges: {:d}, edge 1: {:s}, edge 2 {:s}", repeat_info.genomic_gap,
                        repeat_info.ref_path.size(), curr_env.graph().str(repeat_info.e1),
                        curr_env.graph().str(repeat_info.e2)));
    }

};

class UnresolvedPrinter : public RepeatProcessor {

    void DrawGap(DebruijnEnvironment& curr_env, const vector<EdgeId>& path, string filename, string /*label*/ = "") const {
        visualization::visualization_utils::WriteComponentsAlongPath<Graph>(curr_env.graph(), path, filename, curr_env.coloring(), curr_env.labeler());
        LOG("The pictures is written to " << filename);
    }

public:

    virtual void ProcessUnresolved(DebruijnEnvironment& curr_env, const RepeatInfo& repeat_info) const {
        make_dir(curr_env.folder());
        string pics_folder = curr_env.folder() + "/" + curr_env.GetFormattedPictureCounter()  + "_" + repeat_info.seq_name + "/";
        make_dir(pics_folder);
        string pic_name = std::to_string(repeat_info.local_cnt) + "_" +  std::to_string(repeat_info.genomic_gap) +
                "_" + std::to_string(curr_env.graph().int_id(repeat_info.e1)) + "_" + std::to_string(curr_env.graph().int_id(repeat_info.e2)) + "_";

        DrawGap(curr_env, repeat_info.ref_path, pics_folder + pic_name);
    }

};

class PairedInfoChecker : public RepeatProcessor {

    bool CheckInfo(const omnigraph::de::PairedInfoIndexT<Graph>& clustered_pi_idx, EdgeId e1, EdgeId e2) const {
        //return !clustered_pi_idx.Get(e1, e2).empty();
        //We don't store empty histograms, do we?
        return clustered_pi_idx.contains(e1, e2);
    }

public:

    virtual void ProcessUnresolved(DebruijnEnvironment& curr_env, const RepeatInfo& repeat_info) const {
        const omnigraph::de::PairedInfoIndexT<Graph>& clustered_pi_idx = curr_env.graph_pack().clustered_indices[0];
        const Graph& g = curr_env.graph();
        vector<EdgeId> edges;
        edges.push_back(repeat_info.e1);
        utils::push_back_all(edges, repeat_info.ref_path);
        edges.push_back(repeat_info.e2);
        for (EdgeId e : edges) {
            if (!CheckInfo(clustered_pi_idx, repeat_info.e1, e)) {
                cerr << "NO_PI: " << g.int_id(repeat_info.e1) <<
                        " " << g.int_id(repeat_info.e2) <<
                        " " << g.int_id(e) << endl;
            }
        }
    }

};

class DrawUnresolvedRepeatsCommand : public DrawingCommand {
private:
    const string ref_prefix_;
    const size_t gap_diff_threshold_;
    const double good_mapping_coeff_;
    vector<shared_ptr<RepeatProcessor>> processors_;

    vector<EdgePosition> GatherPositions(const GraphPack& gp, EdgeId e, const string& prefix) const {
        vector<EdgePosition> answer;
        for (EdgePosition pos : gp.edge_pos.GetEdgePositions(e)) {
            if (boost::starts_with(pos.contigId, prefix)) {
                answer.push_back(pos);
            }
        }
        return answer;
    }

//    bool IsNext(MappingRange mr1, MappingRange mr2, size_t max_gap = 5) const {
//       return mr2.initial_range.start_pos >= mr1.initial_range.end_pos && mr2.initial_range.start_pos <= mr1.initial_range.end_pos + max_gap;
//    }
//
//    shared_ptr<MappingRange> FindNextRange(MappingRange curr_range, const set<MappingRange>& ranges) const {
//        cout << "Looking for next range for " << curr_range << endl;
//        for (auto r : ranges) {
//            cout << "Considering range " << r << endl;
//            if (IsNext(curr_range, r)) {
//                cout << "Found next" << endl;
//                return make_shared<MappingRange>(r);
//            }
//        }
//        cout << "Couldn't find suitable range" << endl;
//        return shared_ptr<MappingRange>(0);
//    }
//    
//    shared_ptr<pair<EdgeId, EdgePosition>> NextEdge(const GraphPack& gp, VertexId v, EdgePosition curr_pos) const {
//        for (EdgeId next_e : gp.g.OutgoingEdges(v)) {
//            cout << "Considering " << gp.g.str(next_e) << " as next edge " << endl;
//            set<MappingRange> relevant_ranges = gp.edge_pos.GetEdgePositions(next_e, curr_pos.contigId);
//            auto next_range = FindNextRange(curr_pos.mr, relevant_ranges);
//            cout << "Considered " << relevant_ranges.size() << " relevant ranges" << endl;
//            if (next_range) {
//                cout << "Found next edge" << endl;
//                return make_shared<pair<EdgeId, EdgePosition>>(next_e, EdgePosition(curr_pos.contigId, *next_range));
//            }
//        }
//        cout << "Couldn't find next edge" << endl;
//        return shared_ptr<pair<EdgeId, EdgePosition>>(0);
//    }
//
//    vector<EdgeId> FindReferencePath(const GraphPack& gp, EdgeId e1, EdgeId e2) const {
//        EdgePosition curr_pos = GatherPositions(gp, e1, ref_prefix_).front();
//        VertexId curr_v = gp.g.EdgeEnd(e1);
//        vector<EdgeId> answer;
//        answer.push_back(e1);
//    //    for (size_t i = 0 ; i < 1000 ; ++i) {
//        while(true) {
//            auto next_info = NextEdge(gp, curr_v, curr_pos);
//            if (next_info) {
//                EdgeId next_e = next_info->first;
//                answer.push_back(next_e);
//                if (next_e == e2) {
//                    break;
//                }
//                curr_v = gp.g.EdgeEnd(next_e);
//                curr_pos = next_info->second;
//            } else {
//                return vector<EdgeId>();
//            }
//        }
//        return answer;
//    }

    MappingPath<EdgeId> FindReferencePath(const GraphPack& gp, EdgeId e1, EdgeId e2) const {
        auto e1_poss = GatherPositions(gp, e1, ref_prefix_);
        auto e2_poss = GatherPositions(gp, e2, ref_prefix_);
        VERIFY(e1_poss.size() == 1 && e2_poss.size() == 1);
        EdgePosition e1_pos = e1_poss.front();
        EdgePosition e2_pos = e2_poss.front();
        VERIFY(e1_pos.contigId == e1_pos.contigId);
        Sequence ref = (e1_pos.contigId == "ref0") ? gp.genome.GetSequence() : !gp.genome.GetSequence();
        size_t gap_start = e1_pos.mr.initial_range.end_pos;
        size_t gap_end = e2_pos.mr.initial_range.start_pos + gp.g.k();
        VERIFY(gap_end >= gap_start && gap_end <= ref.size());
        Sequence gap_fragment = ref.Subseq(gap_start, gap_end);
        return debruijn_graph::MapperInstance(gp)->MapSequence(gap_fragment);
    }

    size_t GenomicGap(const GraphPack& gp, EdgeId e1, EdgeId e2) const {
        auto poss1 = GatherPositions(gp, e1, ref_prefix_);
        auto poss2 = GatherPositions(gp, e2, ref_prefix_);
        VERIFY_MSG(poss1.size() == 1, "Positions of first edge " << poss1);
        VERIFY_MSG(poss2.size() == 1, "Positions of second edge " << poss2);
        if (poss1.front().contigId != poss2.front().contigId) {
            WARN("Assembly contig stitching edges from different strains");
            return -1u;
        }

        MappingRange r1 = poss1.front().mr;
        MappingRange r2 = poss2.front().mr;
        if (r2.initial_range.start_pos < r1.initial_range.end_pos) {
            WARN("Wrong order of edges in the contig");
            return -1u;
        }

        return r2.initial_range.start_pos - r1.initial_range.end_pos;
    }

    set<string> GatherNames(const GraphPack& gp, EdgeId e, const string& prefix) const {
        set<string> answer;
        for (auto pos : GatherPositions(gp, e, prefix)) {
            answer.insert(pos.contigId);
        }
        return answer;
    }

    bool IsOfMultiplicityOne(const GraphPack& gp, EdgeId e) const {
        return GatherPositions(gp, e, ref_prefix_).size() == 1;
    }
    
    bool BelongToSameContig(const GraphPack& gp, EdgeId e1, EdgeId e2, string assembly_prefix) const {
        auto names1 = GatherNames(gp, e1, assembly_prefix);
        auto names2 = GatherNames(gp, e2, assembly_prefix);
        //checking non-empty intersection
        for (auto name: names1) {
            if (names2.count(name)) {
                return true;
            }
        }
        return false;
    }

    bool CheckMapping(size_t edge_length, size_t mapped_range_length) const {
        return double(mapped_range_length) > good_mapping_coeff_ * double(edge_length);
    }

    vector<EdgeId> EdgesOfInterest(const GraphPack& gp, const MappingPath<EdgeId>& mapping_path, size_t length_threshold) const {
        vector<EdgeId> answer;
        for (size_t i = 0; i < mapping_path.size(); ++i) {
            EdgeId e = mapping_path[i].first;
            if (gp.g.length(e) >= length_threshold 
                    && IsOfMultiplicityOne(gp, e)
                    && CheckMapping(gp.g.length(e), mapping_path[i].second.mapped_range.size())) {
                answer.push_back(e);
            }
        }
        return answer;
    }

protected:

    bool AnalyzeGaps(DebruijnEnvironment& curr_env,
                                        io::SingleRead contig,
                                        string base_assembly_prefix, 
                                        size_t edge_length,
                                        size_t max_genomic_gap, 
                                        size_t max_gap_cnt = -1u) const {
        GraphPack& gp = curr_env.graph_pack();
        auto mapper_ptr = debruijn_graph::MapperInstance(gp);
        MappingPath<EdgeId> mapping_path = mapper_ptr->MapRead(contig);
        auto pos_handler = gp.edge_pos;
        
        auto long_unique = EdgesOfInterest(gp, mapping_path, edge_length);

        bool found_smth = false;
        size_t cnt = 0;
        for (size_t i = 1; i < long_unique.size(); ++i) {
            if (max_gap_cnt != -1u && cnt >= max_gap_cnt) {
                INFO("Number of gaps exceeded " << max_gap_cnt);
                return found_smth;
            }
            EdgeId e1 = long_unique[i-1];
            EdgeId e2 = long_unique[i];

            size_t contig_gap = mapping_path[i].second.initial_range.start_pos - mapping_path[i-1].second.initial_range.end_pos;
            size_t genomic_gap = GenomicGap(gp, e1, e2);
            if (genomic_gap == -1u) {
                DEBUG("Contig likely misassembled. Unique long regions in wrong order. e1 " << 
                        gp.g.str(e1) << " genome pos : " << GatherPositions(gp, e1, "ref") << " and e2 " << gp.g.str(e2) << 
                        " genome pos : " << GatherPositions(gp, e2, "ref"));
                continue;
            }
            DEBUG("Found genomic gap " << genomic_gap << 
                " between e1 " << gp.g.str(e1) << " genome pos : " << GatherPositions(gp, e1, "ref") << " and e2 " << gp.g.str(e2)
                << " genome pos : " << GatherPositions(gp, e2, "ref"));

            DEBUG("Looking for reference path");
            auto ref_mapping_path = FindReferencePath(gp, e1, e2);
            if (ref_mapping_path.size() == 0) {
                DEBUG("Couldn't find ref path between " << gp.g.str(e1) << " and " << gp.g.str(e2));
                continue;
            } 

            vector<EdgeId> ref_path = curr_env.path_finder().FindReadPath(ref_mapping_path);
            if (ref_path.empty()) {
                DEBUG("Couldn't fix ref path");
                ref_path = ref_mapping_path.simple_path();
            }
            DEBUG("Found ref path between " << gp.g.str(e1) << " and " << gp.g.str(e2));
            DEBUG(ref_path.size() << " edges of cumulative length " << CumulativeLength(gp.g, ref_path));
                    
            if (std::abs(int(genomic_gap) - int(contig_gap)) >= int(gap_diff_threshold_)) {
                DEBUG("Contig likely misassembled. Genomic gap is " << genomic_gap << " while contig gap was " << contig_gap);
                continue;
            }

            if (genomic_gap >= max_genomic_gap) {
                DEBUG("Genomic gap exceeded max_gap value and will be skipped. Gap " << genomic_gap << " max_gap " << max_genomic_gap);
                continue;
            }
            
            RepeatInfo info(e1, e2, genomic_gap, ref_path, contig.name(), cnt++);
            if (!BelongToSameContig(gp, e1, e2, base_assembly_prefix)) {
                DEBUG("Long unique edges not in the same contig of base assembly");
                for (auto processor : processors_) {
                    processor->ProcessUnresolved(curr_env, info);
                }
                found_smth = true;
            } else {
                for (auto processor : processors_) {
                    processor->ProcessResolved(curr_env, info);
                }
            }
        }
        return found_smth;
    }

    DrawUnresolvedRepeatsCommand(const string& command_name,
                                 const vector<shared_ptr<RepeatProcessor>>& processors)
            : DrawingCommand(command_name), ref_prefix_("ref"), gap_diff_threshold_(1000),
              good_mapping_coeff_(0.7), processors_(processors) {
    }

};

class DrawUnresolvedWRTAssemblyCommand : public DrawUnresolvedRepeatsCommand {

protected:
    size_t MinArgNumber() const {
        return 3;
    }

    bool CheckCorrectness(const vector<string>& args) const {
        if (!CheckEnoughArguments(args))
            return false;

        return true;
    }

public:
    DrawUnresolvedWRTAssemblyCommand() : 
        DrawUnresolvedRepeatsCommand("draw_unresolved_wrt_assembly",
                                     vector<shared_ptr<RepeatProcessor>>{make_shared<StructuredFileLogger>()}) {
    }

    string Usage() const {
        string answer;
        answer = answer + "Command `draw_unresolved_wrt_assembly` \n" + "Usage:\n"
                + "> draw_unresolved_wrt_assembly <contigs_file> <prefix_of_base_assembly> <unique_edge_length> [first N contigs to analyze]\n"
                + " Draws pictures of unresolved repeats.";
        return answer;
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        std::string contigs_file = args[1];
        string base_assembly_prefix = args[2];
        size_t edge_length = std::stoll(args[3]);

        if (!CheckFileExists(contigs_file)) {
            LOG("File with contigs " << contigs_file << " not found");
            return;
        }

        size_t contig_cnt = -1u;
        if (args.size() > 4) {
            LOG("Will analyze first " << args[4] << " contigs");
            contig_cnt = std::stoll(args[4]);
        }
        
        io::FileReadStream reader(contigs_file);

        size_t i = 0;
        while (!reader.eof() && i < contig_cnt) {
            io::SingleRead contig;
            reader >> contig;
            LOG("Considering contig " << contig.name());

            if (AnalyzeGaps(curr_env, contig, base_assembly_prefix,
                            edge_length, numeric_limits<size_t>::max())) {
                curr_env.inc_pic_counter();
            }
            ++i;
        }

    }

};

class DrawUnresolvedWRTReferenceCommand : public DrawUnresolvedRepeatsCommand {

protected:
    size_t MinArgNumber() const {
        return 3;
    }

    bool CheckCorrectness(const vector<string>& args) const {
        if (!CheckEnoughArguments(args))
            return false;

        return true;
    }

public:
    DrawUnresolvedWRTReferenceCommand() : 
        DrawUnresolvedRepeatsCommand("draw_unresolved_wrt_reference",
                                     vector<shared_ptr<RepeatProcessor>>{make_shared<StructuredFileLogger>(),
                                     make_shared<PairedInfoChecker>()}) {
    }

    string Usage() const {
        string answer;
        answer = answer + "Command `draw_unresolved_wrt_reference ` \n" + "Usage:\n"
                + "> draw_unresolved_wrt_reference <gap_length> <prefix_of_base_assembly> <unique_edge_length> [first N gaps to analyze]\n"
                + " Draws pictures of unresolved repeats longer then gap_length between unique edges longer than some constant.";
        return answer;
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        size_t max_interesting_gap = std::stoll(args[1]);
        std::string base_assembly_prefix = args[2];
        size_t edge_length = std::stoll(args[3]);

        if (curr_env.graph_pack().genome.size() == 0) {
            LOG("Reference genome hasn't been loaded");
            return;
        }

        size_t gap_cnt = -1u;
        if (args.size() > 4) {
            LOG("Will analyze first " << args[4] << " gaps");
            gap_cnt = std::stoll(args[4]);
        }
        
        io::SingleRead ref_as_read("ref", curr_env.graph_pack().genome.str());
        AnalyzeGaps(curr_env, ref_as_read, base_assembly_prefix,
                    edge_length, max_interesting_gap, gap_cnt);
    }

};

class DrawPoorlyAssembledCommand : public DrawingCommand {
    const double WELL_ASSEMBLED_CONSTANT = 0.7;
private:
    
    bool IsPoorlyAssembled(const GraphPack& gp, io::SingleRead contig, string base_assembly_prefix) const {
        MappingPath<EdgeId> mapping_path = debruijn_graph::MapperInstance(gp)->MapRead(contig);
        auto pos_handler = gp.edge_pos;
        map<string, size_t> base_ctg_2_len;
        for (EdgeId e : mapping_path.simple_path()) {
            auto positions = pos_handler.GetEdgePositions(e);
            for (EdgePosition pos : positions) {
                if (boost::starts_with(pos.contigId, base_assembly_prefix)) {
                    base_ctg_2_len[pos.contigId] += pos.mr.mapped_range.size();
                }
            }
        }
        for (pair<string, size_t> entry : base_ctg_2_len) {
            if (double(entry.second) > double(contig.size()) * WELL_ASSEMBLED_CONSTANT) {
                LOG("Contig " << contig.name() << 
                        " was well covered by contig " << entry.first << " of base assembly")
                return false;
            }
        }

        return true;
    }

protected:
    size_t MinArgNumber() const {
        return 2;
    }

    bool CheckCorrectness(const vector<string>& args) const {
        if (!CheckEnoughArguments(args))
            return false;

        return true;
    }

public:
    string Usage() const {
        string answer;
        answer = answer + "Command `draw_poorly_assembled` \n" + "Usage:\n"
                + "> draw_poorly_assembled <contigs_file> <prefix_of_base_assembly> [first N contigs to analyze]\n"
                + " Draws pictures of contigs that are not well covered with any contig in base assembly.";
        return answer;
    }

    DrawPoorlyAssembledCommand()
            : DrawingCommand("draw_poorly_assembled") {
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        std::string contigs_file = args[1];
        string base_assembly_prefix = args[2];

        if (!CheckFileExists(contigs_file)) {
            LOG("File with contigs " << contigs_file << " not found");
            return;
        }

        size_t contig_cnt = -1u;
        if (args.size() > 3) {
            LOG("Will analyze first " << args[3] << " contigs");
            contig_cnt = std::stoll(args[3]);
        }
        
        io::FileReadStream reader(contigs_file);

        size_t i = 0;
        while (!reader.eof() && i < contig_cnt) {
            io::SingleRead contig;
            reader >> contig;
            LOG("Considering contig " << contig.name());

            if (IsPoorlyAssembled(curr_env.graph_pack(), contig, base_assembly_prefix)) {
                LOG("Was poorly assembled, drawing");
                DrawPicturesAlongContig(curr_env, contig);
            } else {
                LOG("Was well assembled");
            }

            ++i;
        }

    }

};

class DrawCoverageDropsCommand : public DrawingCommand {
    const size_t cov_drop = 25;
    const size_t min_ende_len = 2000;
private:

    bool IsRepeat(const GraphPack& gp, EdgeId e) const {
        auto v1 = gp.g.EdgeStart(e);
        auto v2 = gp.g.EdgeEnd(e);
        return gp.g.IncomingEdgeCount(v1) >= 2 || gp.g.OutgoingEdgeCount(v2) >= 2 ;
    }

    std::vector<std::vector<EdgeId>> Split(const GraphPack& gp, std::vector<EdgeId> mapping_path) const {
        std::vector<std::vector<EdgeId>> answer;
        std::vector<EdgeId> temp;
        for(auto e : mapping_path) {

            if(gp.g.OutgoingEdgeCount(gp.g.EdgeEnd(e)) == 0) {
                temp.push_back(e);
                answer.push_back(temp);
                temp.clear();
                continue;
            }

            if(gp.g.IncomingEdgeCount(gp.g.EdgeStart(e)) == 0) {
                answer.push_back(temp);
                temp.clear();
                temp.push_back(e);
                continue;
            }
            temp.push_back(e);
        }
        if(temp.size() > 0) {
            answer.push_back(temp);
        }
        return answer;
    }

    bool HasCoverageDrops(const GraphPack& gp, std::vector<EdgeId> mapping_path) const {
        double min_coverage = std::numeric_limits<double>::max();
        double max_coverage = 0;

        std::for_each(mapping_path.begin(), mapping_path.end(), [this, &max_coverage, &min_coverage, &gp](EdgeId e){
            if(!IsRepeat(gp, e) && gp.g.length(e) > min_ende_len) {
                min_coverage = std::min(gp.g.coverage(e), min_coverage);
                max_coverage = std::max(gp.g.coverage(e), max_coverage);
            }
        });
        if(max_coverage > min_coverage && max_coverage - min_coverage > cov_drop) {
            return true;
        }
        return false;
    }

protected:
    size_t MinArgNumber() const {
        return 1;
    }

    bool CheckCorrectness(const vector<string>& args) const {
        if (!CheckEnoughArguments(args))
            return false;
        return true;
    }

public:
    string Usage() const {
        string answer;
        answer = answer + "Command `draw_coverage_drops` \n" + "Usage:\n"
                + "> draw_coverage_drops <contigs_file> [first N contigs to analyze]\n"
                + " Draws pictures of contigs that have substantial coverage drops during theirs alignments to the graph.";
        return answer;
    }

    DrawCoverageDropsCommand()
            : DrawingCommand("draw_coverage_drops") {
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        std::string contigs_file = args[1];
        if (!CheckFileExists(contigs_file)) {
            LOG("File with contigs " << contigs_file << " not found");
            return;
        }

        size_t contig_cnt = -1u;
        if (args.size() > 2) {
            LOG("Will analyze first " << args[2] << " contigs");
            contig_cnt = std::stoll(args[2]);
        }


        io::FileReadStream reader(contigs_file);

        size_t i = 0;
        while (!reader.eof() && i < contig_cnt) {
            io::SingleRead contig;
            reader >> contig;
            LOG("Considering contig " << contig.name());

            std::vector<EdgeId> mapping_path = debruijn_graph::MapperInstance(curr_env.graph_pack())->MapRead(contig).simple_path();
            std::vector<std::vector<EdgeId>> splitted_path = Split(curr_env.graph_pack(), mapping_path);
            for(auto subpath : splitted_path) {
                if (HasCoverageDrops(curr_env.graph_pack(), subpath)) {
                    LOG("Has coverage drops, drawing");
                    DrawPicturesAlongPath(curr_env, subpath, contig.name());
                } else {
                    LOG("OK");
                }
            }

            ++i;
        }

    }

};

}
