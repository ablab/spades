//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "io/wrapper_collection.hpp"

namespace online_visualization {

class DrawUnresolvedRepeatsCommand : public DrawingCommand {
private:
    static const size_t LENGTH_THRESHOLD = 1500;
    static const size_t GAP_DIFF_THRESHOLD = 1000;
    const string ref_prefix_;
    const double good_mapping_coeff_ = 0.7;

    void DrawGap(DebruijnEnvironment& curr_env, const vector<EdgeId>& path, string filename, string /*label*/ = "") const {
        omnigraph::visualization::WriteComponentsAlongPath<Graph>(curr_env.graph(), path, filename, curr_env.coloring(), curr_env.labeler());
        LOG("The pictures is written to " << filename);
    }

    vector<EdgePosition> GatherPositions(const GraphPack& gp, EdgeId e, const string& prefix) const {
        vector<EdgePosition> answer;
        for (EdgePosition pos : gp.edge_pos.GetEdgePositions(e)) {
            if (boost::starts_with(pos.contigId, prefix)) {
                answer.push_back(pos);
            }
        }
        return answer;
    }

    bool IsNext(MappingRange mr1, MappingRange mr2, size_t max_gap = 5) const {
       return mr2.initial_range.start_pos >= mr1.initial_range.end_pos && mr2.initial_range.start_pos <= mr1.initial_range.end_pos + max_gap;
    }

    shared_ptr<MappingRange> FindNextRange(MappingRange curr_range, const set<MappingRange>& ranges) const {
        cout << "Looking for next range for " << curr_range << endl;
        for (auto r : ranges) {
            cout << "Considering range " << r << endl;
            if (IsNext(curr_range, r)) {
                cout << "Found next" << endl;
                return make_shared<MappingRange>(r);
            }
        }
        cout << "Couldn't find suitable range" << endl;
        return shared_ptr<MappingRange>(0);
    }
    
    shared_ptr<pair<EdgeId, EdgePosition>> NextEdge(const GraphPack& gp, VertexId v, EdgePosition curr_pos) const {
        for (EdgeId next_e : gp.g.OutgoingEdges(v)) {
            cout << "Considering " << gp.g.str(next_e) << " as next edge " << endl;
            set<MappingRange> relevant_ranges = gp.edge_pos.GetEdgePositions(next_e, curr_pos.contigId);
            auto next_range = FindNextRange(curr_pos.mr, relevant_ranges);
            cout << "Considered " << relevant_ranges.size() << " relevant ranges" << endl;
            if (next_range) {
                cout << "Found next edge" << endl;
                return make_shared<pair<EdgeId, EdgePosition>>(next_e, EdgePosition(curr_pos.contigId, *next_range));
            }
        }
        cout << "Couldn't find next edge" << endl;
        return shared_ptr<pair<EdgeId, EdgePosition>>(0);
    }

    vector<EdgeId> FindReferencePath(const GraphPack& gp, EdgeId e1, EdgeId e2) const {
        EdgePosition curr_pos = GatherPositions(gp, e1, ref_prefix_).front();
        VertexId curr_v = gp.g.EdgeEnd(e1);
        vector<EdgeId> answer;
        answer.push_back(e1);
    //    for (size_t i = 0 ; i < 1000 ; ++i) {
        while(true) {
            auto next_info = NextEdge(gp, curr_v, curr_pos);
            if (next_info) {
                EdgeId next_e = next_info->first;
                answer.push_back(next_e);
                if (next_e == e2) {
                    break;
                }
                curr_v = gp.g.EdgeEnd(next_e);
                curr_pos = next_info->second;
            } else {
                return vector<EdgeId>();
            }
        }
        return answer;
    }

    MappingPath<EdgeId> FindReferencePath2(const GraphPack& gp, EdgeId e1, EdgeId e2) const {
        auto e1_poss = GatherPositions(gp, e1, ref_prefix_);
        auto e2_poss = GatherPositions(gp, e2, ref_prefix_);
        VERIFY(e1_poss.size() == 1 && e2_poss.size() == 1);
        EdgePosition e1_pos = e1_poss.front();
        EdgePosition e2_pos = e2_poss.front();
        VERIFY(e1_pos.contigId == e1_pos.contigId);
        Sequence ref = (e1_pos.contigId == "ref0") ? gp.genome : !gp.genome;
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

    bool AnalyzeGaps(DebruijnEnvironment& curr_env, const GraphPack& gp, io::SingleRead contig, string base_assembly_prefix) const {
        auto mapper_ptr = debruijn_graph::MapperInstance(gp);
        MappingPath<EdgeId> mapping_path = mapper_ptr->MapRead(contig);
        auto pos_handler = gp.edge_pos;

        auto simple_path = mapping_path.simple_path();
        vector<size_t> long_unique_idx;
        for (size_t i = 0; i < simple_path.size(); ++i) {
            EdgeId e = simple_path[i];
            if (gp.g.length(e) >= LENGTH_THRESHOLD 
                    && IsOfMultiplicityOne(gp, e)
                    && CheckMapping(gp.g.length(e), mapping_path[i].second.mapped_range.size())) {
                long_unique_idx.push_back(i);
            }
        }

        bool found_smth = false;
        size_t cnt = 0;
        for (size_t i = 1; i < long_unique_idx.size(); ++i) {
            EdgeId e1 = simple_path[long_unique_idx[i-1]];
            EdgeId e2 = simple_path[long_unique_idx[i]];
            if (!BelongToSameContig(gp, e1, e2, base_assembly_prefix)) {
                INFO("Long unique edges not in the same contig of base assembly");

                vector<EdgeId> contig_gap_path;
                for (size_t j = long_unique_idx[i-1] + 1; j < long_unique_idx[i-1]; ++j) {
                    contig_gap_path.push_back(simple_path[j]);
                }
                //checking correctness of path is too restricting, check genomic distance instead
                size_t genomic_gap = GenomicGap(gp, e1, e2);
                size_t contig_gap = mapping_path[i].second.initial_range.start_pos - mapping_path[i-1].second.initial_range.end_pos;
                if (genomic_gap != -1u && std::abs(int(genomic_gap) - int(contig_gap)) < GAP_DIFF_THRESHOLD) {
                    INFO("Found genomic gap " << genomic_gap << 
                        " between e1 " << gp.g.str(e1) << " genome pos : " << GatherPositions(gp, e1, "ref") << " and e2 " << gp.g.str(e2)
                        << " genome pos : " << GatherPositions(gp, e2, "ref"));
                    INFO("Looking for reference path");
                    auto ref_mapping_path = FindReferencePath2(gp, e1, e2);
                    if (ref_mapping_path.size() == 0) {
                        WARN("Couldn't find ref path between " << gp.g.str(e1) << " and " << gp.g.str(e2));
                    } else {
                        found_smth = true;
                        
                        vector<EdgeId> ref_path = mapper_ptr->FindReadPath(ref_mapping_path);
                        if (ref_path.empty()) {
                            WARN("Couldn't fix ref path");
                            ref_path = ref_mapping_path.simple_path();
                        }
                        INFO("Found ref path between " << gp.g.str(e1) << " and " << gp.g.str(e2));
                        INFO(ref_path.size() << " edges of cumulative length " << CumulativeLength(gp.g, ref_path));
                        
                        make_dir(curr_env.folder());
                        string pics_folder = curr_env.folder() + "/" + curr_env.GetFormattedPictureCounter()  + "_" + contig.name() + "/";
                        make_dir(pics_folder);
                        string pic_name = ToString(cnt++) + "_" +  ToString(GenomicGap(gp, e1, e2)) + "_" + ToString(gp.g.int_id(e1)) + "_" + ToString(gp.g.int_id(e2)) + "_";
                        DrawGap(curr_env, ref_path, pics_folder + pic_name);
                    }
                } else {
                    if (genomic_gap == -1u) {
                        WARN("Contig likely misassembled. Unique long regions in wrong order. e1 " << 
                                gp.g.str(e1) << " genome pos : " << GatherPositions(gp, e1, "ref") << " and e2 " << gp.g.str(e2) << 
                                " genome pos : " << GatherPositions(gp, e2, "ref"));
                    } else {
                        WARN("Contig likely misassembled. Genomic gap is " << genomic_gap << " while contig gap was " << contig_gap);
                    }
                }
            }
        }
        return found_smth;
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
        answer = answer + "Command `draw_unresolved_repeats` \n" + "Usage:\n"
                + "> draw_unresolved_repeats <contigs_file> <prefix_of_base_assembly> [first N contigs to analyze]\n"
                + " Draws pictures of unresolved repeats.";
        return answer;
    }

    DrawUnresolvedRepeatsCommand()
            : DrawingCommand("draw_unresolved_repeats"), ref_prefix_("ref"), good_mapping_coeff_(0.7) {
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        std::string contigs_file = args[1];
        string base_assembly_prefix = args[2];

        if (!CheckFileExists(contigs_file)) {
            LOG("File with contigs " << contigs_file << " not found");
        }

        size_t contig_cnt = -1u;
        if (args.size() > 3) {
            LOG("Will analyze first " << args[3] << " contigs");
            contig_cnt = lexical_cast<size_t>(args[3]);
        }
        
        auto reader = make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(contigs_file));

        size_t i = 0;
        while (!reader->eof() && i < contig_cnt) {
            io::SingleRead contig;
            (*reader) >> contig;
            LOG("Considering contig " << contig.name());

            if (AnalyzeGaps(curr_env, curr_env.graph_pack(), contig, base_assembly_prefix)) {
                curr_env.inc_pic_counter();
            }
            ++i;
        }

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
        }

        size_t contig_cnt = -1u;
        if (args.size() > 3) {
            LOG("Will analyze first " << args[3] << " contigs");
            contig_cnt = lexical_cast<size_t>(args[3]);
        }
        
        auto reader = make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(contigs_file));

        size_t i = 0;
        while (!reader->eof() && i < contig_cnt) {
            io::SingleRead contig;
            (*reader) >> contig;
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
}
