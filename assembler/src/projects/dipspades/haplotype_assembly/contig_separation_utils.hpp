//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../utils/range_utils.hpp"
#include "../utils/path_routines.hpp"
#include "../utils/bulge_utils.hpp"
#include "conservative_regions_storage.hpp"
#include <string>

using namespace debruijn_graph;

namespace dipspades {

enum haplotype {unknown, from_one, from_different};

typedef map<pair<int, int>, haplotype>::iterator signed_label_iter;

class SignedLabels{

    map<pair<int, int>, haplotype> contigs_pairs;

public:
    void Add(int contig_id1, int contig_id2, haplotype new_label){
        pair<int, int> proc_pair(contig_id1, contig_id2);
        if(contigs_pairs.find(proc_pair) == contigs_pairs.end())
            contigs_pairs[proc_pair] = new_label;
        else{
            haplotype old_label = contigs_pairs[proc_pair];
            if(old_label < new_label)
                contigs_pairs[proc_pair] = new_label;
        }
    }

    haplotype GetHaplotypeByPair(int contig_id1, int contig_id2){
        return contigs_pairs[pair<int, int>(contig_id1, contig_id2)];
    }

    signed_label_iter begin(){
        return contigs_pairs.begin();
    }

    signed_label_iter end(){
        return contigs_pairs.end();
    }

    void MergeWith(SignedLabels new_signed_labels){
        for(auto it = new_signed_labels.begin(); it != new_signed_labels.end(); it++){
            Add(it->first.first, it->first.second, it->second);
        }
    }

    size_t Size(){
        return contigs_pairs.size();
    }

    string ToString(){
        stringstream ss;
        for(auto it = contigs_pairs.begin(); it != contigs_pairs.end(); it++)
            ss << "Pair " << it->first.first << ", " << it->first.second << " - " << it->second << ". ";
        return ss.str();
    }

    void WriteToFile(string fname, ContigStoragePtr contig_storage){
        ofstream out(fname.c_str());
        for(auto it= contigs_pairs.begin(); it != contigs_pairs.end(); it++)
            if(it->second == from_different){
                auto contig1 = contig_storage->GetContigById(it->first.first);
                auto contig2 = contig_storage->GetContigById(it->first.second);
                out << contig1->src_file() << ":" << contig1->name() << "\t" <<
                        contig2->src_file() << ":" << contig2->name() << endl;
            }
        out.close();
    }
};

class ContigLabelAllocator{
    ContigStoragePtr contig_storage_;

    Sequence GetSequenceByRange(Sequence seq, pair<size_t, size_t> r){
        return seq.Subseq(r.first, r.second);
    }

    bool AreRangesIntersect(MappingRange mapping_range1, MappingRange mapping_range2){

        Range mapped_range1 = mapping_range1.mapped_range;
        Range mapped_range2 = mapping_range2.mapped_range;

        if(!is_intersection_exist(mapped_range1, mapped_range2))
            return false;

        Range intersection = get_intersection_of_ranges(mapped_range1, mapped_range2);
        return intersection.end_pos - intersection.start_pos > 100;
    }

    haplotype ComputeLabelForPair(MappingRange mapping_range1, Sequence seq1,
            MappingRange mapping_range2, Sequence seq2){

        Range mapped_range1 = mapping_range1.mapped_range;
        Range mapped_range2 = mapping_range2.mapped_range;

        VERIFY(is_intersection_exist(mapped_range1, mapped_range2));

        TRACE("Mapping range1: " << mapped_range1.start_pos << " " <<
                mapped_range1.end_pos);
        TRACE("Mapping range2: " << mapped_range2.start_pos << " " <<
                mapped_range2.end_pos);

        TRACE("Init range1: " << mapping_range1.initial_range.start_pos << " " <<
                mapping_range1.initial_range.end_pos);
        TRACE("Init range2: " << mapping_range2.initial_range.start_pos << " " <<
                mapping_range2.initial_range.end_pos);


        Range intersection = get_intersection_of_ranges(mapped_range1, mapped_range2);

        TRACE("Intersection: " << intersection.start_pos << " " << intersection.end_pos);

        auto new_init_pair1 = project_init_range_to_new(mapped_range1, intersection,
                mapping_range1.initial_range);
        auto new_init_pair2 = project_init_range_to_new(mapped_range2, intersection,
                mapping_range2.initial_range);

        TRACE("1st projection: " << new_init_pair1.first << " " << new_init_pair1.second);
        TRACE("2nd projection: " << new_init_pair2.first << " " << new_init_pair2.second);

        if(!is_range_pair_correct(new_init_pair1) || !is_range_pair_correct(new_init_pair2))
            return unknown;

        Sequence subseq1 = GetSequenceByRange(seq1, new_init_pair1);
        Sequence subseq2 = GetSequenceByRange(seq2, new_init_pair2);

        double relative_align = RelAlignmentOfSequences(subseq1, subseq2);

        TRACE("Seq1 size - " << subseq1.size() << ", seq2 size - " << subseq2.size());
        TRACE("Relative alignment - " << relative_align);

        if(fabs(relative_align) < 0.0001)
            return from_one;
        return from_different;
    }

public:
    ContigLabelAllocator(ContigStoragePtr contig_storage) :
        contig_storage_(contig_storage) { }

    SignedLabels SignLabelsOnEdge(set<size_t> contigs, EdgeId current_edge){

        SignedLabels this_edge_labels;
        vector<int> indexes_of_edge;
        for(auto contig = contigs.begin(); contig != contigs.end(); contig++){
            int index = get_index_of_edge(contig_storage_->GetContigById(*contig)->
                    mapping_path().simple_path(), current_edge);
            VERIFY(index != -1);
            indexes_of_edge.push_back(index);
        }
        vector<int> oppa(contigs.begin(), contigs.end());
        for(size_t cnt1 = 0; cnt1 < oppa.size(); cnt1++) {
            int id1 = oppa[cnt1];
            auto seq1 = contig_storage_->GetContigById(id1)->seq();
            auto mapping_path1 = contig_storage_->GetContigById(id1)->mapping_path();
            MappingRange mapping_range1 = mapping_path1[indexes_of_edge[cnt1]].second;
            for(size_t cnt2 = cnt1 + 1; cnt2 < oppa.size(); cnt2++) {
                int id2 = oppa[cnt2];
                auto seq2 = contig_storage_->GetContigById(id2)->seq();
                auto mapping_path2 = contig_storage_->GetContigById(id2)->mapping_path();
                TRACE("Sign label for " << id1 << " " << id2);
                TRACE("Seq1 size - " << seq1.size() << " , seq2 size - " << seq2.size())
                MappingRange mapping_range2 = mapping_path2[indexes_of_edge[cnt2]].second;
                if(AreRangesIntersect(mapping_range1, mapping_range2)){
                    TRACE("Intersection exists");
                    haplotype label = ComputeLabelForPair(mapping_range1, seq1, mapping_range2, seq2);
                    this_edge_labels.Add(id1, id2, label);
                }
            }
        }

        return this_edge_labels;
    }

private:
    DECL_LOGGER("ContigLabelAllocator");
};

class IndexedPairOfEdges{
    pair<EdgeId, EdgeId> edges;
    pair<size_t, size_t> indexes;
    bool is_correct;

public:
    IndexedPairOfEdges(){
        is_correct = false;
    }

    IndexedPairOfEdges(EdgeId edge1, EdgeId edge2, size_t index1, size_t index2){
        edges.first = edge1;
        edges.second = edge2;

        indexes.first = index1;
        indexes.second = index2;

        is_correct = index1 <= index2;
    }

    EdgeId FirstEdge(){
        VERIFY(is_correct);
        return edges.first;
    }

    EdgeId SecondEdge(){
        VERIFY(is_correct);
        return edges.first;
    }

    size_t FirstIndex(){
        VERIFY(is_correct);
        return indexes.first;
    }

    size_t SecondIndex(){
        VERIFY(is_correct);
        return indexes.first;
    }

    bool IsNull(){
        return !is_correct;
    }
};

enum separation_result {not_identified, separated, diploid_repeat, conservative_region};

class SeparationResultInterpretator{

    bool IsConservativeRegion(SignedLabels labels){

        if(labels.Size() == 0)
            return false;

        for(auto it = labels.begin(); it != labels.end(); it++){
            if(it->second == from_different || it->second == unknown)
                return false;
        }

        return true;
    }

    bool AddNewEdgeIntoBigraph(set<int> &first_part, set<int> &second_part,
            pair<int,int> new_edge){

        int vertex1 = new_edge.first, vertex2 = new_edge.second;

        if(first_part.find(vertex1) != first_part.end() &&
                first_part.find(vertex2) != first_part.end())
            return false;

        if(second_part.find(vertex1) != second_part.end() &&
                second_part.find(vertex2) != second_part.end())
            return false;

        if(first_part.find(vertex1) != first_part.end()){
            second_part.insert(vertex2);
        }
        else{
            if(second_part.find(vertex1) != second_part.end())
                first_part.insert(vertex2);
            else{

                if(first_part.find(vertex2) != first_part.end())
                    second_part.insert(vertex1);
                else{
                    first_part.insert(vertex1);
                    second_part.insert(vertex2);
                }
            }
        }
        return true;
    }

    bool AreSeparatedContigs(SignedLabels labels){

        if(labels.Size() == 0)
            return false;

        set<int> first_part, second_part;

        for(auto it = labels.begin(); it != labels.end(); it++){
            if(it->second == from_different){
                pair<int, int> new_edge = it->first;
                if(!AddNewEdgeIntoBigraph(first_part, second_part, new_edge)){
                    TRACE("Edge doesn't added");
                    return false;
                }
            }
        }

        return true;
    }

public:
    separation_result Interpretate(SignedLabels labels){

        if(labels.Size() == 0){
            TRACE("Result unknown");
            return not_identified;
        }

        if(IsConservativeRegion(labels)){
            TRACE("Conservative region");
            return conservative_region;
        }

        if(AreSeparatedContigs(labels)){
            TRACE("Contigs are separated into two haplotypes");
            return separated;
        }

        TRACE("Diploid repeat");
        return diploid_repeat;
    }
};

class DiploidContigSeparator{

    typedef map<EdgeId, set<size_t> > EdgeContigsMap;

    Graph &g_;
    ContigStoragePtr default_storage_;
    ContigStoragePtr composite_storage_;
    CorrectionResult res_of_corr_cycle_;

    set<size_t> GetOfInnerContigsOf(size_t composite_contig_index){
        MappingContigPtr contig = (*composite_storage_)[composite_contig_index];
        set<size_t> contigs;
        vector<MappingContigPtr> inner_contigs = contig->AllMappingContigs();
        if(inner_contigs.size() == 0)
            contigs.insert(contig->id());
        else{
            for(auto it = inner_contigs.begin(); it != inner_contigs.end(); it++){
                size_t cur_id = (*it)->id();
                contigs.insert(cur_id);
                auto set_red_conts = res_of_corr_cycle_.redundancy_map.GetValuesByKey(cur_id);
                for(auto c = set_red_conts.begin(); c != set_red_conts.end(); c++)
                    contigs.insert(*c);
            }
        }
        return contigs;
    }

    set<EdgeId> GetSetOfEdgesFromPath(vector<EdgeId> path){
        set<EdgeId> res;
        for(auto e = path.begin(); e != path.end(); e++){
            res.insert(*e);
        }
        return res;
    }

    IndexedPairOfEdges DefineStartAndEndEdges(vector<EdgeId> common_path, MappingContigPtr inner_contig){
        MappingPath<EdgeId> map_path = inner_contig->mapping_path();
        VERIFY(map_path.size() > 0);
        EdgeId first_edge;
        size_t first_ind = size_t(-1);
        bool is_1st_found = false;
        for(size_t i = 0; i < map_path.size(); i++){
            for(size_t j = 0; j < common_path.size(); j++)
                if(map_path[i].first == common_path[j]){
                    first_edge = map_path[i].first;
                    first_ind = j;
                    is_1st_found = true;
                    break;
                }
            if(is_1st_found)
                break;
        }

        EdgeId last_edge;
        size_t last_ind = size_t(-1);
        bool is_2nd_found = false;
        for(int i = int(map_path.size() - 1); i >= 0; i--){
            for(int j = int(common_path.size()- 1); j >= 0; j--)
                if(map_path[i].first == common_path[j]){
                    last_edge = map_path[i].first;
                    last_ind = size_t(j);
                    is_2nd_found = true;
                    break;
                }
            if(is_2nd_found)
                break;
        }

        if(first_ind <= last_ind && is_1st_found && is_2nd_found)
            return IndexedPairOfEdges(first_edge, last_edge, first_ind, last_ind);
        else
            return IndexedPairOfEdges();
    }

    set<size_t> DeleteSubsetFromSet(set<size_t> set_, set<size_t> subset_){
        for(auto it = subset_.begin(); it != subset_.end(); it++)
            set_.erase(*it);
        return set_;
    }

    EdgeContigsMap DefineContigsOnEdges(set<size_t> contigs){
        EdgeContigsMap res;
        for(auto contig = contigs.begin(); contig != contigs.end(); contig++){
            auto map_path = default_storage_->GetContigById(*contig)->mapping_path();
            for(size_t i = 0; i < map_path.size(); i++)
                res[map_path[i].first].insert(*contig);
        }
        return res;
    }

    SignedLabels signed_labels_;
    ConservativeRegionStorage cons_regions_stor_;

public:
    DiploidContigSeparator(Graph &g, ContigStoragePtr default_storage,
        ContigStoragePtr composite_storage, CorrectionResult res_of_corr_cycle) :
        g_(g), default_storage_(default_storage), composite_storage_(composite_storage),
        res_of_corr_cycle_(res_of_corr_cycle){
    }

    void SeparateContigs(){

        SignedLabels signed_labels;
        ContigLabelAllocator label_allocator(default_storage_);

        // for each composite contig
        for(size_t i = 0; i < composite_storage_->Size(); i++){

            TRACE("New composite contig");

            // computing set of inner contigs
            set<size_t> inner_contigs = GetOfInnerContigsOf(i);

            TRACE("Number of contigs - " << inner_contigs.size());

            // define which contigs intersect consensus path
            vector<EdgeId> consensus_path = (*composite_storage_)[i]->path_seq();

            set<EdgeId> start_edge_edges_set;
            map<size_t, IndexedPairOfEdges> contig_start_end_map;
            set<size_t> contigs_for_deletion;

            for(auto c = inner_contigs.begin(); c != inner_contigs.end(); c++){
                MappingContigPtr contig = default_storage_->GetContigById(*c);
                auto edges = DefineStartAndEndEdges(consensus_path, contig);
                if(!edges.IsNull()){
                    contig_start_end_map[*c] = edges;
                    start_edge_edges_set.insert(edges.FirstEdge());
                    start_edge_edges_set.insert(edges.SecondEdge());
                }
                else
                    contigs_for_deletion.insert(*c);
            }

            inner_contigs = DeleteSubsetFromSet(inner_contigs, contigs_for_deletion);

            EdgeContigsMap contigs_on_edge = DefineContigsOnEdges(inner_contigs);

            TRACE("Defining labels");
            SeparationResultInterpretator interpret;
            for(auto e = consensus_path.begin(); e != consensus_path.end(); e++){

                TRACE("Edge - " << g_.str(*e) << ", start - " << g_.str(g_.EdgeStart(*e)) <<
                        ", end - " << g_.str(g_.EdgeEnd(*e)));
                auto contigs_ids_on_edge = contigs_on_edge[*e];

                if(contigs_ids_on_edge.size() == 1){
                    cons_regions_stor_.AddPossiblyConservativeRegion(g_.EdgeNucls(*e));
                    TRACE(g_.int_id(*e) << " - possibly conservative region");
                }

                TRACE("Contigs on this edge: " << SetToString<size_t>(contigs_ids_on_edge));

                SignedLabels current_signed_labels = label_allocator.SignLabelsOnEdge(contigs_ids_on_edge, *e);

                TRACE("Signed labels for this edge");
//                current_signed_labels.Print(cout);
                signed_labels_.MergeWith(current_signed_labels);

                TRACE("Interpretation of results");
                auto inpret_res = interpret.Interpretate(current_signed_labels);
                TRACE("------------------------------------------");

                if(inpret_res == conservative_region){
                    cons_regions_stor_.AddConservativeRegion(g_.EdgeNucls(*e));
                    TRACE(g_.int_id(*e) << " - conservative region");
                }
            }
        }

        TRACE("Signed labels:");
        TRACE(signed_labels_.ToString());

    }

    SignedLabels GetSignedLabels(){
        return signed_labels_;
    }

    ConservativeRegionStorage GetConservativeRegionStorage(){
        return cons_regions_stor_;
    }

private:
    DECL_LOGGER("DiploidContigSeparator");

};

}
