//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/alignment/sequence_mapper.hpp"
#include "contig_separation_utils.hpp"

using namespace debruijn_graph;

namespace dipspades {

class ConservativeRegionsSearcher{
    conj_graph_pack & dbl_gp_;
    ContigStoragePtr storage_;
    SignedLabels signed_labels_;
    ConservativeRegionStorage cons_reg_storage_;

    BasicSequenceMapper<conj_graph_pack::graph_t, conj_graph_pack::index_t> mapper_;
    map<int, MappingPath<EdgeId> > contig_map_path_;

    typedef map<int, vector<int> > diff_labeled_contigs;
    diff_labeled_contigs map_of_diff_contigs_;

    MappingPath<EdgeId> GetMappingPath(int contig){
        if(contig_map_path_.find(contig) == contig_map_path_.end()){
            auto seq = storage_->GetContigById(contig)->seq();
            MappingPath<EdgeId> map_path = mapper_.MapSequence(seq);
            contig_map_path_[contig] = map_path;
        }
        return contig_map_path_[contig];
    }

    void ComputeDifferentLabeledContigs(){
        for(auto it = signed_labels_.begin(); it != signed_labels_.end(); it++)
            if(it->second == from_different){
                int contig1 = it->first.first, contig2 = it->first.second;
                map_of_diff_contigs_[contig1].push_back(contig2);
            }
    }

    vector<EdgeId> GetConservativeEdges(vector<MappingPath<EdgeId> > paths,
            vector<int> labels){

        map<EdgeId, set<int> > edge_labels;
        for(size_t i = 0; i < paths.size(); i++){
            MappingPath<EdgeId> path = paths[i];
            int label = labels[i];

            for(size_t j = 0; j < path.size(); j++)
                edge_labels[path[j].first].insert(label);
        }

//        for(auto it = edge_labels.begin(); it != edge_labels.end(); it++){
//            cout << it->first << ". Labels ";
//            PrintSet<int>(cout, it->second);
//        }

        vector<EdgeId> cons_edges;
        for(auto it = edge_labels.begin(); it != edge_labels.end(); it++)
            if(it->second.size() > 1)
                cons_edges.push_back(it->first);

        return cons_edges;
    }

    vector<int> GatherLabelsForSeparatedContigs(vector<int> separated_contigs){
        vector<int> labels;
        labels.push_back(1);
        for(auto c = separated_contigs.begin(); c != separated_contigs.end(); c++){
            labels.push_back(2);
        }
        return labels;
    }

    vector<MappingPath<EdgeId> > GatherMappingPathForContigs(vector<int> contigs){
        vector<MappingPath<EdgeId> > map_paths;
        for(auto c = contigs.begin(); c != contigs.end(); c++)
            map_paths.push_back(GetMappingPath(*c));
        return map_paths;
    }

    void FindTwoColoredEdges(){

        ComputeDifferentLabeledContigs();

        for(auto it = map_of_diff_contigs_.begin(); it != map_of_diff_contigs_.end(); it++){

            auto contig = it->first;
            auto separated_contigs = it->second;

//            cout << contig << ". Separated set - ";
//            PrintVector<int>(cout, separated_contigs);

            auto labels = GatherLabelsForSeparatedContigs(separated_contigs);

            // gather all mapping paths
            auto contig_map_path = GetMappingPath(contig);
            vector<MappingPath<EdgeId> > map_paths = GatherMappingPathForContigs(separated_contigs);
            map_paths.insert(map_paths.begin(), contig_map_path);

            // find two or more colored edges
            auto cur_cons_edges = GetConservativeEdges(map_paths, labels);

            // add them in storage
            for(auto e = cur_cons_edges.begin(); e != cur_cons_edges.end(); e++)
                cons_reg_storage_.AddConservativeRegion(dbl_gp_.g.EdgeNucls(*e));
        }
    }

    void WriteConservativeRegionsStorageToFile(string filename, cons_regions_iterator iter_begin,
            cons_regions_iterator iter_end){
        ofstream fout(filename);
        int cnt = 1;
        for(auto it = iter_begin; it != iter_end; it++){
            Sequence curr_seq = *it;
            fout  << ">" << cnt << "_conservative_region_length_" << curr_seq.size() << endl;
            fout << curr_seq.str() << endl;
            cnt++;
        }
    }

    size_t ComputeSummaryLengthOfRegionInStorage(cons_regions_iterator iter_begin,
            cons_regions_iterator iter_end){
        size_t summary_cons_reg_length = 0;
        for(auto it = iter_begin; it != iter_end; it++){
            summary_cons_reg_length += it->size();
        }
        return summary_cons_reg_length;
    }

public:
    ConservativeRegionsSearcher(conj_graph_pack &dbl_gp,  ContigStoragePtr storage,
            SignedLabels signed_labels, ConservativeRegionStorage cons_reg_storage) :
                dbl_gp_(dbl_gp),
                storage_(storage),
                signed_labels_(signed_labels),
                cons_reg_storage_(cons_reg_storage),
                mapper_(dbl_gp_.g, dbl_gp_.index,
                        dbl_gp_.kmer_mapper) { }

    void Search(){
        FindTwoColoredEdges();
        size_t cons_regions_length = ComputeSummaryLengthOfRegionInStorage(cons_reg_storage_.cons_regions_begin(),
                cons_reg_storage_.cons_regions_end());
        if(cons_regions_length > 0){
            string cons_regions_fname(path::append_path(dsp_cfg::get().io.output_dir,
                    "conservative_regions.fasta").c_str());
            WriteConservativeRegionsStorageToFile(cons_regions_fname, cons_reg_storage_.cons_regions_begin(),
                cons_reg_storage_.cons_regions_end());
            INFO("Conservative regions with total length " << cons_regions_length <<
                    " written in file " << cons_regions_fname);
        }

        size_t poss_cons_regions_length = ComputeSummaryLengthOfRegionInStorage(cons_reg_storage_.poss_cons_regions_begin(),
                cons_reg_storage_.poss_cons_regions_end());
        if(poss_cons_regions_length > 0){
            string poss_cons_regions_fname(path::append_path(dsp_cfg::get().io.output_dir,
                "possibly_conservative_regions.fasta").c_str());
//            INFO("Possibly conservative regions written in file " << poss_cons_regions_fname);
            WriteConservativeRegionsStorageToFile(poss_cons_regions_fname, cons_reg_storage_.poss_cons_regions_begin(),
                cons_reg_storage_.poss_cons_regions_end());
            INFO("Conservative regions with total length " << poss_cons_regions_length <<
                    " written in file " << poss_cons_regions_fname);
        }
    }
};

}

