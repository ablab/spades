//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/reads/io_helper.hpp"

#include "utils/element_printers.hpp"
#include "utils/files_utils.hpp"

#include "contig_correctors/close_gaps_corrector.hpp"
#include "contig_correctors/iterative_redundant_contigs_remover.hpp"
#include "contig_correctors/overlap_searcher.hpp"
#include "contig_correctors/same_edge_deletion_corrector.hpp"
#include "contig_correctors/incorrect_contig_remover.hpp"

using namespace debruijn_graph;

namespace dipspades{

class ConsensusContigsConstructor {
    conj_graph_pack &graph_pack_;
    BaseHistogram<size_t> &bulge_len_hist_;
    BasicSequenceMapper<conj_graph_pack::graph_t, conj_graph_pack::index_t> seq_mapper_;
    VertexPathIndex path_index_;

    CorrectionResult correction_result_;
    ContigStoragePtr default_storage_;
    ContigStoragePtr composite_storage_;

    set<size_t> ind_zero_paths_;

    struct contig_name {
        string fname;
        string name;

        contig_name(string new_fname, string new_name) :
            fname(cut_fname_from_path(new_fname)),
            name(new_name) { }
    };

    typedef pair<contig_name, Sequence> contig;

    vector<contig> ReadContigs(string contig_fname){
        vector<contig> contigs;
        auto fstream = io::SplittingWrap(EasyStream(contig_fname, false));

        while(!fstream->eof()){
            SingleRead single_read;
            (*fstream) >> single_read;
            contigs.push_back(contig(contig_name(contig_fname, single_read.name()),
                    single_read.sequence()));
        }
        INFO(contigs.size() << " contigs from " << contig_fname << " were read");
        return contigs;
    }

    vector<contig> ReadContigsFromFiles(vector<string> contig_fnames){
        vector<contig> contigs;
        for(auto it = contig_fnames.begin(); it != contig_fnames.end(); it++){
            if(fname_valid(*it)){
                auto contigs_from_file = ReadContigs(*it);
                contigs.insert(contigs.end(), contigs_from_file.begin(), contigs_from_file.end());
            }
        }
        return contigs;
    }

    vector<MappingPath<EdgeId> > ConstructMappathsWithoutRC(vector<contig> &contigs){
        vector<MappingPath<EdgeId> > map_paths;
        size_t zero_paths = 0;
        size_t total_length_unmapped = 0;
        for(size_t i = 0; i < contigs.size(); i++){
            map_paths.push_back(seq_mapper_.MapSequence(contigs[i].second));
            if(map_paths[map_paths.size() - 1].size() == 0){
                total_length_unmapped += contigs[i].second.size();
                zero_paths++;
            }
        }
        if(zero_paths != 0)
            INFO(ToString(zero_paths) + " contigs with total length " << total_length_unmapped <<
                    " have mapped path with zero length");
        return map_paths;
    }

    void DefineIndicesOfZeroPaths(vector<MappingPath<EdgeId> > &map_paths){
        for(size_t i = 0; i < map_paths.size(); i++)
            if(map_paths[i].size() == 0)
                ind_zero_paths_.insert(i);
    }

    ContigStoragePtr CreateContigStorage(vector<contig> &contigs,
            vector<MappingPath<EdgeId> > &mapping_paths) {
        ContigStoragePtr default_storage(new SimpleContigStorage());
        for(size_t i = 0; i < mapping_paths.size(); i++){
            if(ind_zero_paths_.find(i) == ind_zero_paths_.end()){
                int i1 = int(i * 2), i2 = int(i * 2 + 1);
                default_storage->Add(MappingContigPtr(
                        new SimpleMappingContig(contigs[i].first.name,
                        contigs[i].first.fname, contigs[i].second,
                        mapping_paths[i], i1, i2)));
            }
        }
        return default_storage;
    }

    void PrimaryContigsProcessing(ContigStoragePtr storage){
        INFO("Removing repetitive edges in contigs mapping starts");
        SameEdgeDeletionCorrector same_edges_corr(graph_pack_.g);
        same_edges_corr.Correct(storage);
//            INFO(storage->Size() << " contigs will be used");
        INFO("Removing repetitive edges in contigs mapping ends");

        INFO("Close gaps in contigs mappings starts")
        CloseGapsCorrector close_gaps_corr(graph_pack_.g);
        close_gaps_corr.Correct(storage);
//           INFO(storage->Size() << " contigs will be used");
        INFO("Close gaps in contigs mappings ends");

        INFO("Removing incorrect contigs")
        RemoveUnconnectContigsCorrector del_unconn_corr(graph_pack_.g);
        del_unconn_corr.Correct(storage);
//           INFO(storage->Size() << " contigs will be used");
    }

    string name_to_rc_name(string name){
        return name + "_RC";
    }

    ContigStoragePtr CreateStorageWithRCContigs(ContigStoragePtr old_storage){
        ContigStoragePtr new_storage(new SimpleContigStorage());
        TRACE("CreateStorageWithRCContigs starts");
        for(size_t i = 0; i < old_storage->Size(); i++){
            auto contig = (*old_storage)[i];
            new_storage->Add(contig);

            MappingContigPtr rc_contig = MappingContigPtr(
                    new SimpleMappingContig(
                            name_to_rc_name(contig->name()),
                            contig->src_file(),
                            !contig->seq(),
                            GetRCToMappingPath(graph_pack_.g, contig->mapping_path(), contig->seq().size()),
                            GetRCToPathSeq(graph_pack_.g, contig->path_seq()),
                            contig->id() + 1, contig->id()));
            new_storage->Add(rc_contig);
        }
        TRACE("CreateStorageWithRCContigs ends");
        INFO("Addition of RC contigs. " << new_storage->Size() << " contigs will be used");
        return new_storage;
    }

    void RemoveRedundantContigs(ContigStoragePtr storage){
        INFO("Redundant contigs remover starts");
        VertexPathIndex path_index(graph_pack_.g);
        IterativeLoopCorrector iter_loop_corr(
                graph_pack_.g,
                graph_pack_.k_value,
                path_index,
                dsp_cfg::get().cc.max_loop_length,
                dsp_cfg::get().cc.min_lcs_size,
                dsp_cfg::get().cc.estimate_tails ?
                        bulge_len_hist_.Quantile(dsp_cfg::get().cc.bulge_len_quantile) :
                        dsp_cfg::get().pbr.max_bulge_nucls_len);
        iter_loop_corr.Correct(storage);
        INFO("Redundant contigs remover ends");
        correction_result_ = iter_loop_corr.Results();
    }

    ContigStoragePtr DefineOverlappingContigs(ContigStoragePtr storage){
        INFO("Overlapping search starts");
        path_index_.Initialize(storage);
        OverlapCorrector over_corr(graph_pack_.g,
                graph_pack_.k_value,
                dsp_cfg::get().cc.min_overlap_size,
                path_index_);
        auto new_storage = over_corr.Correct(storage);
        path_index_.Clear();
        INFO("Overlapping search ends");
        return new_storage;
    }

    void WriteContigsToFile(ContigStoragePtr contigs, string filename){
        size_t total_length = 0;
        ofstream out(filename);
        for(size_t i = 0; i < contigs->Size(); i++){
            vector<EdgeId> contig_path = (*contigs)[i]->path_seq();
            TRACE(i << " path: " << SimplePathWithVerticesToString(graph_pack_.g, contig_path));
            Sequence seq = (*contigs)[i]->seq();
            out << ">" << (*contigs)[i]->name() << endl;
            out << seq.str() << endl;
            total_length += seq.size();
        }
        INFO(contigs->Size() << " with total length " << total_length << " were written in " <<
                filename);
    }

    void WritePairedAndUnpairedContigs(ContigStoragePtr storage){
        ContigStoragePtr double_contigs(new SimpleContigStorage());
        ContigStoragePtr single_contigs(new SimpleContigStorage());
        for(size_t i = 0; i < storage->Size(); i++){
            auto contig = (*storage)[i];
            if(contig->AllMappingContigs().size() == 0){
                if(correction_result_.redundancy_map.GetValuesByKey(contig->id()).size() == 0)
                    single_contigs->Add(contig);
                else
                    double_contigs->Add(contig);
            }
            else
                double_contigs->Add(contig);
        }
        WriteContigsToFile(double_contigs,
                path::append_path(dsp_cfg::get().io.output_dir, "paired_consensus_contigs.fasta").c_str());
        WriteContigsToFile(single_contigs,
                path::append_path(dsp_cfg::get().io.output_dir, "unpaired_consensus_contigs.fasta").c_str());
    }

    void WriteAlignedHaplocontigs(){
        string fname = path::append_path(dsp_cfg::get().io.output_dir, "haplocontigs_alignment");
        ofstream out(fname.c_str());
        INFO("Writing haplocontigs alignment to " << fname);

        for(size_t i = 0; i < composite_storage_->Size(); i++){
            auto composite_contig = (*composite_storage_)[i];
            out << "Consensus contig: " << composite_contig->name() << endl;
            auto haplocontigs = composite_contig->AllMappingContigs();
            if(haplocontigs.size() == 0) // contig is not composite
                haplocontigs.push_back(composite_contig);

            if(haplocontigs.size() > 1){
                out << "\tOverlapped haplocontigs: " << endl;
                for(size_t i = 0; i < haplocontigs.size() - 1; i++)
                    out << "\t\t" << haplocontigs[i]->full_name() << "\t" <<
                    haplocontigs[i + 1]->full_name() << endl;
            }

            out << "\tAligned pairs: " << endl;
            size_t written_pairs = 0;
            for(auto h = haplocontigs.begin(); h != haplocontigs.end(); h++){
                size_t id = (*h)->id();
                auto redundant_contigs = correction_result_.redundancy_map.GetValuesByKey(id);
                for(auto it = redundant_contigs.begin(); it != redundant_contigs.end(); it++){
                    out << "\t\t" << (*h)->full_name() << "\t" <<
                            default_storage_->GetContigById(*it)->full_name() << endl;
                    written_pairs++;
                }
            }

            if(written_pairs == 0)
                out << "\t\tNo pairs" << endl;
        }

/*        for(auto it = correction_result_.redundancy_map.begin();
                it != correction_result_.redundancy_map.end(); it++){
            auto contig1 = default_storage_->GetContigById(it->first);
            auto set_ids = it->second;
            for(auto set_it = set_ids.begin(); set_it != set_ids.end(); set_it++){
                auto contig2 = default_storage_->GetContigById(*set_it);
                out << contig1->src_file() << ":" << contig1->name() << "\t" <<
                        contig2->src_file() << ":" << contig2->name() << endl;
            }
        }*/

    }

public:
    ConsensusContigsConstructor(conj_graph_pack &graph_pack,
            BaseHistogram<size_t> &bulge_len_hist) :
            graph_pack_(graph_pack),
            bulge_len_hist_(bulge_len_hist),
            seq_mapper_(graph_pack.g, graph_pack.index,
                    graph_pack.kmer_mapper, false),
            path_index_(graph_pack.g),
            correction_result_(),
            default_storage_(),
            composite_storage_() { }

    void Run() {
        INFO("Consensus contigs constructor starts");
        auto contigs = ReadContigsFromFiles(GetAllLinesFromFile(dsp_cfg::get().io.haplocontigs));
        INFO("Total: " << contigs.size() << " contigs were read");
        if(contigs.size() == 0)
            return;

        vector<MappingPath<EdgeId> > mapping_paths = ConstructMappathsWithoutRC(contigs);
        VERIFY(mapping_paths.size() == contigs.size());
        DefineIndicesOfZeroPaths(mapping_paths);

        auto preliminary_storage = CreateContigStorage(contigs, mapping_paths);

        TRACE("Preliminary storage:");
        TRACE(preliminary_storage->ToString(graph_pack_.g));

        PrimaryContigsProcessing(preliminary_storage);

        TRACE("Preliminary storage after 1st processing:");
        TRACE(preliminary_storage->ToString(graph_pack_.g));

        auto processed_storage = CreateStorageWithRCContigs(preliminary_storage);
        VERIFY(processed_storage->Size() % 2 == 0);

        default_storage_ = processed_storage->Clone();
        RemoveRedundantContigs(processed_storage);

        TRACE("Storage after removing redundant contigs:");
        TRACE(processed_storage->ToString(graph_pack_.g));

        composite_storage_ = DefineOverlappingContigs(processed_storage);

        string consensus_fname(path::append_path(dsp_cfg::get().io.output_dir, "consensus_contigs.fasta").c_str());
        WriteContigsToFile(composite_storage_, consensus_fname);
        WritePairedAndUnpairedContigs(composite_storage_);

        WriteAlignedHaplocontigs();

        INFO("Consensus contigs constructor ends");
    }

    ContigStoragePtr DefaultContigsStorage() { return default_storage_; }

    ContigStoragePtr CompositeContigsStorage() { return composite_storage_; }

    CorrectionResult RedundancyResult() { return correction_result_; }

private:
    DECL_LOGGER("ConsensusContigsConstructor");
};

}
