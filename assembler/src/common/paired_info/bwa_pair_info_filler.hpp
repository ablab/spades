//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"

#include "pipeline/config_struct.hpp"
#include <paired_info/paired_info.hpp>

#ifndef PROJECT_BWA_PAIR_INFO_FILLER_HPP_H
#define PROJECT_BWA_PAIR_INFO_FILLER_HPP_H

namespace bwa_pair_info {

using debruijn_graph::EdgeId;

typedef omnigraph::de::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> PairedInfoIndexT;
typedef io::SequencingLibrary<debruijn_graph::config::DataSetData> SequencingLibraryT;
typedef std::pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId> EdgePair;
typedef unordered_map<size_t, debruijn_graph::EdgeId> EdgeIdMap;

//Base class for aligned read processor (simple analog of SequenceMapperListener)
class MapperReadT;

class BWAPairedReadProcessor {
public:
    virtual void ProcessPairedRead(const MapperReadT& l, const MapperReadT& r) = 0;
    virtual ~BWAPairedReadProcessor() {}
};

//Class for running BWA, managing and parsing SAM files
class BWAPairInfoFiller {
public:
    DECL_LOGGER("BWAPairInfo");

private:
    const debruijn_graph::Graph& g_;

    string bwa_path_;
    string base_dir_;
    string work_dir_;
    size_t nthreads_;
    string index_base_;
    bool index_constructed_;
    bool remove_tmp_files_;
    unordered_map<size_t, debruijn_graph::EdgeId> edge_id_map_;

private:

    //Save graph in fasta format
    void OutputEdges(const string& filename) const;

    //Construct int_id -> EdgeId map
    void FillEdgeIdMap();

    //Run bwa index
    bool CreateIndex(const string& contigs);

    //Initialize for read aligment (includes all above)
    bool Init();

    //Run bwa mem on single file
    bool RunBWA(const string& reads_file, const string& out_sam_file) const;

    //Process single read library
    bool AlignLib(const SequencingLibraryT& lib,
                      const string& sam_file_base,
                      vector<pair<string, string>>& resulting_sam_files);

    //Parse a pair of same files and analyze alignments with processor
    void ProcessSAMFiles(const string &left_sam, const string &right_sam,
                         BWAPairedReadProcessor& processor);

public:

    BWAPairInfoFiller(const debruijn_graph::Graph& g,
                      const string& bwa_path,
                      const string& work_dir,
                      size_t nthreads = 1,
                      bool remove_tmp = true):
        g_(g), bwa_path_(bwa_path), base_dir_(work_dir), work_dir_(""),
        nthreads_(nthreads), index_base_(""), index_constructed_(false),
        remove_tmp_files_(remove_tmp),
        edge_id_map_() {
    }

    ~BWAPairInfoFiller() {
        if (remove_tmp_files_)
            path::remove_if_exists(work_dir_);
    }

    //Count IS and fill pair info index for the given lib
    bool ProcessLib(size_t lib_index,
                    SequencingLibraryT& lib,
                    PairedInfoIndexT& paired_index,
                    size_t counter_edge_len,
                    size_t index_filler_edge_len);
};

}

#endif //PROJECT_BWA_PAIR_INFO_FILLER_HPP_H
