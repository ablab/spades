//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sam/sam_reader.hpp"
#include "standard.hpp"
#include "debruijn_graph.hpp"
#include "config_struct.hpp"

#include <io/osequencestream.hpp>
#include <de/paired_info.hpp>

#ifndef PROJECT_BWA_PAIR_INFO_FILLER_HPP_H
#define PROJECT_BWA_PAIR_INFO_FILLER_HPP_H

namespace bwa_pair_info {

using namespace sam_reader;

class BWAPairInfoFiller {
public:

    typedef std::pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId> EdgePair;

    typedef omnigraph::de::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> PairedInfoIndexT;

    typedef io::SequencingLibrary<debruijn_graph::debruijn_config::DataSetData> SequencingLibraryT;

private:

    class MapperReadT {
    public:
        MapperReadT(): contig_id_(""), pos_(-1), len_(-1), is_forward_(true),
                       left_hard_clip_(0), right_hard_clip_(0), left_soft_clip_(0), right_soft_clip_(0){}

        MapperReadT(const string& ctg_id, int32_t pos, int32_t len, bool is_forward, const string& cigar):
            contig_id_(ctg_id), pos_(pos), len_(len), is_forward_(is_forward),
            left_hard_clip_(0), right_hard_clip_(0), left_soft_clip_(0), right_soft_clip_(0) {

            ParseCigar(cigar);
        }

        bool IsValid() const {
            return contig_id_ != "";
        }

    private:

        void ParseCigar(const string& cigar);

    public:
        const string &get_contig_id() const {
            return contig_id_;
        }
        int32_t get_pos() const {
            return pos_;
        }
        int32_t get_len() const {
            return len_;
        }
        bool is_forward() const {
            return is_forward_;
        }
        uint32_t get_left_soft_clip() const {
            return left_soft_clip_;
        }
        uint32_t get_right_soft_clip() const {
            return right_soft_clip_;
        }
        uint32_t get_left_hard_clip() const {
            return left_hard_clip_;
        }
        uint32_t get_right_hard_clip() const {
            return right_hard_clip_;
        }

    private:
        string contig_id_;
        int32_t pos_;
        int32_t len_;
        bool is_forward_;
        uint32_t left_hard_clip_:16, right_hard_clip_:16;
        uint32_t left_soft_clip_:16, right_soft_clip_:16;
    };

private:
    const debruijn_graph::Graph& g_;

    string bwa_path_;

    string work_dir_;

    size_t min_contig_len_;

    size_t nthreads_;

    string index_base_;

    bool index_constructed_;

    unordered_map<size_t, debruijn_graph::EdgeId> edge_id_map_;

private:
    EdgePair ConjugatePair(EdgePair ep) const;

    void OutputEdges(const string& filename) const;

    void FillEdgeIdMap();

    bool CreateIndex(const string& contigs);

    bool CreateIndexIfNeeded(const string& contigs);


    bool RunBWA(const string& reads_file, const string& out_sam_file) const;

    bool AlignLib(const SequencingLibraryT& lib,
                      const string& sam_file_base,
                      vector<pair<string, string>>& resulting_sam_files);


    void ProcessPairedRead(const SequencingLibraryT& lib,
                               const MapperReadT& l, const MapperReadT& r,
                               PairedInfoIndexT& paired_index);

    void ParseSAMFiles(const SequencingLibraryT& lib,
                           const string& left_sam, const string& right_sam,
                           PairedInfoIndexT& paired_index);

public:

    BWAPairInfoFiller(const debruijn_graph::Graph& g,
                      const string& bwa_path,
                      const string& work_dir,
                      size_t min_contig_len = 1000,
                      size_t nthreads = 1):
        g_(g), bwa_path_(bwa_path), work_dir_(work_dir),
        min_contig_len_(min_contig_len), nthreads_(nthreads),
        index_base_(""), index_constructed_(false),
        edge_id_map_() {

        Init();
    }

    bool Init();

    ~BWAPairInfoFiller() {
        //TODO: uncomment
        //path::remove_if_exists(work_dir_);
    }

    bool ProcessLib(size_t lib_index,
                        const SequencingLibraryT& lib,
                        PairedInfoIndexT& paired_index);
};



struct ContigOrienter {
    typedef debruijn_graph::EdgeId EdgeId;

    const debruijn_graph::Graph& g_;

    ContigOrienter(const debruijn_graph::Graph& g): g_(g) {}

    virtual pair<EdgeId, EdgeId> Apply(EdgeId e1, EdgeId e2) const = 0;

    virtual ~ContigOrienter() {}
};

struct FFOrienter: public ContigOrienter {
    FFOrienter(const debruijn_graph::Graph& g): ContigOrienter(g) {}

    virtual pair<EdgeId, EdgeId> Apply(EdgeId e1, EdgeId e2) const {
        return make_pair(e1, e2);
    }
};

struct FROrienter: public ContigOrienter {
    FROrienter(const debruijn_graph::Graph& g): ContigOrienter(g) {}

    virtual pair<EdgeId, EdgeId> Apply(EdgeId e1, EdgeId e2) const {
        return make_pair(e1, g_.conjugate(e2));
    }
};

struct RFOrienter: public ContigOrienter {
    RFOrienter(const debruijn_graph::Graph& g): ContigOrienter(g) {}

    virtual pair<EdgeId, EdgeId> Apply(EdgeId e1, EdgeId e2) const {
        return make_pair(g_.conjugate(e1), e2);
    }
};

struct RROrienter: public ContigOrienter {
    RROrienter(const debruijn_graph::Graph& g): ContigOrienter(g) {}

    virtual pair<EdgeId, EdgeId> Apply(EdgeId e1, EdgeId e2) const {
        return make_pair(g_.conjugate(e1), g_.conjugate(e2));
    }
};



inline std::unique_ptr<ContigOrienter> GetContigOrienter(const debruijn_graph::Graph& g, io::LibraryOrientation orientation) {
    ContigOrienter* result;
    switch (orientation) {
        case io::LibraryOrientation::FF:  {
            result = new FFOrienter(g);
            break;
        }
        case io::LibraryOrientation::RR:  {
            result = new RROrienter(g);
            break;
        }
        case io::LibraryOrientation::FR:  {
            result = new FROrienter(g);
            break;
        }
        case io::LibraryOrientation::RF:  {
            result = new RFOrienter(g);
            break;
        }
        default: {
            result = new FFOrienter(g);
            break;
        }
    }
    return std::unique_ptr<ContigOrienter>(result);
}


};

#endif //PROJECT_BWA_PAIR_INFO_FILLER_HPP_H
