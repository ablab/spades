//***************************************************************************
//* Copyright (c) 2015-2017 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "pipeline/graph_pack.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "io/reads/io_helper.hpp"
#include "formats.hpp"

namespace debruijn_graph {

class AnnotationStream {
    std::ifstream inner_stream_;
    std::string line_;

    ContigAnnotation Parse(const std::string& s) const;
public:

    AnnotationStream(const std::string& fn) : inner_stream_(fn) {
        std::getline(inner_stream_, line_);
    }

    bool eof() const {
        return inner_stream_.eof();
    }

    AnnotationStream& operator >>(ContigAnnotation& annotation);

    void close() {
        inner_stream_.close();
    }
};

class AnnotationOutStream {
    std::ofstream inner_stream_;
public:

    AnnotationOutStream(const std::string& fn) : inner_stream_(fn) {
    }

    AnnotationOutStream& operator <<(const ContigAnnotation& annotation);

    void close() {
        inner_stream_.close();
    }
};

class EdgeAnnotation {
    const GraphPack& gp_;
    BinSet bins_of_interest_;
    std::map<EdgeId, BinSet> edge_annotation_;

    template<class BinCollection>
    void InnerStickAnnotation(EdgeId e, const BinCollection& bins) {
        edge_annotation_[e].insert(bins.begin(), bins.end());
    }

public:

    EdgeAnnotation(const GraphPack& gp,
                   const std::set<bin_id>& bins_of_interest) :
                       gp_(gp),
                       bins_of_interest_(bins_of_interest)
    {
    }

    template<class BinCollection>
    void StickAnnotation(EdgeId e, const BinCollection& bins) {
        InnerStickAnnotation(e, bins);
        InnerStickAnnotation(gp_.get<Graph>().conjugate(e), bins);
    }

    void StickAnnotation(EdgeId e, const bin_id& bin) {
        StickAnnotation(e, std::vector<bin_id>{bin});
    }

    template<class EdgeCollection>
    void StickAnnotation(const EdgeCollection& edges, const bin_id& bin) {
        for (EdgeId e : edges) {
            StickAnnotation(e, bin);
        }
    }

    std::vector<bin_id> Annotation(EdgeId e) const;
    std::set<bin_id> RelevantBins(const std::vector<EdgeId> &path) const;
    std::set<EdgeId> EdgesOfBin(bin_id bin, size_t min_length = 0) const;

    size_t size() const {
        return edge_annotation_.size();
    }

    const std::set<bin_id> &interesting_bins() const {
        return bins_of_interest_;
    }
};

class AnnotationFiller {

    const GraphPack& gp_;
    BinSet interesting_bins_;
    std::shared_ptr<SequenceMapper<Graph>> mapper_;

    std::vector<EdgeId> EdgesOfContig(const io::SingleRead& contig) const;

    typedef std::map<contig_id, BinSet> AnnotationMap;
    typedef std::map<bin_id, size_t> ColoringLengths;
    typedef std::map<EdgeId, ColoringLengths> ColoringMap;

    Bins FilterInteresting(const Bins& bins) const;
    AnnotationMap LoadAnnotation(AnnotationStream& splits_annotation_stream) const;
    void ProcessSplit(const io::SingleRead& split, std::set<bin_id> bins,
                      ColoringMap& coloring) const;
    ColoringMap FillColorInfo(io::SingleStream& splits_stream,
                              const AnnotationMap& split_annotation) const;

    static bool IsSpurious(size_t colored_len, size_t full_len);

    void FilterSpuriousInfo(ColoringMap& coloring) const;
    std::set<bin_id> GatherAllBins(const ColoringMap &coloring) const;
    std::set<bin_id> DetermineBins(const std::vector<EdgeId> &path, const ColoringMap &coloring) const;
public:

    AnnotationFiller(const GraphPack& gp, const std::vector<bin_id> &interesting_bins) :
        gp_(gp),
        interesting_bins_(interesting_bins.begin(), interesting_bins.end()),
        mapper_(MapperInstance(gp)) {
    }

    EdgeAnnotation operator() (io::SingleStream& contig_stream,
                     io::SingleStream& splits_stream,
                     AnnotationStream& splits_annotation_stream);

private:
    DECL_LOGGER("AnnotationFiller");
};

}
