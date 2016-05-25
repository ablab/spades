//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "dev_support/standard_base.hpp"
#include "pipeline/graph_pack.hpp"
#include "assembly_graph/graph_alignment/sequence_mapper.hpp"
#include "io/reads_io/io_helper.hpp"
#include "formats.hpp"

namespace debruijn_graph {

class AnnotationStream {
    std::ifstream inner_stream_;
    std::string line_;

    ContigAnnotation Parse(const std::string& s) const {
        ContigAnnotation annotation;
        stringstream ss(s);
        ss >> annotation.first;
        string delim;
        ss >> delim;
        VERIFY(delim == ":");
        while (true) {
            bin_id bin;
            ss >> bin;
            if (ss.fail())
                break;
            annotation.second.push_back(bin);
        }
        return annotation;
    }

public:

    AnnotationStream(const std::string& fn) : inner_stream_(fn) {
        std::getline(inner_stream_, line_);
    }

    bool eof() const {
        return inner_stream_.eof();
    }

    AnnotationStream& operator >>(ContigAnnotation& annotation) {
        VERIFY(!inner_stream_.eof())

        annotation = Parse(line_);
        std::getline(inner_stream_, line_);
        return *this;
    }

    void close() {
        inner_stream_.close();
    }
};

class AnnotationOutStream {
    std::ofstream inner_stream_;
public:

    AnnotationOutStream(const std::string& fn) : inner_stream_(fn) {
    }

    AnnotationOutStream& operator <<(const ContigAnnotation& annotation) {
        inner_stream_ << annotation.first;
        string delim = " : ";
        for (bin_id bin : annotation.second) {
            inner_stream_ << delim << bin;
            delim = " ";
        }
        inner_stream_ << endl;
        return *this;
    }

    void close() {
        inner_stream_.close();
    }
};

class EdgeAnnotation {
    const conj_graph_pack& gp_;
    set<bin_id> bins_of_interest_;
    shared_ptr<SequenceMapper<Graph>> mapper_;
    map<EdgeId, set<bin_id>> edge_annotation_;

    Bins FilterInteresting(const Bins& bins) const {
        Bins answer;
        for (const bin_id& bin : bins) {
            if (bins_of_interest_.count(bin)) {
                answer.push_back(bin);
            }
        }
        return answer;
    }

    template<class BinCollection>
    void InnerStickAnnotation(EdgeId e, const BinCollection& bins) {
        edge_annotation_[e].insert(bins.begin(), bins.end());
    }

public:

    EdgeAnnotation(const conj_graph_pack& gp,
                   const vector<bin_id>& bins_of_interest) :
                       gp_(gp),
                       bins_of_interest_(bins_of_interest.begin(), bins_of_interest.end()),
                       mapper_(MapperInstance(gp)) {
    }

    template<class BinCollection>
    void StickAnnotation(EdgeId e, const BinCollection& bins) {
        InnerStickAnnotation(e, bins);
        InnerStickAnnotation(gp_.g.conjugate(e), bins);
    }

    void StickAnnotation(EdgeId e, const bin_id& bin) {
        StickAnnotation(e, vector<bin_id>{bin});
    }

    template<class EdgeCollection>
    void StickAnnotation(const EdgeCollection& edges, const bin_id& bin) {
        for (EdgeId e : edges) {
            StickAnnotation(e, bin);
        }
    }

    void Fill(io::SingleStream& contigs, AnnotationStream& annotation_stream) {
        INFO("Filling edge annotation");
        set<bin_id> all_bins;

        INFO("Reading (split) contigs annotation");
        map<contig_id, std::set<bin_id>> annotation_map;
        ContigAnnotation contig_annotation;
        while (!annotation_stream.eof()) {
            annotation_stream >> contig_annotation;
            auto bins = contig_annotation.second;
            if (!bins_of_interest_.empty()) {
                bins = FilterInteresting(bins);
            } else {
                insert_all(all_bins, bins);
            }
            if (!bins.empty()) {
                insert_all(annotation_map[GetBaseId(contig_annotation.first)], bins);
            }
        }
        INFO("Annotation available for " << annotation_map.size() << " contigs");
        INFO("Sticking annotation to edges");

        io::SingleRead contig;
        while (!contigs.eof()) {
            contigs >> contig;
            contig_id id = GetBaseId(GetId(contig));
            auto bins = annotation_map.find(id);
            if (bins != annotation_map.end() && !(bins->second.empty())) {
                for (EdgeId e : mapper_->MapRead(contig).simple_path()) {
                    StickAnnotation(e, bins->second);
                }
            }
        }

        if (bins_of_interest_.empty()) {
            INFO("Bins of interest not specified. Marking all bins as bins of interest");
            bins_of_interest_ = all_bins;
        }
        INFO("Edge annotation filled. Annotated " << edge_annotation_.size() << " edges.");
    }

    vector<bin_id> Annotation(EdgeId e) const {
        if (!edge_annotation_.count(e)) {
            return {};
        }
        const auto& annotation = get(edge_annotation_, e);
        return vector<bin_id>(annotation.begin(), annotation.end());
    }

    set<bin_id> RelevantBins(const io::SingleRead& r) const {
        set<bin_id> answer;
        for (EdgeId e : EdgesOfContig(r)) {
            insert_all(answer, Annotation(e));
        }
        return answer;
    }

    vector<EdgeId> EdgesOfContig(const io::SingleRead& contig) const {
        //TODO: memoize mapping
        return mapper_->MapRead(contig).simple_path();
    }

    set<EdgeId> EdgesOfBin(bin_id bin) const {
        set<EdgeId> answer;
        for (auto ann_pair : edge_annotation_) {
            if (ann_pair.second.count(bin)) {
                answer.insert(ann_pair.first);
            }
        }
        return answer;
    }

    const set<bin_id>& interesting_bins() const {
        return bins_of_interest_;
    }

};

}
