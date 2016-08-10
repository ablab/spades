//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "utils/standard_base.hpp"
#include "pipeline/graph_pack.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "io/reads/io_helper.hpp"
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
    map<EdgeId, set<bin_id>> edge_annotation_;

    template<class BinCollection>
    void InnerStickAnnotation(EdgeId e, const BinCollection& bins) {
        edge_annotation_[e].insert(bins.begin(), bins.end());
    }

public:

    EdgeAnnotation(const conj_graph_pack& gp,
                   const set<bin_id>& bins_of_interest) :
                       gp_(gp),
                       bins_of_interest_(bins_of_interest)
    {
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

    vector<bin_id> Annotation(EdgeId e) const {
        if (!edge_annotation_.count(e)) {
            return {};
        }
        const auto& annotation = get(edge_annotation_, e);
        return vector<bin_id>(annotation.begin(), annotation.end());
    }

    set<bin_id> RelevantBins(const vector<EdgeId>& path) const {
        set<bin_id> answer;
        for (EdgeId e : path) {
            insert_all(answer, Annotation(e));
        }
        return answer;
    }

    set<EdgeId> EdgesOfBin(bin_id bin, size_t min_length = 0) const {
        set<EdgeId> answer;
        for (auto ann_pair : edge_annotation_) {
            if (ann_pair.second.count(bin) &&
                    gp_.g.length(ann_pair.first) > min_length) {
                answer.insert(ann_pair.first);
            }
        }
        return answer;
    }

    size_t size() const {
        return edge_annotation_.size();
    }

    const set<bin_id>& interesting_bins() const {
        return bins_of_interest_;
    }

};

class AnnotationFiller {
    const conj_graph_pack& gp_;
    set<bin_id> interesting_bins_;
    shared_ptr<SequenceMapper<Graph>> mapper_;

    vector<EdgeId> EdgesOfContig(const io::SingleRead& contig) const {
        return mapper_->MapRead(contig).simple_path();
    }

    Bins FilterInteresting(const Bins& bins) const {
        if (interesting_bins_.empty()) {
            return bins;
        } else {
            Bins answer;
            for (const bin_id& bin : bins) {
                if (interesting_bins_.count(bin)) {
                    answer.push_back(bin);
                } 
            }
            return answer;
        }
    }

    map<contig_id, std::set<bin_id>> LoadAnnotation(AnnotationStream& splits_annotation_stream) const {
        map<contig_id, std::set<bin_id>> annotation_map;
        INFO("Reading (split) contigs annotation");
        ContigAnnotation contig_annotation;
        size_t cnt = 0;
        while (!splits_annotation_stream.eof()) {
            splits_annotation_stream >> contig_annotation;
            auto bins = FilterInteresting(contig_annotation.second);
            if (!bins.empty()) {
                insert_all(annotation_map[contig_annotation.first], bins);
            }
            ++cnt;
        }
        INFO(cnt << " records read; annotation available for " << annotation_map.size() << " splits");
        return annotation_map;
    };

    void ProcessSplit(const io::SingleRead& split, std::set<bin_id> bins,
                      map<EdgeId, map<bin_id, size_t>>& coloring) const {
        auto mapping_path = mapper_->MapRead(split);
        for (size_t i = 0; i < mapping_path.size(); ++i) {
            auto map_info = mapping_path[i];
            MappingRange mr = map_info.second;
            auto& bin_lens = coloring[map_info.first];
            for (bin_id b : bins) {
                bin_lens[b] += mr.mapped_range.size();
            }
        }
    }

    map<EdgeId, map<bin_id, size_t>> FillColorInfo(io::SingleStream& splits_stream,
                                               const map<contig_id, std::set<bin_id>>& split_annotation) const {
        INFO("Sticking annotation to edges");
        map<EdgeId, map<bin_id, size_t>> answer;
        io::SingleRead split;
        while (!splits_stream.eof()) {
            splits_stream >> split;
            auto id = GetId(split);
            auto bins = split_annotation.find(id);
            if (bins != split_annotation.end() && !(bins->second.empty())) {
                ProcessSplit(split, bins->second, answer);
                //TODO think if it is overkill
                ProcessSplit(!split, bins->second, answer);
            }
        }
        INFO("Color info available for " << answer.size() << " edges");
        return answer;
    };

    void FilterSpuriousInfo(map<EdgeId, map<bin_id, size_t>>& coloring) const {
        for (auto& edge_info : coloring) {
            size_t edge_len = gp_.g.length(edge_info.first);
            for (auto color_it = edge_info.second.begin(); color_it != edge_info.second.end(); ) {
                if (math::ls(double(color_it->second) / double(edge_len), 0.3)) {
                    edge_info.second.erase(color_it++);
                } else {
                    ++color_it;
                }
            }
        }
    }

    set<bin_id> GatherAllBins(const map<EdgeId, map<bin_id, size_t>>& coloring) const {
        set<bin_id> answer;
        for (const auto& edge_info : coloring) {
            for (const auto& bin_info : edge_info.second) {
                answer.insert(bin_info.first);
            }
        }
        return answer;
    }

    set<bin_id> DetermineBins(const vector<EdgeId>& path,
                              const map<EdgeId, map<bin_id, size_t>>& coloring) const {
        map<bin_id, size_t> path_colors;
        size_t total_len = 0;
        for (EdgeId e : path) {
            size_t edge_len = gp_.g.length(e);
            total_len += edge_len;
            auto it = coloring.find(e);
            if (it != coloring.end()) {
                for (auto color_info : it->second) {
                    //TODO think carefully
                    path_colors[color_info.first] += edge_len; //color_info.second;
                }
            }
        }
        set<bin_id> answer;
        for (auto color_info : path_colors) {
            if (math::gr(double(color_info.second) / double(total_len), 0.3)) {
                answer.insert(color_info.first);
            }
        }
        return answer;
    }

public:

    AnnotationFiller(const conj_graph_pack& gp,
                     const vector<bin_id>& interesting_bins) :
        gp_(gp),
        interesting_bins_(interesting_bins.begin(), interesting_bins.end()),
        mapper_(MapperInstance(gp)) {
    }

    EdgeAnnotation operator() (io::SingleStream& contig_stream,
                     io::SingleStream& splits_stream,
                     AnnotationStream& splits_annotation_stream) {
        INFO("Filling edge annotation");
        INFO("Interesting bins " << interesting_bins_);

        auto coloring = FillColorInfo(splits_stream, LoadAnnotation(splits_annotation_stream));
        FilterSpuriousInfo(coloring);

        EdgeAnnotation edge_annotation(gp_, interesting_bins_.empty() ? GatherAllBins(coloring) : interesting_bins_);

        io::SingleRead contig;
        while (!contig_stream.eof()) {
            contig_stream >> contig;
            auto path = mapper_->MapRead(contig).simple_path();
            auto bins = DetermineBins(path, coloring);
            for (EdgeId e : path) {
                edge_annotation.StickAnnotation(e, bins);
            }
        }

        INFO("Edge annotation filled. Annotated " << edge_annotation.size() << " edges.");
        return edge_annotation;
    }
};
}
