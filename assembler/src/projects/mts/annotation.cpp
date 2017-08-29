//***************************************************************************
//* Copyright (c) 2015-2017 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#include "annotation.hpp"

namespace debruijn_graph {

//------------------------------------------------------------------------------

ContigAnnotation AnnotationStream::Parse(const std::string& s) const {
    ContigAnnotation annotation;
    std::istringstream ss(s);
    ss >> annotation.first;
    while (true) {
        bin_id bin;
        ss >> bin;
        if (ss.fail())
            break;
        annotation.second.push_back(bin);
    }
    return annotation;
}

AnnotationStream& AnnotationStream::operator >>(ContigAnnotation& annotation) {
    VERIFY(!inner_stream_.eof())

    annotation = Parse(line_);
    std::getline(inner_stream_, line_);
    return *this;
}

//------------------------------------------------------------------------------

AnnotationOutStream& AnnotationOutStream::operator <<(const ContigAnnotation& annotation) {
    inner_stream_ << annotation.first;
    string delim = "\t";
    for (bin_id bin : annotation.second) {
        inner_stream_ << delim << bin;
        delim = " ";
    }
    inner_stream_ << endl;
    return *this;
}

//------------------------------------------------------------------------------

Bins EdgeAnnotation::Annotation(EdgeId e) const {
    if (!edge_annotation_.count(e)) {
        return {};
    }
    const auto& annotation = utils::get(edge_annotation_, e);
    return vector<bin_id>(annotation.begin(), annotation.end());
}

BinSet EdgeAnnotation::RelevantBins(const vector<EdgeId>& path) const {
    BinSet answer;
    for (EdgeId e : path) {
        utils::insert_all(answer, Annotation(e));
    }
    return answer;
}

set<EdgeId> EdgeAnnotation::EdgesOfBin(bin_id bin, size_t min_length) const {
    set<EdgeId> answer;
    for (auto ann_pair : edge_annotation_) {
        if (ann_pair.second.count(bin) &&
                gp_.g.length(ann_pair.first) > min_length) {
            answer.insert(ann_pair.first);
        }
    }
    return answer;
}

//------------------------------------------------------------------------------

template<typename K, typename V>
std::ostream& operator<<(std::ostream& str, const std::map<K, V>& map) {
    str << "{";
    for (const auto& kv : map)
        str << kv.first << ": " << kv.second << ", ";
    str << "}";
    return str;
}

vector<EdgeId> AnnotationFiller::EdgesOfContig(const io::SingleRead& contig) const {
    return mapper_->MapRead(contig).simple_path();
}

Bins AnnotationFiller::FilterInteresting(const Bins& bins) const {
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

AnnotationFiller::AnnotationMap AnnotationFiller::LoadAnnotation(AnnotationStream& splits_annotation_stream) const {
    AnnotationFiller::AnnotationMap annotation_map;
    INFO("Reading (split) contigs annotation");
    ContigAnnotation contig_annotation;
    size_t cnt = 0;
    while (!splits_annotation_stream.eof()) {
        splits_annotation_stream >> contig_annotation;
        auto bins = FilterInteresting(contig_annotation.second);
        if (!bins.empty()) {
            utils::insert_all(annotation_map[contig_annotation.first], bins);
        }
        ++cnt;
    }
    INFO(cnt << " records read; annotation available for " << annotation_map.size() << " splits");
    return annotation_map;
}

void AnnotationFiller::ProcessSplit(const io::SingleRead& split, std::set<bin_id> bins,
                  ColoringMap& coloring) const {
    auto mapping_path = mapper_->MapRead(split);
    for (const auto& map_info : mapping_path) {
        MappingRange mr = map_info.second;
        auto& bin_lens = coloring[map_info.first];
        for (bin_id b : bins) {
            bin_lens[b] += mr.mapped_range.size();
        }
        TRACE("Initially " << map_info.first << " is colored into " << bin_lens);
    }
}

AnnotationFiller::ColoringMap AnnotationFiller::FillColorInfo(io::SingleStream& splits_stream,
                                               const AnnotationMap& split_annotation) const {
    INFO("Sticking annotation to edges");
    AnnotationFiller::ColoringMap answer;
    io::SingleRead split;
    while (!splits_stream.eof()) {
        splits_stream >> split;
        auto id = GetId(split);
        auto bins = split_annotation.find(id);
        if (bins != split_annotation.end() && !(bins->second.empty())) {
            DEBUG("Split " << id << " is colored into " << bins->second);
            ProcessSplit(split, bins->second, answer);
            //TODO think if it is overkill
            ProcessSplit(!split, bins->second, answer);
        }
    }
    INFO("Color info available for " << answer.size() << " edges");
    return answer;
}

bool AnnotationFiller::IsSpurious(size_t colored_len, size_t full_len) {
    return math::ls(double(colored_len) / double(full_len), 0.3); //FIXME: extract magic constant to config
}

void AnnotationFiller::FilterSpuriousInfo(ColoringMap& coloring) const {
    for (auto& edge_info : coloring) {
        size_t edge_len = gp_.g.length(edge_info.first);
        for (auto color_it = edge_info.second.begin(); color_it != edge_info.second.end(); ) {
            if (IsSpurious(color_it->second, edge_len)) {
                edge_info.second.erase(color_it++);
            } else {
                ++color_it;
            }
        }
    }
}

BinSet AnnotationFiller::GatherAllBins(const ColoringMap& coloring) const {
    set<bin_id> answer;
    for (const auto& edge_info : coloring) {
        for (const auto& bin_info : edge_info.second) {
            answer.insert(bin_info.first);
        }
    }
    return answer;
}

BinSet AnnotationFiller::DetermineBins(const vector<EdgeId>& path,
                                       const ColoringMap& coloring) const {
    ColoringLengths path_colors;
    size_t total_len = 0;
    for (const auto& e : path) {
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
    TRACE("Total len: " << total_len << "; candidates: " << path_colors);
    BinSet answer;
    using ColorInfo = ColoringLengths::value_type;
    auto it = std::max_element(path_colors.begin(), path_colors.end(),
        [](const ColorInfo& p1, const ColorInfo& p2) {
             return p1.second < p2.second;
        });
    //Majority strategy: choose the longest coloring
    //if (it != path_colors.end() && !IsSpurious(it->second, total_len))
    //    answer.insert(it->first);
    for (auto color_info : path_colors) {
        if (!IsSpurious(color_info.second, total_len)) {
            answer.insert(color_info.first);
        }
    }
    return answer;
}

EdgeAnnotation AnnotationFiller::operator() (io::SingleStream& contig_stream,
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
        DEBUG("Filling annotation for contig " << contig.name());
        auto raw_path = mapper_->MapRead(contig);
        std::vector<EdgeId> path;
        path.reserve(raw_path.size());
        for (const auto& ep : raw_path) { //Filter the poorly mapped edges
            EdgeId e = ep.first;
            size_t edge_len = gp_.g.length(e);
            if (math::ge(double(ep.second.mapped_range.size()) / double(edge_len), 0.9)) //FIXME: extract magic constant to config
                path.push_back(e);
        }
        auto bins = DetermineBins(path, coloring);
        for (const auto& e : path) {
            TRACE("Edge " << e.int_id() << " will be colored to " << bins);
            edge_annotation.StickAnnotation(e, bins);
        }
    }

    INFO("Edge annotation filled. Annotated " << edge_annotation.size() << " edges.");
    return edge_annotation;
}

}
