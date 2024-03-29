//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/paths/bidirectional_path_container.hpp"
#include "assembly_graph/components/connected_component.hpp"
#include "io/reads/header_naming.hpp"

namespace path_extend {
using namespace debruijn_graph;

struct ScaffoldInfo {
    std::string sequence;
    const BidirectionalPath* path;
    std::string name;

    ScaffoldInfo(const std::string& sequence, const BidirectionalPath* path) :
        sequence(sequence), path(path) { }

    size_t length() const {
        return sequence.length();
    }

    double coverage() const {
        return path->Coverage();
    }
};

typedef std::vector<ScaffoldInfo> ScaffoldStorage;

class ScaffoldSequenceMaker {
    const Graph &g_;
    const size_t k_;
public:
    ScaffoldSequenceMaker(const Graph& g)
            : g_(g), k_(g_.k()) {}

    std::string MakeSequence(const BidirectionalPath &scaffold) const;
};

//Finds common long edges in paths and joins them into
//Based on disjoint set union
class TranscriptToGeneJoiner {
private:
    const Graph &g_;
    size_t min_edge_len_; //minimal length for joining transcripts into a gene

    std::unordered_map<uint64_t, size_t> path_id_; //path ids
    std::vector<size_t> parents_; //node parents in
    std::vector<size_t> ranks_; //tree depth

    void MakeSet(size_t x);

    void JoinTrees(size_t x, size_t y);

    void Init(const PathContainer &paths);

    DECL_LOGGER("TranscriptToGeneJoiner");
public:
    TranscriptToGeneJoiner(const Graph &g, size_t min_edge_len): g_(g), min_edge_len_(min_edge_len) {}

    size_t FindTree(size_t x);

    size_t GetPathId(const BidirectionalPath &path) const;

    void Construct(const PathContainer &paths);
};

class ContigNameGenerator {
public:
    virtual void Preprocess(const PathContainer& paths) = 0;

    virtual std::string MakeContigName(size_t index, const ScaffoldInfo &scaffold_info) = 0;

    virtual void PrintStats() const {}

    virtual ~ContigNameGenerator() {}
};

class DefaultContigNameGenerator: public ContigNameGenerator {
public:
    void Preprocess(const PathContainer&) override {}

    std::string MakeContigName(size_t index, const ScaffoldInfo &scaffold_info) override {
        return io::MakeContigId(index, scaffold_info.length(), scaffold_info.coverage());
    }
};

class PlasmidContigNameGenerator: public ContigNameGenerator {
    const ConnectedComponentCounter &c_counter_;

public:
    PlasmidContigNameGenerator(const ConnectedComponentCounter &c_counter): c_counter_(c_counter) {}

    void Preprocess(const PathContainer&) override {}

    std::string MakeContigName(size_t index, const ScaffoldInfo &scaffold_info) override {
        return io::AddComponentId(io::MakeContigId(index, scaffold_info.length(), scaffold_info.coverage()),
                                  c_counter_.GetComponent(scaffold_info.path->Front()));
    }
};

class TranscriptNameGenerator: public ContigNameGenerator {
    TranscriptToGeneJoiner transcript_joiner_;

    std::unordered_map<size_t, size_t> isoform_num_;
    std::unordered_map<size_t, size_t> gene_ids_;
    size_t gene_num_;

public:
    TranscriptNameGenerator(const Graph &g, size_t min_edge_len = 300):
        transcript_joiner_(g, min_edge_len),
        isoform_num_(),
        gene_ids_(),
        gene_num_(0) {}

    void Clear() {
        isoform_num_.clear();
        gene_ids_.clear();
        gene_num_ = 0;
    }

    void Preprocess(const PathContainer& paths) override {
        Clear();
        transcript_joiner_.Construct(paths);
    }

    std::string MakeContigName(size_t index, const ScaffoldInfo &scaffold_info) override {
        size_t id = transcript_joiner_.GetPathId(*scaffold_info.path);
        size_t parent_id = transcript_joiner_.FindTree(id);
        DEBUG("Path " << id << " Parent " << parent_id);
        if (gene_ids_.find(parent_id) == gene_ids_.end()) {
            gene_ids_[parent_id] = gene_num_;
            isoform_num_[parent_id] = 0;
            gene_num_++;
        }
        auto contig_id = io::MakeRNAContigId(index, scaffold_info.length(), scaffold_info.coverage(),
                                             gene_ids_[parent_id], isoform_num_[parent_id]);
        isoform_num_[parent_id]++;
        return contig_id;
    }

    size_t gene_count() const {
        return gene_num_;
    }

    size_t isoform_count() const {
        size_t isoform_cnt = 0;
        for (const auto& it: isoform_num_) {
            isoform_cnt += it.second;
        }
        return isoform_cnt;
    }

    void PrintStats() const override {
        INFO("Total number of genes: " << gene_count())
        INFO("Total number of isoforms: " << isoform_count())
    }
};

class ScaffoldBreaker {
private:
    int min_overlap_;
    void SplitPath(const BidirectionalPath &path, PathContainer &result) const;

public:
    ScaffoldBreaker(int min_overlap): min_overlap_(min_overlap) {}
    void Break(const PathContainer &paths, PathContainer &result) const;
};

}
