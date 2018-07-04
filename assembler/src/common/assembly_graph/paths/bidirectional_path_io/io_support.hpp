//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/graph_support/contig_output.hpp"
#include "assembly_graph/components/connected_component.hpp"

namespace path_extend {
using namespace debruijn_graph;

struct ScaffoldInfo {
    std::string sequence;
    BidirectionalPath* path;
    std::string name;

    ScaffoldInfo(const std::string& sequence, BidirectionalPath* path) :
        sequence(sequence), path(path) { }

    size_t length() const {
        return sequence.length();
    }

    double coverage() const {
        return path->Coverage();
    }
};

typedef vector<ScaffoldInfo> ScaffoldStorage;

class ScaffoldSequenceMaker {
    const Graph &g_;
    const size_t k_;
public:
    ScaffoldSequenceMaker(const Graph& g) : g_(g), k_(g_.k()) {
    }

    string MakeSequence(const BidirectionalPath &scaffold) const;
};

//Finds common long edges in paths and joins them into
//Based on disjoint set union
class TranscriptToGeneJoiner {
private:
    const Graph &g_;
    size_t min_edge_len_; //minimal length for joining transcripts into a gene

    map<BidirectionalPath *, size_t, PathComparator> path_id_; //path ids
    std::vector<size_t> parents_; //node parents in
    std::vector<size_t> ranks_; //tree depth


    void MakeSet(size_t x);

    void JoinTrees(size_t x, size_t y);

    void Init(const PathContainer &paths);
public:
    TranscriptToGeneJoiner(const Graph &g, size_t min_edge_len): g_(g), min_edge_len_(min_edge_len) {}

    size_t FindTree(size_t x);

    size_t GetPathId(BidirectionalPath *path);

    void Construct(const PathContainer &paths);
};

class ContigNameGenerator {
public:
    virtual void Preprocess(const PathContainer& paths) = 0;

    virtual std::string MakeContigName(size_t index, const ScaffoldInfo &scaffold_info) = 0;

    virtual ~ContigNameGenerator() {
    }
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

    unordered_map<size_t, size_t> isoform_num_;
    unordered_map<size_t, size_t> gene_ids_;
    size_t gene_num_;

public:
    TranscriptNameGenerator(const Graph &g, size_t min_edge_len = 300):
        transcript_joiner_(g, min_edge_len),
        isoform_num_(),
        gene_ids_(),
        gene_num_(0) {

    }

    void Preprocess(const PathContainer& paths) override {
        transcript_joiner_.Construct(paths);
    }

    std::string MakeContigName(size_t index, const ScaffoldInfo &scaffold_info) override {
        size_t id = transcript_joiner_.GetPathId(scaffold_info.path);
        size_t parent_id = transcript_joiner_.FindTree(id);
        DEBUG("Path " << id << " Parent " << parent_id);
        if (gene_ids_.find(parent_id) == gene_ids_.end()) {
            gene_ids_[parent_id] = gene_num_;
            isoform_num_[parent_id] = 0;
            gene_num_++;
        }
        string contig_id = io::MakeRNAContigId(index,
                                               scaffold_info.length(),
                                               scaffold_info.coverage(),
                                               gene_ids_[parent_id],
                                               isoform_num_[parent_id]);
        isoform_num_[parent_id]++;
        return contig_id;
    }
};


inline std::shared_ptr<ContigNameGenerator> MakeContigNameGenerator(config::pipeline_type mode,
                                                                    const conj_graph_pack &gp) {
    std::shared_ptr<path_extend::ContigNameGenerator> name_generator;
    if (mode == config::pipeline_type::plasmid)
        name_generator = make_shared<PlasmidContigNameGenerator>(gp.components);
    else if (mode == config::pipeline_type::rna)
        name_generator = make_shared<TranscriptNameGenerator>(gp.g);
    else
        name_generator = make_shared<DefaultContigNameGenerator>();
    return name_generator;
}

class ScaffoldBreaker {
private:

    int min_overlap_;

    void SplitPath(const BidirectionalPath& path, PathContainer &result) const;

public:

    ScaffoldBreaker(int min_overlap): min_overlap_(min_overlap) {}

    void Break(const PathContainer &paths, PathContainer &result) const;
};

}
