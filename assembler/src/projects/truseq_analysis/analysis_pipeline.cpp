//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//
// Created by anton on 16.05.15.
//

#include "io/reads/file_reader.hpp"
#include "stages/construction.hpp"
#include "utils/standard_base.hpp"
#include "analysis_pipeline.hpp"

spades::VariationDetectionStage::VariationDetectionStage(string output_file, const Config &config) : AssemblyStage("VariationDetection", "variation_detection"),
                                                                                                     output_file_(output_file), config_(config) {
}

vector<io::SingleRead> spades::VariationDetectionStage::ReadScaffolds(const string &scaffolds_file) {
    io::FileReadStream scaffold_stream(scaffolds_file);
    vector<io::SingleRead> scaffolds;
    while(!scaffold_stream.eof()) {
        io::SingleRead scaffold;
        scaffold_stream >> scaffold;
        scaffolds.push_back(scaffold);
    }
    return scaffolds;
}

void spades::VariationDetectionStage::run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
    using debruijn_graph::EdgeId;
    using alignment_analysis::EdgeRange;
    INFO("Analysis of contigs from " << config_.scaffolds_file);
    vector<io::SingleRead> scaffolds = ReadScaffolds(config_.scaffolds_file);
    vector<io::SingleRead> genome = ReadScaffolds(config_.genome_file);
    auto mapper_ptr = MapperInstance(graph_pack);
//            alignment_analysis::AlignmentAnalyser aa(scaffolds, genome, graph_pack.g, *mapper_ptr);
    const debruijn_graph::DeBruijnGraph &graph = graph_pack.g;
    alignment_analysis::AlignmentAnalyserNew aa(graph, 2 * graph_pack.k_value + 10);

    ofstream os(output_file_);
    for(auto it = genome.begin(); it != genome.end(); ++it) {
        const io::SingleRead &part = *it;
        MappingPath<EdgeId> path = mapper_ptr->MapRead(part);
        vector<alignment_analysis::ConsistentMapping> result = aa.Analyse(path);
        os << "Analysis of part " << part.name() << endl;
        for(size_t i = 0; i < result.size(); i++) {
            alignment_analysis::ConsistentMapping &cm = result[i];
//            os << "Alignment: " << cm.GetInitialRange() << " -> ";
//            const vector<EdgeRange> &mappedPath = cm.GetMappedPath();
//            for(auto pit = mappedPath.begin(); pit != mappedPath.end(); ++pit) {
//                const EdgeRange &er = *pit;
//                os << er << " ";
//            }
//            os << endl;
            size_t diff = cm.GetInitialRange().size() > cm.size() ? cm.GetInitialRange().size() - cm.size() : cm.size() - cm.GetInitialRange().size();
            if(diff > 500)
                os << cm.CompareToReference(part.GetSequenceString()) << endl;
        }
        result = ExtractConsistentMappings(result);
        for(size_t i = 0; i + 1 < result.size(); i++) {
            alignment_analysis::ConsistentMapping &cm = result[i];
            alignment_analysis::ConsistentMapping &next_cm = result[i + 1];
            if (this->CheckEndVertex(graph, cm.EndEdge(), 150 + cm.Back().second.end_pos) &&
                this->CheckEndVertex(graph, graph.conjugate(next_cm.StartEdge()),
                                     150 + graph.length(next_cm.StartEdge()) - next_cm.Front().second.start_pos)) {
//                    os << "Coverage break: " << "[" << cm.GetInitialRange().end_pos << ", " << next_cm.GetInitialRange().start_pos << "]"<< endl;
            } else {
                if(cm.GetInitialRange().size() < 100 || next_cm.GetInitialRange().size() < 100) {
//                    os << "Unreliable alignment event: " << "[" << cm.GetInitialRange().end_pos << ", " <<
//                    next_cm.GetInitialRange().start_pos << "]" << endl;
                } else {
                    os << "Breakpoint: " << "[" << cm.GetInitialRange().end_pos << ", " <<
                    next_cm.GetInitialRange().start_pos << "]" << endl;
                }
            }
        }
    }
    os.close();
    INFO("Analisys results written to " << output_file_);
}

void spades::run_truseq_analysis() {
    INFO("TruSeq analysis started");

    debruijn_graph::conj_graph_pack conj_gp(cfg::get().K,
                                            cfg::get().tmp_dir,
                                            cfg::get().ds.reads.lib_count(),
                                            cfg::get().ds.reference_genome,
                                            cfg::get().flanking_range,
                                            cfg::get().pos.max_mapping_gap,
                                            cfg::get().pos.max_gap_diff);
    StageManager manager({cfg::get().developer_mode,
                          cfg::get().load_from,
                          cfg::get().output_saves});
    manager.add(new debruijn_graph::Construction());
    std::string output_file = cfg::get().output_dir + "analysis_report";
    manager.add(new VariationDetectionStage(output_file, cfg::get().tsa));
    INFO("Output directory: " << cfg::get().output_dir);
    conj_gp.kmer_mapper.Attach();
    manager.run(conj_gp, cfg::get().entry_point.c_str());
    INFO("Scaffold correction finished.");
}

bool spades::VariationDetectionStage::CheckEndVertex(debruijn_graph::DeBruijnGraph const &graph,
                                                     debruijn_graph::EdgeId e, size_t dist) {
    using debruijn_graph::VertexId;
    if(graph.coverage(e) == 0) {
        return true;
    }
    if (graph.length(e) > dist) {
        return false;
    }
    VertexId v = graph.EdgeEnd(e);
    GraphCore<debruijn_graph::DeBruijnDataMaster>::IteratorContainer outgoingEdges = graph.OutgoingEdges(v);
    for(auto it = outgoingEdges.begin(); it != outgoingEdges.end(); ++it) {
        if(!CheckEndVertex(graph, *it, dist - graph.length(e)))
            return false;
    }
    return true;
}

vector <alignment_analysis::ConsistentMapping> spades::VariationDetectionStage::ExtractConsistentMappings(const vector<alignment_analysis::ConsistentMapping> &path) {
    vector <alignment_analysis::ConsistentMapping> result;
    result.push_back(path[0]);
    for (size_t i = 1; i < path.size(); i++) {
        if (result.empty()) {
            result.push_back(path[i]);
        } else {
            alignment_analysis::ConsistentMapping &back = result.back();
            if (back.CheckConnect(path[i])) {
                back.Join(path[i]);
            } else {
                result.push_back(path[i]);
            }
        }
    }
    return result;
}
