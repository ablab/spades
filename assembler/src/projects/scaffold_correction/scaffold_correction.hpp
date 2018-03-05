//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "io/reads/osequencestream.hpp"
#include "io/reads/single_read.hpp"
#include "io/reads/file_reader.hpp"
#include "pipeline/stage.hpp"
#include "pipeline/graph_pack.hpp"
#include "assembly_graph/paths/path_processor.hpp"
#include "stages/construction.hpp"
#include "pipeline/config_struct.hpp"
#include "assembly_graph/dijkstra/dijkstra_algorithm.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"

namespace scaffold_correction {
    typedef debruijn_graph::ConjugateDeBruijnGraph Graph;

    class PathSideSimilarityComparator {
    private:
        typedef Graph::EdgeId EdgeId;
        const Graph &graph_;
        vector<EdgeId> initial_path_;

    public:
        size_t CalculateDifference(vector<EdgeId> path) {
            size_t lside = 0;
            size_t rside = 0;
            size_t min_size = std::min(path.size(), initial_path_.size());
            while(lside < min_size && path[lside] == initial_path_[lside])
                lside++;
            while(rside < min_size && lside + rside < min_size &&
                    path[path.size() - 1 - rside] == initial_path_[initial_path_.size() - 1 - rside])
                rside++;
            size_t result = 0;
            for(size_t i = lside; i < path.size() - rside; i++) {
                result += graph_.length(path[i]);
            }
            return result;
        }

        PathSideSimilarityComparator(const Graph &graph, vector<EdgeId> path) : graph_(graph), initial_path_(path) {
        }


        bool operator()(const vector<EdgeId> &path1, const vector<EdgeId> &path2) {
            return CalculateDifference(path1) < CalculateDifference(path2);
        }
    };

    class CarefulPathFixer {
    private:
        typedef Graph::EdgeId EdgeId;
        typedef Graph::VertexId VertexId;
        const Graph &graph_;
        size_t max_cut_length_;
        size_t max_insert_;

        bool Consistent(EdgeId e1, EdgeId e2) const {
            return graph_.EdgeEnd(e1) == graph_.EdgeStart(e2);
        }

        size_t StepBack(size_t pos, const vector<EdgeId> &edges) const {
            size_t step_size = 0;
            while(pos > 0 && Consistent(edges[pos - 1], edges[pos]) && step_size + graph_.length(edges[pos]) <= max_cut_length_) {
                step_size += graph_.length(edges[pos]);
                pos -= 1;
            }
            return pos;
        }

        size_t StepForward(size_t pos, const vector<EdgeId> &edges) const {
            size_t step_size = 0;
            while(pos  + 1 < edges.size() && Consistent(edges[pos], edges[pos + 1])&& step_size + graph_.length(edges[pos]) <= max_cut_length_) {
                step_size += graph_.length(edges[pos]);
                pos += 1;
            }
            return pos;
        }


        void PrintPath(vector<EdgeId> path) const {
            for(size_t i = 0; i < path.size(); i++) {
                TRACE(graph_.EdgeNucls(path[i]));
            }
        }


        vector<EdgeId> TryCloseGap(VertexId v1, VertexId v2, const vector<EdgeId> &path) const {
            if (v1 == v2)
                return vector<EdgeId>();
            TRACE(
                    "Trying to close gap between v1=" << graph_.int_id(v1) << " and v2=" << graph_.int_id(v2));
            typedef omnigraph::DijkstraHelper<Graph>::PathIgnoringDijkstraSettings DS;
            size_t max_path_length = max_insert_ + 2 * max_cut_length_;
            DS ds(DS::LC(graph_, path), DS::VPrC(max_path_length), DS::VPuC(max_path_length), DS::NIF(graph_));
            omnigraph::Dijkstra<Graph, DS> dj(graph_, ds);
            dj.Run(v1);
            if(dj.DistanceCounted(v2) && dj.GetDistance(v2) <= max_insert_) {
                vector<EdgeId> result = dj.GetShortestPathTo(v2);
                VERIFY(graph_.EdgeStart(result.front()) == v1);
                VERIFY(graph_.EdgeEnd(result.back()) == v2);
                TRACE("Gap closed");
                TRACE("Cumulative closure length is " << CumulativeLength(graph_, result));
                TRACE("Difference from initial path: " << dj.GetDistance(v2));
                return result;
            } else {
                TRACE("Failed to find closing path");
                return vector<EdgeId>();
            }
/*
            PathSideSimilarityComparator comparator(garph_, path);
            omnigraph::BestPathStorage<Graph, PathSideSimilarityComparator> storage(graph_, comparator);
//            omnigraph::PathStorageCallback<Graph> path_store(graph_);
            //todo reduce value after investigation
            omnigraph::PathProcessor<Graph> path_processor(graph_, 0, max_insert_, v1, v2, 1000000, storage);
            path_processor.Process();
            TRACE(graph_.VertexNucls(v1));
            TRACE(graph_.VertexNucls(v2));
            size_t error_code = path_processor.Process();
            TRACE("Error code: " << error_code);

            if (storage.size() == 0) {
                TRACE("Failed to find closing path");
                return vector<EdgeId>();
            } else if (storage.size() == 1) {
                TRACE("Unique closing path found");
            } else {
                TRACE("Several closing paths found(" << path_store.paths().size() << "), first chosen");
            }
            auto tmp = path_store.paths();
            TRACE("Number of paths: " << tmp.size());
//            for(auto it = tmp.begin(); it != tmp.end(); ++it) {
//                TRACE(ConstructSequence(*it));
//            }
            vector<EdgeId> answer = storage.BestPath();
            TRACE("Gap closed");
            TRACE( "Cumulative closure length is " << CumulativeLength(graph_, answer));
            TRACE( "Difference from initial path: " <<  comparator.CalculateDifference(answer));
            return answer;
*/
        }

    public:
        CarefulPathFixer(const Graph &graph, size_t max_cut_length, size_t max_insert)
                : graph_(graph), max_cut_length_(max_cut_length), max_insert_(max_insert) {
        }

        vector<EdgeId> TryFixPath(const vector<EdgeId>& edges) const {
            vector<EdgeId> result;
            if (edges.empty()) {
                return vector<EdgeId>();
            }
            result.push_back(edges[0]);
            for (size_t i = 1; i < edges.size(); ++i) {
                if (!Consistent(result.back(), edges[i])) {
                    size_t lindex = StepBack(result.size() - 1, result);
                    size_t rindex = StepForward(i, edges);
                    vector<EdgeId> current_path(result.begin() + lindex + 1, result.end());
                    current_path.insert(current_path.end(), edges.begin() + i, edges.begin() + rindex);
                    vector<EdgeId> closure = TryCloseGap(graph_.EdgeEnd(result[lindex]), graph_.EdgeStart(edges[rindex]), current_path);
                    if(closure.size() != 0 || Consistent(result[lindex], edges[rindex])) {
                        result.resize(lindex + 1);
                        VERIFY(closure.size() == 0 || Consistent(result.back(), closure.front()));
                        result.insert(result.end(), closure.begin(), closure.end());
                        i = rindex;
                        VERIFY(Consistent(result.back(), edges[i]));
                    }
                }
                result.push_back(edges[i]);
            }
            return result;
        }
        DECL_LOGGER("CarefulPathFixer")
    };

    class ScaffoldCorrector {
        typedef debruijn_graph::conj_graph_pack graph_pack;
    private:

        const graph_pack& gp_;
        const CarefulPathFixer &fixer_;


        bool CheckPath(const vector<Graph::EdgeId> &path) const {
            if (!path.size())
                return false;

            for (size_t i = 1; i < path.size(); i++) {
                if (gp_.g.EdgeEnd(path[i - 1]) != gp_.g.EdgeStart(path[i]))
                    return false;
            }
            return true;
        }

        Sequence ConstructSequence(const vector<Graph::EdgeId> &path) const {
            Sequence result = gp_.g.EdgeNucls(path[0]);
            for(size_t i = 1; i < path.size(); i++) {
                result = result +  gp_.g.EdgeNucls(path[i]).Subseq(gp_.k_value);
            }
            return result;
        }

    public:
        ScaffoldCorrector(const graph_pack &gp, const CarefulPathFixer &fixer)
                : gp_(gp), fixer_(fixer) {}

        Sequence correct(const vector<Sequence> &scaffold) const {
            auto mapper = debruijn_graph::MapperInstance(gp_);
            MappingPath <debruijn_graph::EdgeId> path;
            for (auto it = scaffold.begin(); it != scaffold.end(); ++it) {
                path.join(mapper->MapSequence(*it));
            }
            vector<Graph::EdgeId> corrected_path = fixer_.TryFixPath(path.simple_path());
            if (CheckPath(corrected_path))
                return ConstructSequence(corrected_path);

            return Sequence();
        }
    };
}

namespace spades {
    class ScaffoldCorrectionStage : public AssemblyStage {
    public:
        typedef debruijn_graph::config::debruijn_config::scaffold_correction Config;
    private:
        std::string output_file_;
        const Config &config_;
    public:
        ScaffoldCorrectionStage(string output_file,
                const Config &config) :
                AssemblyStage("ScaffoldCorrection", "scaffold_correction"),
                output_file_(output_file), config_(config) {
        }

        vector<Sequence> CollectScaffoldParts(const io::SingleRead &scaffold) const {
            vector<Sequence> result;
            for (size_t i = 0; i < scaffold.size(); i++) {
                size_t j = i;
                while (j < scaffold.size() && is_nucl(scaffold.GetSequenceString()[j])) {
                    j++;
                }
                if (j > i) {
                    result.push_back(scaffold.Substr(i, j).sequence());
                    i = j - 1;
                }
            }
            return result;
        }

        void OutputResults(const vector<io::SingleRead> &results) {
            io::OFastaReadStream oss(output_file_);
            for (size_t i = 0; i < results.size(); i++) {
                string sequence = results[i].GetSequenceString();
                if (sequence != "") {
                    oss << io::SingleRead(results[i].name(), sequence);
                }
            }
        }

        vector<io::SingleRead> ReadScaffolds(const string &scaffolds_file) {
            io::FileReadStream scaffold_stream(scaffolds_file);
            vector<io::SingleRead> scaffolds;
            while (!scaffold_stream.eof()) {
                io::SingleRead scaffold;
                scaffold_stream >> scaffold;
                scaffolds.push_back(scaffold);
            }
            return scaffolds;
        }

        vector<io::SingleRead> RunParallelCorrection(const vector<io::SingleRead> &scaffolds, const scaffold_correction::ScaffoldCorrector &corrector) {
            vector<io::SingleRead> results(scaffolds.size());
#pragma omp parallel for
            for(size_t i = 0; i < scaffolds.size(); i++) {
                auto scaffold = scaffolds[i];
                std::string name = scaffold.name();
                vector<Sequence> part_list = CollectScaffoldParts(scaffold);
                TRACE("Correcting scaffold " << name);
                TRACE("Number of parts: " << part_list.size());
                Sequence result = corrector.correct(part_list);
                if (result.size() != 0) {
                    TRACE("Correction successful");
                    results[i] = io::SingleRead(name, result.str());
                } else if (config_.output_unfilled) {
                    TRACE("Correction unsuccessful. Using uncorrected scaffold");
                    results[i] = scaffold;
                }
            }
            return results;
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
            INFO("Correcting scaffolds from " << config_.scaffolds_file);
            scaffold_correction::CarefulPathFixer fixer(graph_pack.g, config_.max_cut_length, config_.max_insert);
            scaffold_correction::ScaffoldCorrector corrector(graph_pack, fixer);
            vector<io::SingleRead> scaffolds = ReadScaffolds(config_.scaffolds_file);
            graph_pack.EnsureIndex();
            vector<io::SingleRead> results = RunParallelCorrection(scaffolds, corrector);
            OutputResults(results);
            INFO(scaffolds.size() << " reads processed");
            INFO("Corrected scaffolds written to " << output_file_);
        }
        DECL_LOGGER("ScaffoldCorrectionStage")
    };

    void run_scaffold_correction() {
        INFO("Scaffold correction started");

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
        manager.add<debruijn_graph::Construction>()
               .add<ScaffoldCorrectionStage>(cfg::get().output_dir + "corrected_scaffolds.fasta", *cfg::get().sc_cor);
        INFO("Output directory: " << cfg::get().output_dir);
        conj_gp.kmer_mapper.Attach();
        manager.run(conj_gp, cfg::get().entry_point.c_str());
        INFO("Scaffold correction finished.");
    }

}
