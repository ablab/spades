#include "molecule_extraction_stage.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "common/modules/path_extend/path_filter.hpp"
#include "common/paired_info/pair_info_filler.hpp"
#include "projects/spades/distance_estimation.hpp"
#include "common/modules/path_extend/paired_library.hpp"
#include "common/modules/path_extend/weight_counter.hpp"
#include "common/modules/path_extend/path_extender.hpp"
#include "common/modules/path_extend/pe_resolver.hpp"
#include "modules/path_extend/pipeline/launch_support.hpp"
#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"

namespace debruijn_graph {
    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> lib_t;
    typedef std::vector<io::SequencingLibrary<debruijn_graph::config::LibraryData>> lib_vector_t;

    bool Has10XRNAReads(config::dataset ds) {
        bool has_10X_RNA = false;
        for (const auto& lib: ds.reads) {
            if (lib.type() == io::LibraryType::RNA10x) {
                has_10X_RNA = true;
            }
        }
        return has_10X_RNA;
    }

    class SimplePathExtractor {
    private:
        const debruijn_graph::conj_graph_pack &gp_;
        std::set<EdgeId> used_edges_;
        PathStorage<Graph> temp_set_;
        std::vector<EdgeId> getForwardIntersect(EdgeId e, const std::unordered_set<EdgeId> &edge_set) const {
            auto edges = gp_.g.IncidentEdges(gp_.g.EdgeEnd(e));
            std::vector<EdgeId> filtered;
            for (auto temp_e : edges) {
                if (temp_e != e && edge_set.count(temp_e) && !used_edges_.count(temp_e)) {
                    filtered.push_back(temp_e);
                }
            }
            return filtered;
        }

        std::vector<EdgeId> getReverseIntersect(EdgeId e, const std::unordered_set<EdgeId> &edge_set) const {
            auto edges = gp_.g.IncidentEdges(gp_.g.EdgeStart(e));
            std::vector<EdgeId> filtered;
            for (auto temp_e : edges) {
                if (temp_e != e && edge_set.count(temp_e) && !used_edges_.count(temp_e)) {
                    filtered.push_back(temp_e);
                }
            }
            return filtered;
        }

        std::vector<EdgeId> getForwardIntersect(EdgeId e, const GraphComponent<Graph> &comp) const {

            auto edges = gp_.g.IncidentEdges(gp_.g.EdgeEnd(e));
            DEBUG("Edges " << edges);
            std::vector<EdgeId> filtered;
            for (auto temp_e : edges) {
                if (temp_e != e && comp.edges().count(temp_e)) {
                    filtered.push_back(temp_e);
                }
            }
            DEBUG("Filtered " << filtered);
            return filtered;
        }

        std::vector<EdgeId> getReverseIntersect(EdgeId e, const GraphComponent<Graph> &comp) const {
            auto edges = gp_.g.IncidentEdges(gp_.g.EdgeStart(e));
            std::vector<EdgeId> filtered;
            for (auto temp_e : edges) {
                if (temp_e != e && comp.edges().count(temp_e)) {
                    filtered.push_back(temp_e);
                }
            }
            return filtered;
        }

        void extendForward(const std::unordered_set<EdgeId> &edge_set,
                           std::deque<EdgeId> &linear_path, EdgeId e) {
            auto extensions = getForwardIntersect(e, edge_set);
            for (auto next_edge : extensions) {
                linear_path.push_back(next_edge);
                used_edges_.insert(next_edge);
                used_edges_.insert(gp_.g.conjugate(next_edge));
                extendForward(edge_set, linear_path, next_edge);
                extendBackward(edge_set, linear_path, next_edge);
            }
        }


        void extendBackward(const std::unordered_set<EdgeId> &edge_set,
                            std::deque<EdgeId> &linear_path, EdgeId e) {
            auto extensions = getReverseIntersect(e, edge_set);
            for (auto prev_edge : extensions) {
                linear_path.push_front(prev_edge);
                used_edges_.insert(prev_edge);
                used_edges_.insert(gp_.g.conjugate(prev_edge));
                extendForward(edge_set, linear_path, prev_edge);
                extendBackward(edge_set, linear_path, prev_edge);
            }
        }

        void extendForward(const GraphComponent<Graph> &comp,
                           std::deque<EdgeId> &linear_path, std::set<EdgeId> &used_edges) {
            EdgeId e = linear_path.back();
            if (comp.VertexInDegree(gp_.g.EdgeEnd(e)) != 1 ||
                    comp.VertexOutDegree(gp_.g.EdgeEnd(e)) != 1 ) {
                return;
            }
            auto extensions = getForwardIntersect(e, comp);
            VERIFY(extensions.size() == 1);
            EdgeId next_edge = extensions[0];
            if (used_edges.count(next_edge)) {
                return;
            }
            linear_path.push_back(next_edge);
            used_edges.insert(next_edge);
            used_edges.insert(gp_.g.conjugate(next_edge));
            extendForward(comp, linear_path, used_edges);
        }


        void extendBackward(const GraphComponent<Graph> &comp,
                            std::deque<EdgeId> &linear_path, std::set<EdgeId> &used_edges) {
            EdgeId e = linear_path.front();
            if (comp.VertexInDegree(gp_.g.EdgeStart(e)) != 1 ||
                    comp.VertexOutDegree(gp_.g.EdgeStart(e)) != 1 ) {
                return;
            }
            auto extensions = getReverseIntersect(e, comp);
            VERIFY(extensions.size() == 1);
            EdgeId prev_edge = extensions[0];
            if (used_edges.count(prev_edge)) {
                return;
            }
            linear_path.push_front(prev_edge);
            used_edges.insert(prev_edge);
            used_edges.insert(gp_.g.conjugate(prev_edge));
            extendBackward(comp, linear_path, used_edges);
        }

        //Here we suppose that this is single connected component
        bool IsSimplePath(const GraphComponent<Graph> &comp) const {
            size_t one_degree = 0;
            for (auto v : comp.vertices()) {
                auto edges = comp.IncidentEdges(v);
                VERIFY(edges.size() > 0);
                if (edges.size() >= 3) {
                    return false;
                }
                if (edges.size() == 1) {
                    one_degree++;
                }
            }
            return one_degree == 4;
        }

        //Here we suppose that this is single connected component
        bool IsSimpleCycle(const GraphComponent<Graph> &comp) const {
            return comp.vertices().size() == 2 && comp.edges().size() == 2;
        }

        size_t SplitComponentsOnSimplePaths(const GraphComponent<Graph> &comp) {
            size_t result = 0;
            //check self-conjugate
            std::set<EdgeId> used_edges;
            for (auto e : comp.edges()) {
                if (gp_.g.conjugate(e) == e) {
                    used_edges.insert(e);
                }
                if (gp_.g.EdgeEnd(e) == gp_.g.EdgeStart(e)) {
                    used_edges.insert(e);
                }
            }

            DEBUG(comp.edges());
            for (auto e : comp.edges()) {
                if (!used_edges.count(e)) {
                    DEBUG("Current edge " << e);
                    used_edges.insert(e);
                    used_edges.insert(gp_.g.conjugate(e));
                    std::deque<EdgeId> linear_path;
                    linear_path.push_back(e);
                    extendForward(comp, linear_path, used_edges);
                    extendBackward(comp, linear_path, used_edges);
                    if (gp_.g.EdgeStart(linear_path.front()) != gp_.g.EdgeEnd(linear_path.back())) {
                        std::vector<EdgeId> linear_path_vector;
                        for (auto e : linear_path) {
                            linear_path_vector.push_back(e);
                        }
                        temp_set_.AddPath(linear_path_vector, 1, true);
                        result++;
                    }
                }
            }
            return result;
        }


        bool SelfConjugate(const std::vector<EdgeId> &path1, const std::vector<EdgeId> &path2) {
            if (path1.size() != path2.size())
                return false;
            for (size_t i = 0; i < path1.size(); ++i) {
                if (path1[i] != gp_.g.conjugate(path2[path2.size() - 1 - i])) {
                    return false;
                }
            }
            return true;
        }

        void RemoveBulges(GraphComponent<Graph> &comp) {
            SplitComponentsOnSimplePaths(comp);
            std::vector<std::vector<EdgeId>> paths;
            for (auto p : temp_set_) {
                paths.push_back(p.path());
            }
            DEBUG(paths);
            std::vector<size_t> to_remove;
            for (size_t i = 0; i < paths.size(); ++i) {
                for (size_t j = i + 1; j < paths.size(); ++j) {
                    if (gp_.g.EdgeStart(paths[i].front()) == gp_.g.EdgeStart(paths[j].front()) &&
                            gp_.g.EdgeEnd(paths[i].back()) == gp_.g.EdgeEnd(paths[j].back())) {
                        if (SelfConjugate(paths[i], paths[j])) {
                            DEBUG("Self-conjugate bulges");
                            continue;
                        }

                        if (comp.Coverage(paths[i]) > comp.Coverage(paths[j])) {
                            to_remove.push_back(j);
                        } else {
                            to_remove.push_back(i);
                        }
                    }
                }
            }

            for (auto i : to_remove) {
                for (size_t j = 0; j < paths[i].size(); ++j) {
                    comp.RemoveEdge(paths[i][j]);
                }
            }
        }


    public:
        SimplePathExtractor(const debruijn_graph::conj_graph_pack &gp)
        : gp_(gp), temp_set_(gp.g) {
        }


        std::unordered_set<EdgeId> simplifyComponent(std::unordered_set<EdgeId> &edge_set, const std::string &barcode) {
            temp_set_.Clear();
            std::set<EdgeId> bad_edges;
            auto initial_component = GraphComponent<Graph>::FromEdges(gp_.g, edge_set, true);
            initial_component.ChangeCoverageProvider(gp_.barcode_coverage[0].GetBarcodeMap(barcode));
            DEBUG("Initial component");
            for (size_t i = 0; i < 2; ++i) {
                initial_component.RemoveIsolated();
                DEBUG("After remove isolated");
                initial_component.ClipTips();
                DEBUG("After clip tips");
                initial_component.FillGaps(30 * (i + 1));
                DEBUG("After fill gaps");
                RemoveBulges(initial_component);
                DEBUG("After remove bulges");
                initial_component.RemoveLowCoveredJunctions();
                temp_set_.Clear();
            }
            DEBUG("Fixed initial component");
            auto out_edge_set = initial_component.edges();
            return out_edge_set;
        }

        PathStorage<Graph>& getLongReads(std::unordered_set<EdgeId> &edge_set, const std::string &barcode) {
            temp_set_.Clear();
            std::set<EdgeId> bad_edges;
            auto initial_component = GraphComponent<Graph>::FromEdges(gp_.g, edge_set, true);
            initial_component.ChangeCoverageProvider(gp_.barcode_coverage[0].GetBarcodeMap(barcode));
            DEBUG("Initial component");
            for (size_t i = 0; i < 2; ++i) {
                initial_component.RemoveIsolated();
                DEBUG("After remove isolated");
                initial_component.ClipTips();
                DEBUG("After clip tips");
                initial_component.FillGaps(30 * (i + 1));
                DEBUG("After fill gaps");
                RemoveBulges(initial_component);
                DEBUG("After remove bulges");
                initial_component.RemoveLowCoveredJunctions();
                temp_set_.Clear();
            }
            DEBUG("Fixed initial component");
            edge_set = initial_component.edges();
            for (auto e : edge_set) {
                if (!used_edges_.count(e)) {
                    std::deque<EdgeId> linear_path;
                    linear_path.push_back(e);
                    used_edges_.insert(e);
                    used_edges_.insert(gp_.g.conjugate(e));
                    extendForward(edge_set, linear_path, e);
                    extendBackward(edge_set, linear_path, e);
                    auto component = GraphComponent<Graph>::FromEdges(gp_.g, linear_path, true);
                    component.ClipTips();
                    if (IsSimplePath(component)) {
                        DEBUG("Component is a simple path");
                        std::vector<EdgeId> path;
                        for (auto e : linear_path) {
                            path.push_back(e);
                        }
                        temp_set_.AddPath(path, 1, true);
                    } else if (IsSimpleCycle(component)) {
                        DEBUG("Component is a simple cycle");
                    } else {
                        DEBUG("Component is not a simple path");
                        size_t result = SplitComponentsOnSimplePaths(component);
                        DEBUG("Component is split on " << result << " paths");
                    }
                }
            }
            used_edges_.clear();
            return temp_set_;
        }
        protected:
            DECL_LOGGER("MoleculeExtraction");
    };

    class LongReadsCreator {
    private:
        const debruijn_graph::conj_graph_pack &gp_;
        SimplePathExtractor extractor_;
    public:
        LongReadsCreator(debruijn_graph::conj_graph_pack &gp)
        : gp_(gp), extractor_(gp) {
        }

        void extractLongReadsPE(debruijn_graph::conj_graph_pack &graph_pack, path_extend::PathContainer &path_set, const std::vector<EdgeId> &out_edges, const std::string &barcode, const lib_t& lib_10x) {


            path_extend::PathExtendParamsContainer params(cfg::get().ds,
                                              cfg::get().pe_params,
                                              cfg::get().ss,
                                              cfg::get().output_dir,
                                              cfg::get().mode,
                                              cfg::get().uneven_depth,
                                              cfg::get().avoid_rc_connections,
                                              cfg::get().use_scaffolder);

            DEBUG("graph_pack.barcode_indices[0].size() - " << graph_pack.barcode_indices[0].size());
            estimate_distance_molecule_extraction(graph_pack, lib_10x, graph_pack.barcode_indices[0], graph_pack.barcode_clustered_indices[0], out_edges);
            DEBUG("graph_pack.barcode_clustered_indices[0].size() - " << graph_pack.barcode_clustered_indices[0].size());

            if (graph_pack.barcode_clustered_indices[0].size() == 0)
                return;
            auto paired_lib = path_extend::MakeNewLib(graph_pack.g, lib_10x, graph_pack.barcode_clustered_indices[0]);
            std::shared_ptr<path_extend::CoverageAwareIdealInfoProvider> iip = std::make_shared<path_extend::CoverageAwareIdealInfoProvider>(graph_pack.g, paired_lib, lib_10x.data().unmerged_read_length);
            auto wc = std::make_shared<path_extend::PathCoverWeightCounter>(graph_pack.g, paired_lib, params.pset.normalize_weight,
                                                          params.pset.extension_options.single_threshold,
                                                          iip);

            auto opts = params.pset.extension_options;
            auto extension_chooser_1 = std::make_shared<path_extend::TrivialExtensionChooser>(graph_pack.g);

            auto extension_chooser_2 = std::make_shared<path_extend::SimpleExtensionChooser>(graph_pack.g, wc,
                                                                         opts.weight_threshold,
                                                                         opts.priority_coeff);

            path_extend::GraphCoverageMap cover_map(graph_pack.g, out_edges.size());

            path_extend::UniqueData unique_data;
            path_extend::UsedUniqueStorage used_unique_storage(unique_data.main_unique_storage_, graph_pack.g);

            std::set<EdgeId> edge_set(out_edges.begin(), out_edges.end());
            auto extender_1 = std::make_shared<path_extend::GoodEdgeExtender>(graph_pack, cover_map,
                                                                              used_unique_storage,
                                                                              extension_chooser_1,
                                                                              paired_lib->GetISMax(),
                                                                              edge_set,
                                                                              true, /*investigate loops*/
                                                                              true /*use short loop coverage resolver*/,
                                                                              opts.weight_threshold);
            auto extender_2 = std::make_shared<path_extend::GoodEdgeExtender>(graph_pack, cover_map,
                                                                            used_unique_storage,
                                                                            extension_chooser_2,
                                                                            paired_lib->GetISMax(),
                                                                            edge_set,
                                                                            true, /*investigate loops*/
                                                                            true /*use short loop coverage resolver*/,
                                                                            opts.weight_threshold);
            path_extend::PathExtendResolver resolver(graph_pack.g);
            auto seeds = resolver.MakeSeedsFromEdgeSet(out_edges);
            seeds.SortByLength();

            typedef std::vector<std::shared_ptr<path_extend::PathExtender>> Extenders;
            Extenders extenders {extender_1, extender_2};
            path_extend::CompositeExtender composite_extender(graph_pack.g, cover_map,
                                                              used_unique_storage,
                                                              extenders);
            path_set = resolver.ExtendSeeds(seeds, composite_extender);
            DEBUG("path_set.size() - " << path_set.size());
        }

        void extractEdges(std::vector<MappingPath<EdgeId>> &paths, std::vector<EdgeId> &out_edges, const std::string &barcode) {
            if (paths.size() < cfg::get().pe_params.param_set.rna_10x.min_cloud_size)
                return;
            std::unordered_set<EdgeId> edge_set;
            for (auto const& path : paths) {
                std::vector<EdgeId> edges = path.simple_path();
                for (auto e : edges) {
                    edge_set.insert(e);
                    edge_set.insert(gp_.g.conjugate(e));
                }
            }
            auto out_edges_set = extractor_.simplifyComponent(edge_set, barcode);
            for (auto e : out_edges_set) {
                out_edges.push_back(e);
            }
        }

        bool extractLongReads(std::vector<MappingPath<EdgeId>> &paths, PathStorage<Graph> &path_set, const std::string &barcode) {
            if (paths.size() < cfg::get().pe_params.param_set.rna_10x.min_cloud_size)
                return false;
            std::map<EdgeId, int> edge_map;
            std::unordered_set<EdgeId> edge_set;
            for (auto const& path : paths) {
                std::vector<EdgeId> edges = path.simple_path();
                for (auto e : edges) {
                    edge_set.insert(e);
                    edge_set.insert(gp_.g.conjugate(e));
                    if (edge_map.find(e) != edge_map.end()) {
                        edge_map.emplace(e, 0);
                        edge_map.emplace(gp_.g.conjugate(e), 0);
                    }
                    edge_map[e] += 1;
                    edge_map[gp_.g.conjugate(e)] += 1;
                }
            }

            std::string edges = "";
            auto &temp_set = extractor_.getLongReads(edge_set, barcode);
            DEBUG("BARCODE: " << barcode << ", number of reads - " << paths.size() << ", number of edges - " << edge_map.size() / 2 << " number of paths " << temp_set.size());
            for (auto em : edge_map) {
                auto e = em.first;
                if (e <= gp_.g.conjugate(e)) {
                    edges += std::to_string(gp_.g.int_id(e)) + ", ";
                }
            }

            DEBUG(temp_set.size() << " long reads extracted");
            for (const auto& path : temp_set) {
                path.add_barcodes(barcode);
                DEBUG("Cut from beginning " << gp_.barcode_coverage[0].GetLeftMostPosition(path.path()[0], barcode));
                DEBUG("Cut from end " << gp_.barcode_coverage[0].GetRightMostPosition(path.path()[path.path().size() - 1], barcode));
                DEBUG("Add path with barcodes " << barcode);
                path_set.AddPath(path.path(), path.weight(), false, path.get_barcodes(), gp_.barcode_coverage[0].GetLeftMostPosition(path.path()[0], barcode), gp_.barcode_coverage[0].GetRightMostPosition(path.path()[path.path().size() - 1], barcode));
            }
            DEBUG("BARCODE END");
            temp_set.Clear();
            return true;
        }
    protected:
        DECL_LOGGER("LongReadCreator");
    };



    std::string GetTenXBarcodeFromRead(const io::PairedRead &read) {
        std::string delimeter = "BX:Z:";
        size_t start_pos = read.first().name().find(delimeter);
        std::string barcode = "";
        if (start_pos == std::string::npos) {
            return barcode;
        }
        for (int i = start_pos; i < read.first().name().length(); ++i) {
            if (read.first().name()[i] != '1') {
                barcode.push_back(read.first().name()[i]);
            } else {
                barcode.push_back(read.first().name()[i]);
                break;
            }
        }
        TRACE(barcode);
        return barcode;
    }

    void Process10XLibrary(debruijn_graph::conj_graph_pack &graph_pack, const lib_t& lib_10x) {
        auto mapper = MapperInstance(graph_pack);
        auto stream = io::paired_easy_reader(lib_10x, false, 0);
        std::set<std::string> used_barcodes;
        std::string current_barcode = "";
        io::PairedRead read;
        size_t counter = 0;
        size_t passed_counter = 0;
        size_t failed_counter = 0;
        std::vector<MappingPath<EdgeId>> paths;
        path_extend::PathContainer& long_reads = graph_pack.contig_paths;

        LatePairedIndexFiller::WeightF weight = [&](const std::pair<EdgeId, EdgeId> &,
                     const MappingRange&, const MappingRange&) {
            return 1.;
        };
        unsigned round_thr = 0;
        LatePairedIndexFiller pif(graph_pack.g,
                                  weight, round_thr,
                                  graph_pack.barcode_indices[0]);

        pif.StartProcessLibrary(0);
        PathStorage<Graph> long_reads_temp_storage(graph_pack.g);
        path_extend::PathContainer long_reads_temp_container;
        LongReadsCreator extractor(graph_pack);
        while (!stream.eof()) {
            stream >> read;
            std::string barcode_string = GetTenXBarcodeFromRead(read);
            if (barcode_string != "") {
                if (barcode_string != current_barcode && !paths.empty()) {
                    DEBUG("graph_pack.barcode_indices[0].size() - " << graph_pack.barcode_indices[0].size());
                    INFO("Processing barcode " << current_barcode);
                    pif.StopProcessLibrary();
                    std::vector<EdgeId> good_edges;

                    extractor.extractEdges(paths, good_edges, current_barcode);
                    if (paths.size() > cfg::get().pe_params.param_set.rna_10x.min_cloud_size)
                        extractor.extractLongReadsPE(graph_pack, long_reads_temp_container, good_edges, current_barcode, lib_10x);

                    for (auto path : long_reads_temp_container) {
                        DEBUG(path.first->ToVector()[0]);
                        DEBUG(path.first->ToVector()[path.first->ToVector().size() - 1]);
                        DEBUG(graph_pack.barcode_coverage[0].GetLeftMostPosition(path.first->ToVector()[0], current_barcode));
                        DEBUG(graph_pack.barcode_coverage[0].GetRightMostPosition(path.first->ToVector()[path.first->ToVector().size() - 1], current_barcode));
                        long_reads_temp_storage.AddPath(path.first->ToVector(), 1, true, current_barcode, graph_pack.barcode_coverage[0].GetLeftMostPosition(path.first->ToVector()[0], current_barcode),
                                                        graph_pack.barcode_coverage[0].GetRightMostPosition(path.first->ToVector()[path.first->ToVector().size() - 1], current_barcode));
                    }

                    if (long_reads_temp_container.size())
                        ++passed_counter;
                    else
                        ++failed_counter;
                    long_reads_temp_container.clear();
                    if (used_barcodes.count(current_barcode) && current_barcode != "") {
                        WARN("Path with " << current_barcode << " barcode was previously extracted");
                    }
                    else {
                        used_barcodes.insert(current_barcode);
                        if (used_barcodes.size() % 10000 == 0) {
                            INFO(used_barcodes.size() << " barcodes processed");
                        }
                        graph_pack.barcode_coverage[0].Clear(current_barcode);
                    }
                    graph_pack.barcode_indices[0].Clear();
                    graph_pack.barcode_clustered_indices[0].Clear();
                    paths.clear();
                    pif.StartProcessLibrary(0);
                }
                current_barcode = barcode_string;
                const auto &path1 = mapper->MapRead(read.first());
                const auto &path2 = mapper->MapRead(read.second());
                paths.push_back(path1);
                paths.push_back(path2);
                pif.ProcessPairedRead(0, read, path1, path2);
            }
            counter++;
            VERBOSE_POWER_T2(counter, 100, "Processed " << counter << " reads.");
        }
        DEBUG("Processing barcode " << current_barcode);
        std::vector<EdgeId> good_edges;
        extractor.extractEdges(paths, good_edges, current_barcode);
        extractor.extractLongReadsPE(graph_pack, long_reads_temp_container, good_edges, current_barcode, lib_10x);
        if (long_reads_temp_container.size())
            ++passed_counter;
        else
            ++failed_counter;

        for (auto path : long_reads_temp_container) {
            long_reads_temp_storage.AddPath(path.first->ToVector(), 1, true, current_barcode, graph_pack.barcode_coverage[0].GetLeftMostPosition(path.first->ToVector()[0], current_barcode),
                                            graph_pack.barcode_coverage[0].GetRightMostPosition(path.first->ToVector()[path.first->ToVector().size() - 1], current_barcode));
        }
        long_reads_temp_container.clear();

        INFO("Starting to save all paths found");
        std::vector<PathInfo<Graph>> debug_path;
        long_reads_temp_storage.SaveAllPaths(debug_path);

        for (auto p : debug_path) {
            path_extend::BidirectionalPath *path = new path_extend::BidirectionalPath(graph_pack.g);
            auto conj = new path_extend::BidirectionalPath(graph_pack.g);
            long_reads.AddPair(path, conj);
            for (auto e : p.path()) {
                path->PushBack(e);
            }
            path->SetBarcode(p.get_barcodes());
            conj->SetBarcode(p.get_barcodes());
            path->SetWeight(p.weight());
            conj->SetWeight(p.weight());
            path->SetCutFromBeginning(p.get_cut_from_begin());
            path->SetCutFromEnd(p.get_cut_from_end());
            conj->SetCutFromBeginning(p.get_cut_from_end());
            conj->SetCutFromEnd(p.get_cut_from_begin());
            DEBUG(p.get_barcodes());
            DEBUG(p.get_cut_from_begin());
            DEBUG(p.get_cut_from_end());
            //path->PrintINFO();
        }

        INFO("Total reads processed: " << counter);
        INFO("Total barcodes passed: " << passed_counter << ", failed " << failed_counter);
        INFO(long_reads.size() << " paths totally extracted");
        long_reads.FilterPaths(path_extend::LengthPathCondition(300));
        INFO(long_reads.size() << " paths after filtering");
        path_extend::ContigWriter writer(graph_pack.g, std::make_shared<path_extend::BarcodeContigNameGenerator>());
        writer.OutputPathsRNA(long_reads, cfg::get().output_dir + "/extracted.fasta");
    }

    void MoleculeExtractionStage::run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
        INFO("Molecule extraction for RNA stage started");
        config::dataset& dataset_info = cfg::get_writable().ds;
        if (!Has10XRNAReads(dataset_info)) {
            INFO("Read cloud libraries have not been found. Skipping molecule extraction stage.")
            return;
        }

        graph_pack.InitRRIndices();
        graph_pack.EnsureBasicMapping();

        for (auto& lib : dataset_info.reads) {
            if (lib.type() == io::LibraryType::RNA10x) {
                Process10XLibrary(graph_pack, lib);
                lib.data().single_reads_mapped = true;
            }
        }
        cfg::get_writable().use_single_reads = false;
        INFO("Molecule extraction for RNA stage ended");
    }

} //debruijn_graph

