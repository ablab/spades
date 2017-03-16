////***************************************************************************
////* Copyright (c) 2011-2014 Saint-Petersburg Academic University
////* All Rights Reserved
////* See file LICENSE for details.
////****************************************************************************
//
//#pragma once
//
//namespace cap {
//
//template<class Graph>
//class LabelFilter : public GraphComponentFilter<Graph> {
//    typedef typename Graph::VertexId VertexId;
//    const EdgesPositionHandler<Graph> &edge_pos_;
//public:
//    LabelFilter(const Graph &graph, const EdgesPositionHandler<Graph> &edge_pos) :
//            GraphComponentFilter<Graph>(graph), edge_pos_(edge_pos) {
//    }
//
//    virtual bool Check(const vector<typename Graph::VertexId> &component) const {
//        set<VertexId> cset(component.begin(), component.end());
//        for(auto vit = component.begin(); vit != component.end(); ++vit) {
//            auto out = this->graph().OutgoingEdges(*vit);
//            for(auto eit = out.begin(); eit != out.end(); ++eit) {
//                if(cset.count(this->graph().EdgeEnd(*eit)) > 0) {
//                    auto labels = edge_pos_.GetEdgePositions(*eit);
//                    for(auto it = labels.begin(); it != labels.end(); ++it) {
//                        if(it->first == "ref_0" || it->first == "ref_1")
//                            return true;
//                    }
//                }
//            }
//        }
//        return false;
//    }
//};
//
////todo как обобщить на сравнение с геномом???
//template<class gp_t>
//class IDBADiffAnalyzer {
//private:
//    typedef typename gp_t::graph_t Graph;
//    typedef typename gp_t::index_t Index;
//    typedef typename gp_t::seq_t Kmer;
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//    typedef io::SingleRead Contig;
//    typedef io::IReader<io::SingleRead> ContigStream;
//    typedef    io::MultifileReader<io::SingleRead> CompositeStream;
//    typedef debruijn_graph::BasicSequenceMapper<Graph, Index> Mapper;
//
//    const gp_t& gp_;
//    const ColorHandler<Graph>& coloring_;
//    Mapper mapper_;
//
//    string good_assembly_prefix_;
//    string bad_assembly_prefix_;
//    string dir_;
//    map<string, Sequence> contigs_map_;
//
//    bool StartsWith(const string& s, const string& prefix) {
//        return boost::starts_with(s, prefix);
//    }
//
//    void CollectContigs(ContigStream& good_assembly, ContigStream& bad_assembly) {
//        CompositeStream composite(good_assembly, bad_assembly);
//        composite.reset();
//        while(!composite.eof()) {
//            Contig c;
//            composite >> c;
//            contigs_map_[c.name()] = c.sequence();
//        }
//        composite.reset();
//    }
//
//    set<string> CollectBadContigIdsAlongPath(const vector<EdgeId>& path) {
//        DEBUG("Collecting intersecting contigs from position marks");
//        set<string> answer;
//        for (auto it = path.begin(); it != path.end(); ++it) {
//            vector<EdgePosition> positions = gp_.edge_pos.GetEdgePositions(*it);
//            for (auto pos_it = positions.begin(); pos_it != positions.end(); ++pos_it) {
//                string id = pos_it->contigId;
//                if (StartsWith(id, bad_assembly_prefix_)) {
//                    answer.insert(id);
//                }
//            }
//        }
//        DEBUG("Collected " << answer.size() << " contigs");
//        return answer;
//    }
//
//    size_t Intersection(const set<EdgeId>& s1, const set<EdgeId>& s2) {
//        size_t ans = 0;
//        for (auto it = s1.begin(); it != s1.end(); ++it) {
//            if (s2.count(*it) > 0) {
//                ans++;
//            }
//        }
//        return ans;
//    }
//
//    set<EdgeId> AsSet(vector<EdgeId> v) {
//        return set<EdgeId>(v.begin(), v.end());
//    }
//
//    vector<EdgeId> MappingEdgeVector(const string& contig_id) {
//        VERIFY(contigs_map_.find(contig_id) != contigs_map_.end());
//        return mapper_.MapSequence(contigs_map_[contig_id]).simple_path();
//    }
//
//    set<EdgeId> MappingEdgeSet(const string& contig_id) {
//        return AsSet(MappingEdgeVector(contig_id));
//    }
//
//    string FindBestBadContigWRTPath(const set<string>& contigs, const vector<EdgeId>& path) {
//        DEBUG("Looking for best contig")
//        set<EdgeId> path_edges(path.begin(), path.end());
//        size_t best_intersection = 0;
//        string best_contig = "";
//        for (auto it = contigs.begin(); it != contigs.end(); ++it) {
//            size_t intersect = Intersection(MappingEdgeSet(*it), path_edges);
//            if (intersect > best_intersection) {
//                best_intersection = intersect;
//                best_contig = *it;
//            }
//        }
//        DEBUG("Best contig is " << best_contig);
//        return best_contig;
//    }
//
//    bool InnerVertex(VertexId v, const vector<EdgeId>& path) {
//        if (path.empty())
//            return false;
//        for (size_t i = 0; i < path.size() - 1; ++i) {
//            if (gp_.g.EdgeEnd(path[i]) == v)
//                return true;
//        }
//        return false;
//    }
//
//    VertexId FirstBranchingVertex(const vector<EdgeId>& bad_path, const vector<EdgeId>& good_path) {
//        for (auto it = bad_path.begin(); it != bad_path.end(); ++it) {
//            if (InnerVertex(gp_.g.EdgeStart(*it), good_path)) {
//                return gp_.g.EdgeStart(*it);
//            }
//            if (InnerVertex(gp_.g.EdgeEnd(*it), good_path)) {
//                return gp_.g.EdgeEnd(*it);
//            }
//        }
//        return (VertexId) NULL;
//    }
//
//    VertexId LastBranchingVertex(const vector<EdgeId>& bad_path, const vector<EdgeId>& good_path) {
//        for (auto it = bad_path.rbegin(); it != bad_path.rend(); ++it) {
//            if (InnerVertex(gp_.g.EdgeEnd(*it), good_path)) {
//                return gp_.g.EdgeEnd(*it);
//            }
//            if (InnerVertex(gp_.g.EdgeStart(*it), good_path)) {
//                return gp_.g.EdgeStart(*it);
//            }
//        }
//        return (VertexId) NULL;
//    }
//
//    bool SingleAssemblyEdge(EdgeId e, const string& prefix) {
//        vector<EdgePosition> positions = gp_.edge_pos.GetEdgePositions(e);
//        for (auto it = positions.begin(); it != positions.end(); ++it) {
//            if (!StartsWith(it->contigId, prefix)) {
//                return false;
//            }
//        }
//        return true;
//    }
//
//    bool ContainsSingleAssemblyEdge(const vector<EdgeId> edges, const string& prefix) {
//        for (auto it = edges.begin(); it != edges.end(); ++it)
//            if (SingleAssemblyEdge(*it, prefix))
//                return true;
//        return false;
//    }
//
//    vector<EdgeId> SingleAssemblyEdges(const vector<EdgeId> edges, const string& prefix) {
//        vector<EdgeId> answer;
//        for (auto it = edges.begin(); it != edges.end(); ++it)
//            if (SingleAssemblyEdge(*it, prefix))
//                answer.push_back(*it);
//        return answer;
//    }
//
//    vector<EdgeId> IncidentEdges(VertexId v) {
//        vector<EdgeId> ans;
//        utils::push_back_all(ans, gp_.g.IncomingEdges(v));
//        utils::push_back_all(ans, gp_.g.OutgoingEdges(v));
//        return ans;
//    }
//
//    vector<EdgeId> IncidentEdgesInPath(VertexId v, const vector<EdgeId>& good_contig_path) {
//        vector<EdgeId> ans;
//        vector<EdgeId> adj = IncidentEdges(v);
//        for (size_t i = 0; i < adj.size(); ++i)
//            if (find(good_contig_path.begin(), good_contig_path.end(), adj[i]) != good_contig_path.end())
//                ans.push_back(adj[i]);
//        return ans;
//    }
//
//    void ReportLocality(VertexId v, const vector<EdgeId>& good_contig_path, const string& best_contig, const Contig& c, const string& folder) {
//        using namespace visualization;
//        make_dir(folder);
//        LengthIdGraphLabeler<Graph> basic_labeler(gp_.g);
//        EdgePosGraphLabeler<Graph> pos_labeler(gp_.g, gp_.edge_pos);
//
//        CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
//
//        LabelFilter<typename gp_t::graph_t> lf(gp_.g, gp_.edge_pos);
//        string file_name = folder + c.name() + "_|_" + best_contig + ".dot";
//        EdgeId edge = IncidentEdgesInPath(v, good_contig_path).front();
//        GraphComponent<Graph> component = omnigraph::EdgeNeighborhood(gp_.g, edge);
//        if(lf.Check(vector<VertexId>(component.v_begin(), component.v_end()))) {
//            WriteComponent(component, file_name, BorderDecorator<Graph>::GetInstance(component, coloring_.GetInstance()), labeler);
//        }
//    }
//
//    size_t LengthForward(EdgeId e, const vector<EdgeId>& good_contig_path) {
//        size_t ans = 0;
//        bool passed = false;
//        for (auto it = good_contig_path.begin(); it != good_contig_path.end(); ++it) {
//            if (*it == e)
//                passed = true;
//            if (passed)
//                ans += gp_.g.length(*it);
//        }
//        return ans;
//    }
//
//    size_t LengthBackward(EdgeId e, const vector<EdgeId>& good_contig_path) {
//        size_t ans = 0;
//        bool passed = false;
//        for (auto it = good_contig_path.rbegin(); it != good_contig_path.rend(); ++it) {
//            if (*it == e)
//                passed = true;
//            if (passed)
//                ans += gp_.g.length(*it);
//        }
//        return ans;
//    }
//
//    size_t LengthToEndOnGoodContig(VertexId v, const vector<EdgeId>& good_contig_path, const set<EdgeId>& best_alt_contig_path) {
//        vector<EdgeId> adj_in_path = IncidentEdgesInPath(v, good_contig_path);
//        for (auto it = adj_in_path.begin(); it != adj_in_path.end(); ++it) {
//            if (best_alt_contig_path.count(*it) == 0) {
//                if (v == gp_.g.EdgeStart(*it)) {
//                    return LengthForward(*it, good_contig_path);
//                } else {
//                    VERIFY(v == gp_.g.EdgeEnd(*it));
//                    return LengthBackward(*it, good_contig_path);
//                }
//            }
//        }
////        VERIFY(false);
//        WARN("Something strange for vertex " << v);
//        return 0;
//    }
//
//    void ClassifyAndReportBreak(VertexId v, const vector<EdgeId>& good_contig_path, const string& best_contig, const Contig& c) {
//        DEBUG("Trying to classify break");
//        if (gp_.g.EdgeStart(good_contig_path.front()) == v || gp_.g.EdgeEnd(good_contig_path.back()) == v) {
//            DEBUG("Vertex was an end of initial contig");
//            return;
//        }
//        if (ContainsSingleAssemblyEdge(IncidentEdgesInPath(v, good_contig_path), good_assembly_prefix_)) {
//            DEBUG("Vertex has adjacent \"good\" assembly edge");
//            if (SingleAssemblyEdges(IncidentEdgesInPath(v, good_contig_path), good_assembly_prefix_).size() == 1) {
//                EdgeId e = SingleAssemblyEdges(IncidentEdgesInPath(v, good_contig_path), good_assembly_prefix_).front();
//                if (e == good_contig_path.front() || e == good_contig_path.back()) {
//                    DEBUG("This edge is at the end of initial contig");
//                    DEBUG("Skipping");
//                    return;
//                }
//            }
//            DEBUG("Reporting locality of vertex " << gp_.g.str(v) << " as possible coverage gap");
//            ReportLocality(v, good_contig_path, best_contig, c, dir_ + "/coverage_gaps/");
//            return;
//        }
//        if (ContainsSingleAssemblyEdge(IncidentEdges(v), bad_assembly_prefix_)
//                && !ContainsSingleAssemblyEdge(IncidentEdgesInPath(v, good_contig_path), good_assembly_prefix_)) {
//            DEBUG("Reporting locality of vertex " << gp_.g.str(v) << " as possible EC problem");
//            ReportLocality(v, good_contig_path, best_contig, c, dir_ + "/ec_problem/");
//            return;
//        }
//        if (!ContainsSingleAssemblyEdge(IncidentEdges(v), bad_assembly_prefix_)
//                && !ContainsSingleAssemblyEdge(IncidentEdgesInPath(v, good_contig_path), good_assembly_prefix_)) {
//            DEBUG("Possible RR problem. Checking remaining length");
//            if (LengthToEndOnGoodContig(v, good_contig_path, MappingEdgeSet(best_contig)) > 10000) {
//                DEBUG("Check ok");
//                DEBUG("Reporting locality of vertex " << gp_.g.str(v) << " as possible RR problem");
//                ReportLocality(v, good_contig_path, best_contig, c, dir_ + "/rr_problem/");
//                return;
//            } else {
//                DEBUG("Check fail, won't report");
//            }
//        }
//        DEBUG("Unclassified problem type");
//    }
//
//    void AnalyzeBadContigsWRTPath(const set<string>& contigs, const vector<EdgeId>& path, const Contig& c) {
//        string best_contig = FindBestBadContigWRTPath(contigs, path);
//        if (best_contig == "")
//            return;
//        vector<EdgeId> path_edges = MappingEdgeVector(best_contig);
//
//        DEBUG("Best contig mapped to path: " << gp_.g.str(path_edges));
//
//        if (path_edges.empty() || !CheckContiguous(gp_.g, path_edges)) {
//            WARN("Path for best contig " << best_contig << " wasn't continuous");
//            return;
//        }
//        DEBUG("Looking for first branching vertex");
//        VertexId first = FirstBranchingVertex(path_edges, path);
//        if(first != VertexId(NULL)) {
//            DEBUG("First branching vertex is " << gp_.g.str(first));
//            ClassifyAndReportBreak(first, path, best_contig, c);
//        } else {
//            DEBUG("Failed to find first branching vertex");
//        }
//        DEBUG("Looking for last branching vertex");
//        VertexId last = LastBranchingVertex(path_edges, path);
//        if(last != VertexId(NULL)) {
//            DEBUG("Last branching vertex is " << gp_.g.str(last));
//            ClassifyAndReportBreak(last, path, best_contig, c);
//        } else {
//            DEBUG("Failed to find last branching vertex");
//        }
//    }
//
//    void AnalyzeGoodContig(const Contig& c) {
//        DEBUG("Analyzing contig " << c.name());
//
//        vector<EdgeId> path_edges = MappingEdgeVector(c.name());
//        DEBUG("Contig mapped to path: " << gp_.g.str(path_edges));
//
//        if (path_edges.empty() || !CheckContiguous(gp_.g, path_edges)) {
//            WARN("Path for good contig " << c.name() << " wasn't continuous");
//            return;
//        }
//
//        set<string> bad_contig_ids = CollectBadContigIdsAlongPath(path_edges);
//        AnalyzeBadContigsWRTPath(bad_contig_ids, path_edges, c);
//    }
//
//public:
//    IDBADiffAnalyzer(const gp_t& gp,
////            const EdgesPositionHandler<Graph>& pos,
//                   const ColorHandler<Graph>& coloring,
//                   const string& good_assembly_prefix,
//                   const string& bad_assembly_prefix,
//                   const string& dir)
//    : gp_(gp)/*, pos_(pos)*/,
//      coloring_(coloring),
//      mapper_(gp.g, gp.index, gp.kmer_mapper, gp.k_value + 1),
//      good_assembly_prefix_(good_assembly_prefix),
//      bad_assembly_prefix_(bad_assembly_prefix),
//      dir_(dir) {
//        DEBUG("\"Good\" assembly prefix " << good_assembly_prefix);
//        DEBUG("\"Bad\" assembly prefix " << bad_assembly_prefix);
//    }
//
//    void Analyze(ContigStream& good_assembly, ContigStream& bad_assembly) {
//        CollectContigs(good_assembly, bad_assembly);
//        while (!good_assembly.eof()) {
//            Contig c;
//            good_assembly >> c;
//            AnalyzeGoodContig(c);
//        }
//    }
//
//    DECL_LOGGER("IDBADiffAnalyzer");
//};
//
////investigates if red edges can close gaps in blue assembly
//template<class Graph>
//class GapComparativeAnalyzer {
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//
//    const Graph& g_;
//    const ColorHandler<Graph>& coloring_;
//    const EdgesPositionHandler<Graph>& pos_;
//
//    bool PurpleOrRed(EdgeId e) {
//        return coloring_.Color(e) == kRedColorSet
//                || coloring_.Color(e) == kVioletColorSet;
//    }
//
//    bool CheckVertexCondition(VertexId v) {
//        return g_.CheckUniqueOutgoingEdge(v) && g_.CheckUniqueIncomingEdge(v)
//                && PurpleOrRed(g_.GetUniqueOutgoingEdge(v))
//                && PurpleOrRed(g_.GetUniqueIncomingEdge(v));
//    }
//
//    void ReportEdge(EdgeId e, const string& folder) {
//        using namespace visualization;
//        INFO(
//                "Can close gap between edges " << g_.str(g_.GetUniqueIncomingEdge(g_.EdgeStart(e))) << " and " << g_.str(g_.GetUniqueOutgoingEdge(g_.EdgeEnd(e))) << " with edge " << g_.str(e));
//        LengthIdGraphLabeler<Graph> basic_labeler(g_);
//        EdgePosGraphLabeler<Graph> pos_labeler(g_, pos_);
//
//        CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
//        GraphComponent<Graph> component = omnigraph::EdgeNeighborhood(g_, e);
//        auto colorer = coloring_.ConstructColorer(component);
//        visualization::visualization_utils::WriteComponent(component, folder + std::to_string(g_.int_id(e)) + "_loc.dot", colorer, labeler);
//    }
//
////    bool CheckEdges(const vector<EdgeId>& edges) {
////        set<TColorSet> colors;
////        for (auto it = edges.begin(); it != edges.end(); ++it) {
////            colors.insert(coloring_.Color(*it));
////        }
////        return edges.size() == 1
////                || (edges.size() == 2 && colors.count(kBlueColor) == 1);
////    }
//
//public:
//    GapComparativeAnalyzer(const Graph& g, const ColorHandler<Graph>& coloring,
//            const EdgesPositionHandler<Graph>& pos) :
//            g_(g), coloring_(coloring), pos_(pos) {
//    }
//
//    void ReportPotentialGapsCloses(const string& folder) {
//        make_dir(folder);
//        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
//            if (coloring_.Color(*it) == kRedColorSet
//                    && CheckVertexCondition(g_.EdgeStart(*it))
//                    && CheckVertexCondition(g_.EdgeEnd(*it))) {
//                ReportEdge(*it, folder);
//            }
//        }
//    }
//
////    void ReportDeepGapsCloses(const string& folder) {
////        make_dir(folder);
////        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
////            if (coloring_.Color(*it) == kRedColorSet
////                    && CheckEdges(g_.OutgoingEdges(g_.EdgeStart(*it)))
////                    && CheckEdges(g_.IncomingEdges(g_.EdgeEnd(*it)))
////                    && ContainsTip()) {
////
////            }
////        }
////    }
//};
//
//
//}
