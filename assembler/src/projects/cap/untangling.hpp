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
//struct bp_graph_pack {
//    typedef Graph graph_t;
//    typedef string contig_id_t;
//    typedef typename Graph::EdgeId EdgeId;
//    Graph g;
//    ColorHandler<Graph> coloring;
//    map<contig_id_t, vector<EdgeId>> red_paths;
//    map<contig_id_t, vector<EdgeId>> blue_paths;
//    EdgesPositionHandler<Graph> edge_pos;
//
//    bp_graph_pack(size_t k) :
//            g(k), coloring(g), edge_pos(g) {
//
//    }
//};
//
//template<class gp_t>
//class UntangledGraphContigMapper {
//    typedef typename gp_t::graph_t Graph;
//    typedef typename Graph::EdgeId EdgeId;
//    const bp_graph_pack<Graph>& bp_gp_;
//
//public:
//    UntangledGraphContigMapper(const bp_graph_pack<Graph>& bp_gp) :
//            bp_gp_(bp_gp) {
//
//    }
//
//    MappingPath<EdgeId> MapRead(const io::SingleRead &read) const {
//        auto it = bp_gp_.red_paths.find(read.name());
//        if (it != bp_gp_.red_paths.end()) {
//            return TrivialMappingPath(it->second);
//        }
//        it = bp_gp_.blue_paths.find(read.name());
//        if (it != bp_gp_.blue_paths.end()) {
//            return TrivialMappingPath(it->second);
//        }VERIFY(false);
//        return MappingPath<EdgeId>();
//    }
//
//};
//
//template<class gp_t>
//class UntangledGraphConstructor {
//private:
//    typedef typename gp_t::graph_t Graph;
//    typedef typename Graph::VertexId VertexId;
//    typedef typename Graph::EdgeId EdgeId;
//
//    static const size_t k = gp_t::k_value;
//
//    const gp_t& old_gp_;
//    const ColorHandler<Graph>& old_coloring_;
//    bp_graph_pack<Graph>& new_gp_;
//    restricted::map<EdgeId, EdgeId> purple_edge_mapping_;
//    restricted::map<VertexId, VertexId> vertex_mapping_;
//    //todo draw in different color!
//    restricted::set<VertexId> artificial_vertices_;
//    set<string> processed_contigs_;
//
//    //todo test that!!!
//    string ConjugateContigId(const string& contig_id) {
//        string answer;
//        if (contig_id.substr(contig_id.size() - 3, 3) == "_RC")
//            answer = contig_id.substr(0, contig_id.size() - 3);
//        else
//            answer = contig_id + "_RC";
//        DEBUG("Conjugate to " << contig_id << " is " << answer);
//        return answer;
//    }
//
//    void AddToProcessed(const string& contig_id) {
//        processed_contigs_.insert(contig_id);
//        processed_contigs_.insert(ConjugateContigId(contig_id));
//    }
//
//    VertexId GetStartVertex(const Path<EdgeId> &path, size_t i) {
//        if (i != 0 || path.start_pos() == 0)
//            return vertex_mapping_[old_gp_.g.EdgeStart(path[i])];
//        else {
//            //todo discuss with Anton!!!
//            VertexId art_v = new_gp_.g.AddVertex();
//            WARN("Art vertex added")
////            VERIFY(false);
//            artificial_vertices_.insert(art_v);
//            return art_v;
//        }
//    }
//
//    VertexId GetEndVertex(const Path<EdgeId> &path, size_t i) {
//        if (i != path.size() - 1 || path.end_pos() == old_gp_.g.length(path[i]))
//            return vertex_mapping_[old_gp_.g.EdgeEnd(path[i])];
//        else {
//            //todo discuss with Anton!!!
//            VertexId art_v = new_gp_.g.AddVertex();
//            WARN("Art vertex added")
////            VERIFY(false);
//            artificial_vertices_.insert(art_v);
//            return art_v;
//        }
//    }
//
//    void Untangle(ContigStream& stream, TColorSet color) {
//        io::SingleRead read;
//        stream.reset();
//        set<string> processed;
//        while (!stream.eof()) {
//            stream >> read;
//            //todo can look at new_gp_.*_paths keys
//            if (processed.count(read.name()) > 0)
//                continue;
//            processed.insert(read.name());
//            processed.insert(ConjugateContigId(read.name()));
//
//            Untangle(read.sequence(), read.name(), color);
//        }
//    }
//
//    void Untangle(const Sequence& contig, const string& name, TColorSet color) {
//        VERIFY(color == kRedColorSet || color == kBlueColorSet);
//        DEBUG("Untangling contig " << name);
//        Path<EdgeId> path = MapperInstance(old_gp_).MapSequence(contig).path();
//        vector<EdgeId> new_path;
//        DEBUG("Mapped contig" << name);
//        for (size_t i = 0; i < path.size(); i++) {
//            EdgeId next;
//            if (old_coloring_.Color(path[i]) != kVioletColorSet) {
//                DEBUG("Next edge is not purple");
//                size_t j = i;
//                vector<EdgeId> to_glue;
//                while (j < path.size()
//                        && old_coloring_.Color(path[j]) != kVioletColorSet) {
//                    to_glue.push_back(path[j]);
//                    j++;
//                }
//                Sequence new_edge_sequence = MergeSequences(old_gp_.g, to_glue);
//                next = new_gp_.g.AddEdge(GetStartVertex(path, i),
//                        GetEndVertex(path, j - 1), new_edge_sequence);
//                DEBUG(
//                        "Added shortcut edge " << new_gp_.g.int_id(next) << " for path " << old_gp_.g.str(to_glue));
//                i = j - 1;
//            } else {
//                DEBUG("Next edge is purple");
//                next = purple_edge_mapping_[path[i]];
//            }
//            new_path.push_back(next);
//            DEBUG("Coloring new edge and complement");
//            PaintEdgeWithVertices(next, color);
//        }
//        if (color == kRedColorSet) {
//            VERIFY(new_gp_.red_paths.find(name) == new_gp_.red_paths.end());
//            new_gp_.red_paths[name] = new_path;
//            new_gp_.red_paths[ConjugateContigId(name)] = ConjugatePath(
//                    new_gp_.g, new_path);
//        } else {
//            VERIFY(new_gp_.blue_paths.find(name) == new_gp_.blue_paths.end());
//            new_gp_.blue_paths[name] = new_path;
//            new_gp_.blue_paths[name] = ConjugatePath(new_gp_.g, new_path);
//        }
//    }
//
//    vector<EdgeId> ConjugatePath(const Graph& g, const vector<EdgeId> path) {
//        vector<EdgeId> answer;
//        for (int i = path.size() - 1; i >= 0; i--) {
//            answer.push_back(g.conjugate(path[i]));
//        }
//        return answer;
//    }
//
//    template<class T>
//    void ColorWithConjugate(T t, TColorSet color) {
//        new_gp_.coloring.Paint(t, color);
//        new_gp_.coloring.Paint(new_gp_.g.conjugate(t), color);
//    }
//
//    void PaintEdgeWithVertices(EdgeId e, TColorSet color) {
//        DEBUG(
//                "Coloring edges " << new_gp_.g.int_id(e) << " and " << new_gp_.g.int_id(new_gp_.g.conjugate(e)));
//        ColorWithConjugate(e, color);
//        ColorWithConjugate(new_gp_.g.EdgeStart(e), color);
//        ColorWithConjugate(new_gp_.g.EdgeEnd(e), color);
//    }
//
//public:
//    UntangledGraphConstructor(const gp_t &old_gp,
//            const ColorHandler<Graph> &old_coloring,
//            bp_graph_pack<Graph>& new_gp, io::IReader<io::SingleRead> &stream1,
//            io::IReader<io::SingleRead> &stream2) :
//            old_gp_(old_gp), old_coloring_(old_coloring), new_gp_(new_gp) {
//        const Graph& old_graph = old_gp.g;
//        //adding vertices
//        restricted::set<VertexId> processed_purple_v;
//        for (auto it = old_graph.begin(); it != old_graph.end(); ++it) {
//            if (processed_purple_v.count(*it) > 0)
//                continue;
//            processed_purple_v.insert(*it);
//            processed_purple_v.insert(old_graph.conjugate(*it));
//            vertex_mapping_[*it] = new_gp_.g.AddVertex();
//            vertex_mapping_[old_graph.conjugate(*it)] = new_gp_.g.conjugate(
//                    vertex_mapping_[*it]);
//            DEBUG(
//                    "Adding purple vertex " << new_gp_.g.int_id(vertex_mapping_[*it]) << " corresponding to " << old_graph.int_id(*it) << " and conjugates")
//        }
//
//        restricted::set<EdgeId> processed_purple;
//        //propagating purple color to new graph
//        for (auto it = old_graph.SmartEdgeBegin(); !it.IsEnd(); ++it) {
//            if (processed_purple.count(*it) > 0)
//                continue;
//            processed_purple.insert(*it);
//            processed_purple.insert(old_graph.conjugate(*it));
//
//            if (old_coloring.Color(*it) == kVioletColorSet) {
//                EdgeId new_edge = new_gp_.g.AddEdge(
//                        vertex_mapping_[old_graph.EdgeStart(*it)],
//                        vertex_mapping_[old_graph.EdgeEnd(*it)],
//                        old_graph.EdgeNucls(*it));
//                DEBUG(
//                        "Adding purple edge " << new_gp_.g.int_id(new_edge) << " corresponding to " << old_graph.int_id(*it) << " and conjugate")
//                purple_edge_mapping_[*it] = new_edge;
//                purple_edge_mapping_[old_graph.conjugate(*it)] =
//                        new_gp_.g.conjugate(new_edge);
//                PaintEdgeWithVertices(new_edge, kVioletColorSet);
//            }
//        }
//
//        VERIFY(new_gp_.red_paths.empty());
//        VERIFY(new_gp_.blue_paths.empty());
//
//        Untangle(stream1, 0);
//        Untangle(stream2, 1);
//
//        UntangledGraphContigMapper<bp_graph_pack<Graph>> contig_mapper(new_gp_);
//        visualization::position_filler::FillPos(new_gp_.g, contig_mapper, new_gp_.edge_pos, stream1);
//        visualization::position_filler::FillPos(new_gp_.g, contig_mapper, new_gp_.edge_pos, stream2);
//    }
//private:
//    DECL_LOGGER("UntangledGraphConstructor")
//    ;
//};
//
////Currently works for conjugate graphs only
//template<class Graph>
//class RestrictedOneManyResolver {
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//
//    Graph& g_;
//    const ColorHandler<Graph>& coloring_;
//    TColorSet restricting_color_;
//
//    bool CheckColor(const vector<EdgeId>& edges) {
//        DEBUG("Checking color")
//        for (auto it = edges.begin(); it != edges.end(); ++it) {
//            if (coloring_.Color(*it) != restricting_color_) {
//                DEBUG("fail")
//                return false;
//            }
//        }DEBUG("ok")
//        return true;
//    }
//
//    bool CheckColor(VertexId v) {
//        return CheckColor(g_.IncomingEdges(v))
//                && CheckColor(g_.OutgoingEdges(v));
//    }
//
//    bool CheckSimple(const vector<EdgeId>& edges) {
//        DEBUG("Checking simple")
//        for (auto it = edges.begin(); it != edges.end(); ++it) {
//            if (g_.EdgeStart(*it) == g_.EdgeEnd(*it)
//                    || g_.EdgeStart(*it) == g_.conjugate(g_.EdgeEnd(*it))) {
//                DEBUG("fail")
//                return false;
//            }
//        }DEBUG("ok")
//        return true;
//    }
//
//    bool CheckSimple(VertexId v) {
//        return CheckSimple(g_.IncomingEdges(v))
//                && CheckSimple(g_.OutgoingEdges(v));
//    }
//
//    bool CheckVertex(VertexId v) {
//        return CheckSimple(v) && CheckColor(v);
//    }
//
//    void SplitVertex(VertexId v) {
//        DEBUG("Splitting vertex " << g_.str(v))
//        EdgeId incoming_edge = g_.GetUniqueIncomingEdge(v);
//        vector<EdgeId> outgoing_edges = g_.OutgoingEdges(v);
//        DEBUG("Going to create " << outgoing_edges.size() << " new edges")
//        for (auto it = outgoing_edges.begin(); it != outgoing_edges.end();
//                ++it) {
//            VertexId copy_vertex = g_.AddVertex(g_.data(v));
//            EdgeId e1 = g_.AddEdge(g_.EdgeStart(incoming_edge), copy_vertex,
//                    g_.data(incoming_edge));
//            g_.FireProject(incoming_edge, e1);
//            EdgeId e2 = g_.AddEdge(copy_vertex, g_.EdgeEnd(*it), g_.data(*it));
//            //todo think of better way!!! now not stable and awful because of th order of information transfer!!!
//            g_.FireProject(*it, e2);
//            EdgeId e = g_.MergePath(vector<EdgeId> { e1, e2 });
//            DEBUG("Created edge " << g_.str(e))
//        }
//        g_.ForceDeleteVertex(v);
//    }
//
//public:
//    RestrictedOneManyResolver(Graph& g, const ColorHandler<Graph>& coloring,
//            TColorSet restricting_color) :
//            g_(g), coloring_(coloring), restricting_color_(restricting_color) {
//
//    }
//
//    void Resolve() {
//        INFO("Running one-many resolve");
//        for (auto it = g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
//            DEBUG("Checking vertex " << g_.str(*it) << " for split.")
//            if (g_.IncomingEdgeCount(*it) == 1 && CheckVertex(*it)) {
//                DEBUG("Condition was satisfied.")
//                SplitVertex(*it);
//            } else {
//                DEBUG("Condition was not satisfied.")
//            }
//        }INFO("Finished one-many resolve");
//    }
//
//private:
//    DECL_LOGGER("RestrictedOneManyResolver")
//    ;
//};
//
//}
