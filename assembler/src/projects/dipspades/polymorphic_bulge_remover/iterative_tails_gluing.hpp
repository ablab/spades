//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../utils/bulge_utils.hpp"
#include "../dipspades_config.hpp"

using namespace debruijn_graph;

namespace dipspades {

class IterativeTailGluing {
    Graph &graph_;
    double rel_align_;

    typedef VertexId::type::edge_const_iterator edge_const_iterator;
    typedef pair<edge_const_iterator, edge_const_iterator> edge_iters;
    typedef boost::optional<EdgeId> OptEdgeId;


    OptEdgeId GetEdgeForGlue(EdgeId edge, edge_iters iters){
        double best_rel_align = 1;
        OptEdgeId res;
        for(auto it = iters.first; it != iters.second; it++){
            if(edge != *it){
                pair<Sequence, Sequence> seqs;
                if(graph_.length(edge) <= graph_.length(*it)){
                    seqs.first = graph_.EdgeNucls(edge);
                    seqs.second = graph_.EdgeNucls(*it).Subseq(0, seqs.first.size());
                }
                else{
                    seqs.first = graph_.EdgeNucls(*it);
                    seqs.second = graph_.EdgeNucls(edge).Subseq(0, seqs.first.size());
                }
                double rel_align = RelAlignmentOfSequences(seqs.first, seqs.second);
                if(rel_align <= rel_align_ && rel_align <= best_rel_align){
                    best_rel_align = rel_align;
                    res = *it;
                }
            }
        }
        return res;
    }

    bool ProcessTail(EdgeId edge, edge_iters iters){
        auto edge_for_glue = GetEdgeForGlue(edge, iters);
        if(edge_for_glue.is_initialized()){

            TRACE("Edge for glue " << graph_.str(edge_for_glue.get()));
            TRACE("Edges lengths" << graph_.length(edge) << " - " << graph_.length(edge_for_glue.get()));

            size_t min_len = min<size_t>(graph_.length(edge), graph_.length(edge_for_glue.get()));
            if(min_len == graph_.length(edge) && min_len == graph_.length(edge_for_glue.get())){
                graph_.GlueEdges(edge, edge_for_glue.get());
            }
            else{
                if(min_len == graph_.length(edge)){
                    pair<EdgeId, EdgeId> new_edges = graph_.SplitEdge(edge_for_glue.get(), min_len);
                    graph_.GlueEdges(edge, new_edges.first);
                }
                else {
                    auto new_edges = graph_.SplitEdge(edge, min_len);
                    graph_.GlueEdges(new_edges.first, edge_for_glue.get());
                }
            }
            return true;
        }
        return false;
    }

    bool IsTailIncoming(EdgeId edge){
        return graph_.IncomingEdgeCount(graph_.EdgeStart(edge)) == 0 &&
                graph_.OutgoingEdgeCount(graph_.EdgeStart(edge)) == 0;
    }

    bool ProcessTail(EdgeId edge){
        if(IsTailIncoming(edge))
            return ProcessTail(edge,
                    edge_iters(graph_.IncomingEdges(graph_.EdgeEnd(edge)).begin(),
                            graph_.IncomingEdges(graph_.EdgeEnd(edge)).end()));
        return ProcessTail(edge,
                edge_iters(graph_.OutgoingEdges(graph_.EdgeStart(edge)).begin(),
                        graph_.OutgoingEdges(graph_.EdgeStart(edge)).end()));
    }

    bool EdgeIsTail(EdgeId edge) {
        return (graph_.IncomingEdgeCount(graph_.EdgeStart(edge)) == 0 &&
                graph_.OutgoingEdgeCount(graph_.EdgeStart(edge)) == 1) ||
                        (graph_.IncomingEdgeCount(graph_.EdgeEnd(edge)) == 1 &&
                                        graph_.OutgoingEdgeCount(graph_.EdgeEnd(edge)) == 0);
    }

    bool EdgeIsIsolate(EdgeId edge){
        return (graph_.IncomingEdgeCount(graph_.EdgeStart(edge)) == 0 &&
                graph_.OutgoingEdgeCount(graph_.EdgeStart(edge)) == 1) &&
                        (graph_.IncomingEdgeCount(graph_.EdgeEnd(edge)) == 1 &&
                                        graph_.OutgoingEdgeCount(graph_.EdgeEnd(edge)) == 0);    }

    size_t ProcessTails(){
        size_t num_glued_tails = 0;
        for(auto edge = graph_.SmartEdgeBegin(); !edge.IsEnd(); ++edge)
            if(EdgeIsTail(*edge) && !EdgeIsIsolate(*edge)){
                TRACE("Processing edge " << graph_.str(*edge));
                if(ProcessTail(*edge))
                    num_glued_tails++;
            }
        return num_glued_tails;
    }

public:
    IterativeTailGluing(Graph &graph, double rel_align) :
        graph_(graph),
        rel_align_(rel_align) { }

    size_t IterativeProcessTails(){
        size_t num_glued_tails = 1;
        size_t num_iter = 1;
        while(num_glued_tails > 0){
            num_glued_tails = ProcessTails();
            INFO(num_iter << " iteration : " << num_glued_tails << " tails were glued");
            num_iter++;
        }
        return num_glued_tails;
    }
};

}
