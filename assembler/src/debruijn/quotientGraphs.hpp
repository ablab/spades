#ifndef QUOTIENTGRAPHS_H
#define QUOTIENTGRAPHS_H
#include "omni/dijkstra.hpp"
#include "omni/omni_utils.hpp"


namespace debruijn_graph {
template<class Graph>
class IQuotientGraphs{
    private: 
        typedef typename Graph::EdgeId EdgeId;
        typedef typename Graph::VertexId VertexId;
    public:
        virtual ~IQuotientGraphs()
        {
        }
        /*
         * Add a new quotient graph to the store
         */
        virtual bool AddGraph(const Graph &graph, std::multimap<EdgeId, EdgeId>  &edgeMaps) =0;
        /* Check if two edges are \emph{adjacent} in 
         * all the quotient graphs. Note that edge here should belongs to the first 
         * quotient graph added. By Adjacent I mean edge2 follows edge1 
         */
        virtual bool IsAdjacent(const EdgeId& edge1,const EdgeId& edge2) const = 0;
        /*
         * Return the shortest distance between two edges. (A path that passes through fromEdge and reach toEdge
         */
        virtual size_t Distance(EdgeId & fromEdge, EdgeId &toEdge) = 0;
};
/*
 * Here, just one quotient graph of the original sequence is supported
 */
template<class Graph>
class SingleQuotientGraph: public IQuotientGraphs<Graph>{
    Graph &graph_; 
    DistanceCounter<Graph> distanceTool_;
    public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    SingleQuotientGraph(Graph& graph): graph_(graph), distanceTool_(graph) {}
    //Not support for singleQuotientGraph
    bool AddGraph(const Graph &graph, std::multimap<EdgeId,EdgeId> &edgeMap){return true;}
    bool IsAdjacent(const EdgeId& edge1, const EdgeId& edge2) const
    {
        if (graph_.EdgeEnd(edge1) == graph_.EdgeStart(edge2)) {
                    return true;
        }
        return false;
    }
    //TODO
    size_t Distance(EdgeId &fromEdge, EdgeId  &toEdge) 
    {
        VertexId from = graph_.EdgeEnd(fromEdge);
        VertexId to = graph_.EdgeStart(toEdge);
        if (distanceTool_.IsReachable(from,to)) 
        {
            return distanceTool_.Distance(from,to); 
        }
        return SIZE_MAX;
    }
};

template<class Graph>
class MultipleQuotientGraphs: public IQuotientGraphs<Graph>{
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    Graph &graph_; 
    DistanceCounter<Graph> distanceTool_;
    vector<Graph&> otherGraphs_;
    vector< DistanceCounter<Graph>* > otherDistanceTools_; 
    vector< multimap<EdgeId, EdgeId> &> edgeMaps_;
    public:

    MultipleQuotientGraphs(Graph& graph): graph_(graph),distanceTool_(graph){}
    //Not support for singleQuotientGraph
    bool AddGraph(const Graph &graph, multimap<EdgeId,EdgeId> &edgeMap){
        otherGraphs_.push_back(graph);
        DistanceCounter<Graph> newDistanceTool = new DistanceCounter<Graph>(graph);
        otherDistanceTools_.push_back(newDistanceTool);
        edgeMaps_.push_back(edgeMap);
        return true;

    }

    /* If for every quotient graph, there exists at least 
     * one instance of pair of edges that are adjacent, return true
     * else, return false
     */
    bool IsAdjacent(const EdgeId& edge1, const EdgeId& edge2) const
    {
        if (graph_.EdgeEnd(edge1) != graph_.EdgeStart(edge2)) {
                    return false;
        }
        for(size_t i = 0 ; i < otherGraphs_.size() ; ++i)
        {
            Graph& currentGraph = otherGraphs_[i];
            typename std::multimap<EdgeId,EdgeId> &currentEdgeMap = edgeMaps_[i];
            typename std::multimap<EdgeId,EdgeId>::iterator it, itlow, itup;
            itlow = currentEdgeMap.lower_bound(edge1);
            itup = currentEdgeMap.upper_bound(edge1);
            restricted::set<VertexId> firstEnds;
            for(it = itlow ; it != itup ; ++it)
            {
                firstEnds.insert(currentGraph.EdgeEnd( it->second)); 
            }
            int originalSize = firstEnds.size();
            int updatedSize = originalSize;
            itlow = currentEdgeMap.lower_bound(edge2);
            itup = currentEdgeMap.upper_bound(edge2);
            for(it = itlow; it!= itup ; ++it)
            {
                if (firstEnds.find(it->second) != firstEnds.end() ) {
                            updatedSize --;
                }
            }
            if (updatedSize == originalSize) {
                        return false;
            }
        }

        return true;
    }
    //TODO
    //in all the quotient graphs.
    size_t Distance(EdgeId &fromEdge, EdgeId  &toEdge) 
    {
        VertexId from = graph_.EdgeEnd(fromEdge);
        VertexId to = graph_.EdgeStart(toEdge);
        if (distanceTool_.IsReachable(from,to)) 
        {
            return distanceTool_.Distance(from,to); 
        }
        return SIZE_MAX;
    }
};

} /* debruijn */
#endif /* end of include guard: QUOTIENTGRAPHS_H */
