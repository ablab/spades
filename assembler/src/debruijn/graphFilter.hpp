#ifndef GRAPHFILTER_H

#define GRAPHFILTER_H

//This graph will be a crucial component in multi-library
//TODO tdd this in the next iteration
template<class Graph>
class IGraphFilter
{
public:
    virtual ~IGraphFilter(){};
    /*
     * Filter all \weakly non eulerian paths.
     */
    virtual void Filter(Graph& graph) = 0;
};


template <class Graph>
class SimpleGraphFilter: public IGraphFilter<Graph>
{
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
public:
    SimpleGraphFilter(){}
    void Filter(Graph& graph) ;
    ~SimpleGraphFilter(){};
};

template<class Graph>
void SimpleGraphFilter<Graph>::Filter(Graph& graph)
{
    for (auto iter = graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
        EdgeId edge = *iter;
        VertexId start = graph.EdgeStart(edge);
        VertexId end = graph.EdgeEnd(edge);
        if((graph.OutgoingEdgeCount(start) == 1) && (graph.IncomingEdgeCount(start) == 0) && (graph.IncomingEdgeCount(end) == 1) && (graph.OutgoingEdgeCount(end) == 0))
        {
            graph.DeleteEdge(edge);
            graph.DeleteVertex(start);
            graph.DeleteVertex(end);
        }
    }
}
#endif /* end of include guard: GRAPHFILTER_H */
