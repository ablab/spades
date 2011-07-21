#ifndef GRAPHCOPY_H

#define GRAPHCOPY_H
#include <cmath>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>

#include "logging.hpp"
#include "paired_info.hpp"
#include "config.hpp"
#include "omni_utils.hpp"

#include "omni_tools.hpp"
#include "omnigraph.hpp"

#include "ID_track_handler.hpp"
#include "dijkstra.hpp"


using omnigraph::SmartVertexIterator;
using omnigraph::Compressor;
using omnigraph::PairedInfoIndex;
using omnigraph::PairInfoIndexData;


namespace debruijn_graph{

template<class Graph>
class GraphCopy {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    Graph &graph_;
public:
    GraphCopy(Graph &graph) :
        graph_(graph) {
        }
    template<class CopyGraph>
        void Copy(CopyGraph& new_graph, PairInfoIndexData<EdgeId> &old_index, PairInfoIndexData<typename CopyGraph::EdgeId> &new_index)
        {
            typedef typename CopyGraph::EdgeId NewEdgeId;
            typedef typename CopyGraph::VertexId NewVertexId;



            map<VertexId, NewVertexId> copy;
            map<EdgeId, NewEdgeId> edgeMap;
            //Copy Nodes
            for (auto iter = graph_.begin(); iter != graph_.end(); ++iter) {
                NewVertexId new_vertex = new_graph.AddVertex(graph_.data(*iter));
                copy.insert(make_pair(*iter, new_vertex));
            }
            //Copy Edges
            for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
                EdgeId edge = *iter;
                NewEdgeId newEdge = new_graph.AddEdge(copy[graph_.EdgeStart(edge)],
                        copy[graph_.EdgeEnd(edge)], graph_.data(edge));
                edgeMap.insert(make_pair(edge, newEdge));
            }
            //Copy PairInfoIndexData
            for(auto iter = old_index.begin() ; iter != old_index.end() ; ++iter)
            {
                PairInfo<NewEdgeId> *newPairInfo = new PairInfo<NewEdgeId>(edgeMap[iter->first] , edgeMap[iter->second], iter->d, iter->weight);
                new_index.AddPairInfo(*newPairInfo,1);
            }
        }




    void Copy(Graph& new_graph)
    {
        map<VertexId, VertexId> copy;
        map<EdgeId, EdgeId> edgeMap;
        //Copy Nodes
        for (auto iter = graph_.begin(); iter != graph_.end(); ++iter) {
            VertexId new_vertex = new_graph.AddVertex(graph_.data(*iter));
            copy.insert(make_pair(*iter, new_vertex));
        }
        //Copy Edges
        for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            EdgeId edge = *iter;
            new_graph.AddEdge(copy[graph_.EdgeStart(edge)],
                    copy[graph_.EdgeEnd(edge)], graph_.data(edge));
        }
    }
};


}
#endif /* end of include guard: GRAPHCOPY_H */

