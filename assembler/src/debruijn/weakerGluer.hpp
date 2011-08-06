#ifndef WEAKERGLUER_H

#define WEAKERGLUER_H

#include <cmath>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include "common/io/single_read.hpp"
#include "logging.hpp"
#include "config.hpp"
#include "omni_utils.hpp"
#include "omni_tools.hpp"
#include "omnigraph.hpp"
#include "ID_track_handler.hpp"
#include "union.h"   
#include "paired_info.hpp"
#include "new_debruijn.hpp"
#include "quotientGraphs.hpp"
#include "dijkstra.hpp"
#include "distanceEstimator.hpp"
#include "graphFilter.hpp"
#include "graphCopy.hpp"
#include <time.h>
namespace debruijn_graph {


using omnigraph::SmartVertexIterator;
using omnigraph::Compressor;
using omnigraph::PairedInfoIndex;
using omnigraph::PairInfoIndexData;
using debruijn_graph::NonconjugateDeBruijnGraph;



enum Direction{FORWARD, BACKWARD};



template<class Graph>
class WeakerGluer {
public:



	//	typedef SmartVertexIterator<Graph> VertexIter;
	//	typedef omnigraph::SmartEdgeIterator<Graph> EdgeIter;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef omnigraph::PairInfo<EdgeId> PairInfo;
	typedef omnigraph::PairInfoIndexData<typename Graph::EdgeId> PIIndex;
	typedef vector<PairInfo> PairInfos;
    

	typedef map<VertexId, set<EdgeId> > NewVertexMap;
	typedef map<VertexId, set<VertexId> > VertexIdMap;
    typedef map<PairInfo, size_t> PairInfoVertexMap;


    WeakerGluer(Graph &debruijn, IQuotientGraphs<Graph> &quotientGraphs, IGraphFilter<Graph> &graphFilter, IDistanceEstimator<Graph> &distanceEstimator, 
            PairInfoIndexData<typename Graph::EdgeId> &index, int errorDistance,
             Graph &weakerRedGraph, Graph &weakerBlueGraph)
        :debruijn_(debruijn), quotientGraphs_(quotientGraphs), graphFilter_(graphFilter), distanceEstimator_(distanceEstimator),
        pairIndex_(index), errorDistance_(errorDistance), weakerRedGraph_(weakerRedGraph),
        weakerBlueGraph_(weakerBlueGraph) {
                
            


        //Want to print all the pairinfos
       distanceEstimator_.MatePairFilter(pairIndex_);
    }
    WeakerGluer(){}
    bool WeakerGlueProcess()
    {
        //TODO Iterate until converge
        Graph redGraph(debruijn_.k());
        Graph blueGraph(debruijn_.k());
        //Clear all intermediateData in the previous iteration
        ClearIntermediateData();
        //Generate the red graph
        GenerateRedGraph(redGraph);
        //Generate the blue graph
        GenerateBlueGraph(blueGraph);
        //add red and blue graphs to the quotient store to use later 
        quotientGraphs_.AddGraph(redGraph,redGraphEdgeMaps_);
        quotientGraphs_.AddGraph(blueGraph, blueGraphEdgeMaps_);
        //Filtered the graph because of noise in distance
        graphFilter_.Filter(redGraph, debruijn_);
        graphFilter_.Filter(blueGraph, debruijn_);
        //store red and blue graphs so that the caller can get the information
        GenerateResultGraphs(redGraph, blueGraph);
        return true;
    }
private:
    Graph &debruijn_; 
    IQuotientGraphs<Graph> &quotientGraphs_;
    IGraphFilter<Graph> &graphFilter_;
    IDistanceEstimator<Graph> &distanceEstimator_;
    PIIndex &pairIndex_;
    int errorDistance_;
    Graph &weakerRedGraph_;
    Graph &weakerBlueGraph_;
    multimap<EdgeId, EdgeId> redGraphEdgeMaps_;
    multimap<EdgeId, EdgeId> blueGraphEdgeMaps_;
    

    void GenerateRedGraph(Graph& redGraph);
    void GenerateBlueGraph(Graph& blueGraph);
    void GenerateResultGraphs(Graph& redGraph, Graph& blueGraph);
    void ClearIntermediateData();
    size_t Glue(vector<PairInfo> &leftPairInfos, vector<PairInfo> &rightPairInfos, Direction direction,
        map<PairInfo, size_t>& leftPairInfosToNewNodesID,map<PairInfo, size_t>& rightPairInfosToNewNodesID,size_t maxNodeId);
     void OneSideGluing(vector<PairInfo> &pairInfos, vector<pair<int,int> > &links, Direction direction); 
     void LeftToRightGlue(vector<PairInfo> &leftPairInfos, vector<PairInfo> &rightPairInfos, vector<pair<int,int> > &leftToRightGluings, Direction direction);
     void ChooseMajorEdges(const PairInfo &pairInfo, const Direction direction, EdgeId& majorEdge, EdgeId& minorEdge);

    /*
     * We present the simplest algorithm to find the 
     * distance between 2 edges
     * 1. (e_1, e_1, d, w): if 0 < d < length(e_1): replace by (e_1, e_1,0, maxweight)
     * others will be added as moving on
     */

    void GetPairInfos(const vector<EdgeId> &edges, Direction direction, vector<PairInfo> &pairInfos);
};


template<class Graph>
void WeakerGluer<Graph>::ChooseMajorEdges(const PairInfo &pairInfo, const Direction direction, EdgeId& majorEdge, EdgeId& minorEdge)
{
    if(direction == FORWARD)
        {
            majorEdge = pairInfo.first;  
            minorEdge = pairInfo.second;
        }
        else
        {
            majorEdge = pairInfo.second ;
            minorEdge = pairInfo.first;
        }
}

/*
 * Here we glue all rectangles (or adjacent pairinfos - in the context of Dima and Shurik) on 
 * that form a vertical stack.
 */
//TODO Check the correctness for the reverse direction.
template <class Graph>
void WeakerGluer<Graph>::OneSideGluing(vector<PairInfo> &pairInfos, vector<pair<int,int> > &links, Direction direction){
    EdgeId firstMajorEdge, secondMajorEdge;
    EdgeId firstMinorEdge, secondMinorEdge;
    int firstDistance, secondDistance;
   //long code, but legible
    for(size_t i = 0 ;  i < pairInfos.size() ; ++i)
    {
        ChooseMajorEdges(pairInfos[i], direction, firstMajorEdge, firstMinorEdge);
        firstDistance = abs(pairInfos[i].d);
        for(size_t j = i+1 ; j < pairInfos.size() ; ++j)
        {
            ChooseMajorEdges(pairInfos[j], direction, secondMajorEdge, secondMinorEdge);
            secondDistance = abs(pairInfos[j].d);
            if(firstMajorEdge != secondMajorEdge)
                continue;
            if((firstMinorEdge == secondMinorEdge) &&  abs( firstDistance - secondDistance ) < errorDistance_)
            {
                links.push_back(make_pair(i,j));
                continue;
            }
            int firstMinorEdgeLength = debruijn_.length(firstMinorEdge);
            int secondMinorEdgeLength = debruijn_.length(secondMinorEdge);
            int varDistanceFirstSecond = quotientGraphs_.Distance(firstMinorEdge, secondMinorEdge) ;
            if((varDistanceFirstSecond < errorDistance_)  && ( abs(firstDistance + firstMinorEdgeLength  - (secondDistance - varDistanceFirstSecond)) < errorDistance_  ))
            {
                links.push_back(make_pair(i,j));
                continue;
            }
            int  varDistanceSecondFirst = quotientGraphs_.Distance(secondMinorEdge, firstMinorEdge);
            if( ( varDistanceSecondFirst < errorDistance_ ) && ( abs( secondDistance + secondMinorEdgeLength - firstDistance + varDistanceSecondFirst ) < errorDistance_ ))
            {
                links.push_back(make_pair(i,j));
                continue;
            }
        }
    }
}
/*
 * We glue all side adjacent rectangles
 */
template <class Graph>
void WeakerGluer<Graph>::LeftToRightGlue(vector<PairInfo> &leftPairInfos, vector<PairInfo> &rightPairInfos, vector<pair<int,int> > &leftToRightGluings, Direction direction){
    
    EdgeId firstMajorEdge, firstMinorEdge; 
    EdgeId secondMajorEdge, secondMinorEdge;
    int firstDistance, secondDistance;
    int firstMajorEdgeLength , firstMinorEdgeLength;
    int  secondMinorEdgeLength;
    for(size_t i = 0 ; i < leftPairInfos.size() ; ++i)
    {
        ChooseMajorEdges(leftPairInfos[i], direction, firstMajorEdge,firstMinorEdge );
        firstDistance = abs(leftPairInfos[i].d); 
        firstMajorEdgeLength = debruijn_.length(firstMajorEdge);
        for(size_t j = 0  ; j < rightPairInfos.size() ; ++j)
        {
            secondDistance = abs(rightPairInfos[j].d);
            secondDistance = rightPairInfos[j].d;
            ChooseMajorEdges(rightPairInfos[j], direction, secondMajorEdge, secondMinorEdge);

            /*
            if(firstMajorEdge == secondMajorEdge)
            {
                continue;
            }
            */
            if(!quotientGraphs_.IsAdjacent(firstMajorEdge, secondMajorEdge))
            {
                continue;
            }
            if( (firstMinorEdge == secondMinorEdge) && ( abs( firstDistance  -  firstMajorEdgeLength  - secondDistance  ) < errorDistance_ ) ) 
            {
                leftToRightGluings.push_back(make_pair(i, j + leftPairInfos.size() ));
                continue;
            }
             int distanceFirstSecond = quotientGraphs_.Distance(firstMinorEdge, secondMinorEdge);
             firstMinorEdgeLength = debruijn_.length(firstMinorEdge);
            if((distanceFirstSecond < errorDistance_) && ( abs(firstDistance - firstMajorEdgeLength + firstMinorEdgeLength - (secondDistance - distanceFirstSecond) ) < errorDistance_)   )
            {
                leftToRightGluings.push_back(make_pair(i, j + leftPairInfos.size()));
                continue;
            }
            int distanceSecondFirst = quotientGraphs_.Distance(secondMinorEdge, firstMinorEdge);
            secondMinorEdgeLength = debruijn_.length(secondMinorEdge);

            if( (distanceSecondFirst < errorDistance_) && ( abs(secondDistance + secondMinorEdgeLength - (firstDistance - distanceSecondFirst - firstMajorEdgeLength )  )  < errorDistance_  ) )
            {
                leftToRightGluings.push_back(make_pair(i, j + leftPairInfos.size()));
                continue;
            }
        }
    }
}
template <class Graph>
void WeakerGluer<Graph>::GenerateRedGraph(Graph &redGraph){
    //The debruijn graph is a result of a strong glue
    //We define a weaker glue of each red nodes.
    map<PairInfo,size_t> leftPairInfosToNewNodesID;// we have to store this map to add the outgoing edges back to the graph
    map<PairInfo,size_t> rightPairInfosToNewNodesID;//we have to store this map to add the incoming edges back to the graph
    map<size_t, VertexId> nodeIntToVertexId;
    size_t maxNodeId = 0;   
    for (auto iter = debruijn_.SmartVertexBegin(); !iter.IsEnd(); ++iter) 
    {
        vector<EdgeId> outGoingEdges = debruijn_.OutgoingEdges(*iter);
        vector<EdgeId> inComingEdges = debruijn_.IncomingEdges(*iter);
        vector<PairInfo> leftPairInfos ;
        GetPairInfos(inComingEdges, FORWARD, leftPairInfos);
        vector<PairInfo> rightPairInfos;
        GetPairInfos(outGoingEdges,FORWARD, rightPairInfos);
        size_t currentMaxNodeId = Glue(leftPairInfos, rightPairInfos, FORWARD, leftPairInfosToNewNodesID, rightPairInfosToNewNodesID, maxNodeId);

        for(size_t currentNodeId = maxNodeId ; currentNodeId != currentMaxNodeId ; ++currentNodeId)
        {
            VertexId id = redGraph.AddVertex();
            nodeIntToVertexId.insert(make_pair(currentNodeId,id));
        }
        maxNodeId = currentMaxNodeId;
    }
    multimap< pair<VertexId, VertexId> , string> edgesMaps;//get rid of parallel edges in a dirty way.Actually this is incorrect. Should use EdgeNucls instead.
    for(auto iter = leftPairInfosToNewNodesID.begin(); iter != leftPairInfosToNewNodesID.end(); ++iter)
    {
        PairInfo pairInfo = iter->first;
        VertexId endVertex = nodeIntToVertexId[iter->second];
        VertexId startVertex = nodeIntToVertexId[rightPairInfosToNewNodesID[pairInfo]];
        auto range = edgesMaps.equal_range(make_pair(startVertex, endVertex));
        bool isAlreadyAdded = false;
        for(auto iter = range.first ; iter != range.second ; ++iter)
        {
            if(iter->second == debruijn_.str(pairInfo.first))
            {
                isAlreadyAdded = true;
                break;
            }

        }
        if(isAlreadyAdded)
            continue;
        edgesMaps.insert(make_pair( make_pair(startVertex, endVertex ), (debruijn_.str(pairInfo.first) ) ));
        redGraph.AddEdge(startVertex, endVertex, debruijn_.data(pairInfo.first));
    }
}


template<class Graph>
size_t WeakerGluer<Graph>::Glue(vector<PairInfo> &leftPairInfos, vector<PairInfo> &rightPairInfos, Direction direction,
        map<PairInfo, size_t>& leftPairInfosToNewNodesID,map<PairInfo, size_t>& rightPairInfosToNewNodesID, size_t maxNodeId)
{

     vector<pair<int,int> > leftGluings; 
     OneSideGluing(leftPairInfos, leftGluings, direction); 
     vector<pair<int,int> > rightGluings ; 
     OneSideGluing(rightPairInfos, rightGluings, direction);
     vector<pair<int,int> > leftToRightGluings;
     LeftToRightGlue(leftPairInfos, rightPairInfos, leftToRightGluings, direction);
     size_t leftSize = leftPairInfos.size();
     size_t rightSize = rightPairInfos.size();
     UnionFindClass myunion(leftSize + rightSize);
     for(size_t  i =  0  ; i < leftGluings.size() ; ++i)
     {
         myunion.unionn(leftGluings[i].first, leftGluings[i].second);
     }
     for(size_t i = 0 ; i < rightGluings.size() ; ++i)
     {
         myunion.unionn(rightGluings[i].first  + leftSize, rightGluings[i].second + leftSize);
     }
     for(size_t i = 0 ; i < leftToRightGluings.size() ; ++i)
     {
         myunion.unionn(leftToRightGluings[i].first, leftToRightGluings[i].second);
     }
     vector<vector<int> > resultsGluing ;
     myunion.get_classes(resultsGluing);
     size_t nodeID= maxNodeId;
     for(size_t i = 0 ; i < resultsGluing.size() ; i++)
     {
         nodeID += i ;
         for(auto iter = resultsGluing[i].begin() ; iter !=  resultsGluing[i].end() ; ++iter)
         {
             if((size_t)*iter < leftSize) //TODO refractor the Union
             {
                 leftPairInfosToNewNodesID.insert(make_pair(leftPairInfos[*iter],nodeID));
             }
             else
             {
                 rightPairInfosToNewNodesID.insert(make_pair(rightPairInfos[(*iter) - leftPairInfos.size()], nodeID));
             }
         }
     }
     if(resultsGluing.size() != 0)
         return nodeID+1;
     return nodeID;
}
/*
 * Clear all  intermediate variables in the WeakerGluer
 */
template<class Graph>
void WeakerGluer<Graph>::ClearIntermediateData()
{
}
/*
 * Get a vector of pairinfos for a corresponding set of edges.
 * If direction == FORWARD, this is for the red graph
 * If direction == BACKWARD, this is for the blue graph
 * pairInfos is a ref type for output
 */
template<class Graph>
void WeakerGluer<Graph>::GetPairInfos(const vector<EdgeId> &edges, Direction direction, vector<PairInfo> &pairInfos)
{
    for(auto e_iter = edges.begin() ; e_iter != edges.end() ; ++e_iter)
    {
        vector<PairInfo> rawPairs =  pairIndex_.GetEdgeInfos(*e_iter);
        for(auto iter = rawPairs.begin() ; iter != rawPairs.end() ; ++iter)
        {
            if(direction == FORWARD)
            {
                if(iter->d >= 0)
                    pairInfos.push_back(*iter);
            }
            else
                if (iter-> d <= 0 ) 
                    pairInfos.push_back(*iter);
        }
    }
}
template <class Graph>
void WeakerGluer<Graph>::GenerateBlueGraph(Graph &redGraph){}

template <class Graph>
void WeakerGluer<Graph>::GenerateResultGraphs(Graph& redGraph, Graph& blueGraph){

    GraphCopy<Graph> redCopier(redGraph);
    redCopier.Copy(weakerRedGraph_);
    GraphCopy<Graph> blueCopier(blueGraph);
    blueCopier.Copy(weakerBlueGraph_);
        return;
    }

    /* debruijn_graph */
}
#endif 
