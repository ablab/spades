#ifndef DISTANCEESTIMATOR_H

#define DISTANCEESTIMATOR_H



#define maxweight 1000
#define middleweight 500
#define variation 20
#define smallEdgeLength 10
#define smallCov 5
#define minEdgeCoverageRate 0.010
/*
 * The interface for Matepair re-adjustment!
 * TODO tdd this class
 */
template<class Graph>
class IDistanceEstimator{
    public:
        virtual ~IDistanceEstimator(){}
        virtual void MatePairFilter(PairInfoIndexData<typename Graph::EdgeId>& pairInfos) = 0;
};


//TODO tdd this class from scratch
template <class Graph>
class SimpleDistanceEstimator: public IDistanceEstimator<Graph>
{
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef omnigraph::PairInfo<EdgeId> PairInfo;

    private:
        Graph &graph_;
        
        DistanceCounter<Graph> distanceTool_;
    public:
        SimpleDistanceEstimator(Graph &graph):graph_(graph), distanceTool_(graph){}
    
        virtual ~SimpleDistanceEstimator(){}
        virtual void MatePairFilter(PairInfoIndexData<typename Graph::EdgeId> &pairInfos) 
        {

            PairInfoIndexData<typename Graph::EdgeId>  *filteredIndex = new PairInfoIndexData<typename Graph::EdgeId>();

            for(auto iter = pairInfos.begin(); iter != pairInfos.end() ; ++iter)
            {
                //if distance < 0 we just ignore, since there should be a reverse pair with possitive distance
                if(iter->d <0 )
                    continue;
                size_t lengthFirstEdge = graph_.length(iter->first);
                size_t lengthSecondEdge = graph_.length(iter->second);
                // if the weight is too small, ignore it
                if((min (lengthFirstEdge, lengthSecondEdge) > smallEdgeLength) && (iter->weight < smallCov)){ 
                    continue;
                }
                if(min(lengthFirstEdge, lengthSecondEdge)* minEdgeCoverageRate > iter->weight)
                {
                    continue;
                }

                // if the same edges, and the distance estimated is shorter than the length of the first edge, it should also be zero.
                // one problem occurs when the sequence traverses the loop multiple time, but this case is rare and the 
                // effect of this error is not substantial
                if(iter->first == iter->second)
                {

                    if(graph_.length(iter->first) > iter->d){
                        PairInfo *pairInfo = new PairInfo(iter->first, iter->first, 0, maxweight, 0.);
                        filteredIndex->AddPairInfo(*pairInfo,1);
                    }
                    else
                        filteredIndex->AddPairInfo(*iter);
                }
                //If they are adjacent edges
                else if(graph_.EdgeStart(iter->second) == graph_.EdgeEnd(iter->first))
                {
                    if( iter->d  < lengthFirstEdge + lengthSecondEdge )  {
                        PairInfo *pairInfo = new PairInfo(iter->first, iter->second, lengthFirstEdge, maxweight, 0.);
                        filteredIndex->AddPairInfo(*pairInfo,1);
                    }
                    else
                        filteredIndex->AddPairInfo(*iter);

                }
                else 
                {
                    //in general, if the distance is smaller than the shortest distance between two vertices, it should be changed
                    //to at least the min distance
                    VertexId endOfFirstEdge = graph_.EdgeEnd(iter->first); 
                    VertexId startOfSecondEdge = graph_.EdgeStart(iter->second); 
                    size_t graphDistance;
                    //If they are shorter than the shortest distance
                    if((distanceTool_.IsReachable(endOfFirstEdge, startOfSecondEdge)) &&    ((graphDistance = distanceTool_.Distance(endOfFirstEdge, startOfSecondEdge)) + graph_.length(iter->first) > iter->d) )
                    {
                        PairInfo *pairInfo = new PairInfo(iter->first, iter->second, graphDistance + graph_.length(iter->first) , middleweight, 0.);
                        filteredIndex->AddPairInfo(*pairInfo,1);
                    }
                    //If they are not reachable
                    else if (!distanceTool_.IsReachable(endOfFirstEdge, startOfSecondEdge))
                    {
                        continue;
                    }
                    //In other cases, just put it it
                    else{
                        filteredIndex->AddPairInfo(*iter);
                    }
                }
            }
            pairInfos.clear();
            pairInfos  = *filteredIndex;

        }

};



#endif /* end of include guard: DISTANCEESTIMATOR_H */

