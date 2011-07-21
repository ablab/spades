#ifndef RECTANGLEREPEATRESOLVER_H

#define RECTANGLEREPEATRESOLVER_H

#include "weakerGluer.hpp"
#include "graphCopy.hpp"
#include "distanceEstimator.hpp"
namespace debruijn_graph {


template<class Graph>
class RectangleRepeatResolver {
private:
    Graph &graph_;
    Graph &resultGraph_;
    PairInfoIndexData<typename Graph::EdgeId> &index_;
    size_t errorDistance_;
public:

    RectangleRepeatResolver(Graph &graph, PairInfoIndexData<typename Graph::EdgeId> &index, Graph &resultGraph, size_t errorDistance) :
         graph_(graph), resultGraph_(resultGraph), index_(index), errorDistance_(errorDistance) {}

    void Process()
    {
        NonconjugateDeBruijnGraph myRed(graph_.k());
        NonconjugateDeBruijnGraph myBlue(graph_.k());

        SimpleDistanceEstimator<Graph>  simpleDistanceEstimator(graph_);
        SingleQuotientGraph<NonconjugateDeBruijnGraph> myquotientGraph(graph_);
        SimpleGraphFilter<Graph> simpleGraphFilter;
        WeakerGluer<NonconjugateDeBruijnGraph> weakerGluer(graph_, myquotientGraph, simpleGraphFilter, 
                simpleDistanceEstimator,  index_, errorDistance_, myRed, myBlue );
        weakerGluer.WeakerGlueProcess();

        //TODO combine red and blue graph for better contigs.
        GraphCopy<Graph> resultCopy(myRed);
        resultCopy.Copy(resultGraph_);
    }
};

}
#endif /* end of include guard: RECTANGLEREPEATRESOLVER_H */
