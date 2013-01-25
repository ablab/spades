//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * gap_jumper.hpp
 *
 *  Created on: Jul 17, 2012
 *      Author: tanunia
 */

#ifndef PATH_GAPJUMPER_HPP_
#define PATH_GAPJUMPER_HPP_

#include "bidirectional_path.hpp"
#include <cstdlib>

using namespace path_extend;

class GapJumper{
	Graph& g_;
	WeightCounter* wc_;
	size_t gapLength_;
	double threshold_;

public:
    GapJumper(Graph& g, WeightCounter* wc, size_t gapLength = 0,double threshold = 0.0):
    	g_(g), wc_(wc), gapLength_(gapLength), threshold_(threshold) {
    }

    double CountWeight(BidirectionalPath* left, BidirectionalPath* right)
    {
    	double w = 0.0;
        for (size_t i = 0; i < right->Size(); ++ i) {
            w += wc_->CountWeight(*left,right->At(i));
        }
        return w;
    }

    void CreatePathFrom(BidirectionalPath * left, BidirectionalPath* right, BidirectionalPath& result)
    {
        for (size_t i = 0; i < left->Size(); ++i) {
        	result.PushBack(left->At(i));
        }

        result.PushBack(right->At(0), gapLength_);

        for (size_t i = 1; i < right->Size(); ++i) {
			result.PushBack(right->At(i));
		}
    }

    bool isSource(EdgeId e)
    {
        auto count = g_.IncomingEdgeCount(g_.EdgeStart(e));
    	return !count;
    }

    bool isSink(EdgeId e)
    {
    	auto count = g_.OutgoingEdgeCount(g_.EdgeEnd(e));
    	return !count;
    }

    void Jumps(PathContainer& paths)
    {
    	std::vector<BidirectionalPath*> newPaths;
    	int count = 0;
    	for (size_t i = 0; i < paths.size(); ++ i) {
			 for (size_t j = 0; j < paths.size(); ++ j) {
			      if (isSink(paths.Get(i)->At(paths.Get(i)->Size() - 1)) && isSource(paths.Get(j)->At(0))){
					  double weigth = CountWeight(paths.Get(i), paths.Get(j));
					  if (weigth > threshold_){
						 count ++;
						 BidirectionalPath newPath(g_);
						 //CreatePathFrom(paths.get(i),paths.get(j), newPath);
						 //newPaths.push_back(&newPath);
					  }
			      }
			 }
		}
    	DEBUG("Paths with gaps: " << count);
    }



};

#endif /* PATH_GAPJUMPER_HPP_ */






