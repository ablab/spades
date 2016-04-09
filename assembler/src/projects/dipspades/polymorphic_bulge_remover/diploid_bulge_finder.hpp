//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "../utils/bulge_utils.hpp"
#include "bulge_paths_searcher.hpp"

using namespace debruijn_graph;

namespace dipspades {

class DiploidBulgeFinder {
    typedef vector<vector<EdgeId> > paths;

    Graph &graph_;
    double rel_length_;
    double rel_align_;

    bool RelativeLengthGood(double rel_length){
        return rel_length >= rel_length_;
    }

    bool RelativeAlignGood(double rel_align){
        return rel_align <= rel_align_;
    }

    bool BulgePathsIntersected(vector<EdgeId> &path1, vector<EdgeId> &path2){
        for(auto e1 = path1.begin(); e1 != path1.end(); e1++)
            for(auto e2 = path2.begin(); e2 != path2.end(); e2++)
                if(*e1 == *e2)
                    return true;
        return false;
    }

    vector<pair<size_t, size_t> > ChoosePairsWithoutIntersection(paths &bulge_paths){
        vector<pair<size_t, size_t> > correct_pairs;
        for(size_t i = 0; i < bulge_paths.size(); i++)
            for(size_t j = i + 1; j < bulge_paths.size(); j++)
                if(!BulgePathsIntersected(bulge_paths[i], bulge_paths[j]))
                    correct_pairs.push_back(make_pair(i, j));
        return correct_pairs;
    }

    vector<pair<size_t, size_t> > DefineLenSatisfiedPairs(vector<size_t> &lens,
            vector<pair<size_t, size_t> > pairs){
        vector<pair<size_t, size_t> > good_pairs;
        for(auto it = pairs.begin(); it != pairs.end(); it++)
            if(RelativeLengthGood(RelativeLengthEquality(lens[it->first], lens[it->second])))
                good_pairs.push_back(*it);
        return good_pairs;
    }

    vector<pair<size_t, size_t> > ChooseSeqSatisfiedPairs(vector<Sequence> &seqs,
            vector<pair<size_t, size_t> > pairs){
        vector<pair<size_t, size_t> > good_pairs;
        for(auto it = pairs.begin(); it != pairs.end(); it++)
            if(RelativeAlignGood(RelAlignmentOfSequences(seqs[it->first], seqs[it->second])))
                good_pairs.push_back(*it);
        return good_pairs;
    }

    vector<Sequence> GetSequences(paths &bulge_paths){
        vector<Sequence> seqs;
        for(auto it = bulge_paths.begin(); it != bulge_paths.end(); it++)
            seqs.push_back(GetSequenceByPath(graph_, graph_.k(), *it));
        return seqs;
    }

    vector<size_t> GetLengths(paths &bulge_paths){
            vector<size_t> lens;
            for(auto it = bulge_paths.begin(); it != bulge_paths.end(); it++)
                lens.push_back(GetPathLength(graph_, *it));
            return lens;
    }

public:
    DiploidBulgeFinder(Graph &graph, double rel_length, double rel_align) :
        graph_(graph),
        rel_length_(rel_length),
        rel_align_(rel_align) { }

    shared_ptr<BaseBulge> Find(paths &bulge_paths){
        if(bulge_paths.size() <= 1)
            return shared_ptr<BaseBulge>(new Bulge(graph_));

        auto good_pairs = ChoosePairsWithoutIntersection(bulge_paths);
        vector<Sequence> seqs = GetSequences(bulge_paths);
        vector<size_t> lens = GetLengths(bulge_paths);
        good_pairs = DefineLenSatisfiedPairs(lens, good_pairs);
        good_pairs = ChooseSeqSatisfiedPairs(seqs, good_pairs);

        if(good_pairs.size() == 0)
            return shared_ptr<BaseBulge>(new Bulge(graph_));
        return shared_ptr<BaseBulge>(new Bulge(graph_, graph_.k(), bulge_paths[good_pairs[0].first],
                bulge_paths[good_pairs[0].second]));
    }
};

}
