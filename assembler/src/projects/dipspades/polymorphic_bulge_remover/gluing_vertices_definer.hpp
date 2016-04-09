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

size_t abs_diff(size_t a, size_t b){
    if(a > b)
        return a - b;
    return b - a;
}

typedef map<size_t, size_t> BulgeGluingVertices;
typedef map<size_t, size_t>::iterator BulgeGluingVerticesIter;

class GluingVericesDefinerResults {
    map<size_t, size_t> gluing_pairs_;

    bool IsNewPairGood(pair<size_t, size_t> new_pair){
        if(gluing_pairs_.size() == 0)
            return true;

        auto upper = gluing_pairs_.upper_bound(new_pair.first);
        auto lower = gluing_pairs_.lower_bound(new_pair.first);

        // map doesnt contain element with greater 1st vertex
        if(upper == gluing_pairs_.end()){
            // map doesnt contain element with the same 1st vertex
            if(lower == upper){
                // go to the previous element (with less 1st vertex)
                // if its 2nd vertex precedes 2nd vertex of new pair, add new pair
                // otherwise do not add new pair
                lower--;
                return lower->second < new_pair.second;
            }
            // map contains element with the same key.
            // addition of new pair is incorrect
            return false;
        }
        // map contains element with greater 1st vertex
        // and corresponding 2nd vertex is <= 2nd vertex of new pair
        // addition is incorrect
        if(upper->second <= new_pair.second)
            return false;
        // map contains element with greater 1st vertex
        // and doesnt contain the same 1st vertex
        if(lower == upper){
            // if there are no other elements
            // add new pair
            if(lower == gluing_pairs_.begin())
                return true;
            // there are other elements exist
            // go to the previous element
            // and check for preceding its 2nd vertex to new pair 2nd vertex
            lower--;
            return lower->second < new_pair.second;
        }
        return false;
    }

public:
    void AddNewPair(pair<size_t, size_t> new_pair){
        if(IsNewPairGood(new_pair)){
            TRACE("New pair was inserted: " << new_pair.first << " " << new_pair.second);
            gluing_pairs_.insert(new_pair);
        }
    }

    BulgeGluingVerticesIter begin() { return gluing_pairs_.begin(); }

    BulgeGluingVerticesIter end() { return gluing_pairs_.end(); }

    size_t size() { return gluing_pairs_.size(); }

private:
    DECL_LOGGER("GluingVericesDefinerResults");
};

class GluingVericesDefiner {
    Graph &graph_;
    double rel_length_threshold_;

    typedef map<pair<size_t, size_t>, double> PairSubpaths;
    PairSubpaths gluing_candidate_;

    double RelativeSimilarityOfLength(size_t len1, size_t len2){
        return double(abs_diff(len1, len2)) / double(min<size_t>(len1, len2));
    }

    size_t StartSpathLength(const vector<size_t> &lens, size_t index){
        VERIFY(index < lens.size());
        return lens[index];
    }

    size_t EndSpathLength(const vector<size_t> &lens, size_t index){
        VERIFY(index < lens.size());
        return lens[lens.size() - 1] - lens[index];
    }

    void ChooseGluingCandidate(shared_ptr<BaseBulge> bulge){

        TRACE("Choosing gluing candidates");
        vector<size_t> part_lens1 = CalculatePathPartLens(graph_, bulge->path1());
        vector<size_t> part_lens2 = CalculatePathPartLens(graph_, bulge->path2());

        for(size_t i = 0; i < part_lens1.size() - 1;  i++)
            for(size_t j = 0; j < part_lens2.size() - 1; j++){
                double rel_len_start_spaths = RelativeSimilarityOfLength(StartSpathLength(part_lens1, i),
                        StartSpathLength(part_lens2, j));
                double rel_len_end_spaths = RelativeSimilarityOfLength(EndSpathLength(part_lens1, i),
                        EndSpathLength(part_lens2, j));

                if(rel_len_start_spaths <= rel_length_threshold_ &&
                        rel_len_end_spaths <= rel_length_threshold_){
                    TRACE("New gluing candidate - " << i << ", " << j);
                    TRACE("rel_len_start_spaths - " << rel_len_start_spaths);
                    TRACE("rel_len_end_spaths - " << rel_len_end_spaths);
                    gluing_candidate_[make_pair(i,j)] = max<double>(rel_len_start_spaths, rel_len_end_spaths);
                }
            }
    }

    pair<size_t, size_t> GetBestPair(){
        double min = 1;
        pair<size_t, size_t> best_res;

        for(auto it = gluing_candidate_.begin(); it != gluing_candidate_.end(); it++)
            if(it->second < min){
                best_res = it->first;
                min = it->second;
            }
        return best_res;
    }

    GluingVericesDefinerResults ChooseGluingPairs(shared_ptr<BaseBulge> bulge){
        gluing_candidate_.clear();
        ChooseGluingCandidate(bulge);
        GluingVericesDefinerResults gluing_pairs;
                while(gluing_candidate_.size() != 0){
            auto best_pair = GetBestPair();
            gluing_pairs.AddNewPair(best_pair);
            gluing_candidate_.erase(best_pair);
        }
        return gluing_pairs;
    }

public:
    GluingVericesDefiner(Graph &graph, double rel_length_threshold) :
        graph_(graph),
        rel_length_threshold_(rel_length_threshold) {    }

    GluingVericesDefinerResults Run(shared_ptr<BaseBulge> bulge){
        return ChooseGluingPairs(bulge);
    }

private:
    DECL_LOGGER("GluingVericesDefiner");
};

}
