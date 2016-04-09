//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "abstract_contig_corrector.hpp"

using namespace debruijn_graph;

namespace dipspades {

class LoopBulgeDeletionCorrector : public AbstractContigCorrector{

protected:
    size_t k_value_;

    size_t max_loop_length_;
    size_t max_tail_length_;
    size_t min_lcs_length_;
    VertexPathIndex &path_index_;

    ostream &out_;

public:
    vector<VertexId> GetListOfVertices(vector<EdgeId> path){
        return get_list_of_vertices_in_path(g_, path);
    }

    size_t ConvInd(size_t ind, size_t size){
        if(ind == size)
            return ind - 1;
        return ind;
    }

    bool IsLoopCorrect(vector<EdgeId> path, size_t i1, size_t i2){
        VERIFY(i1 < path.size());
        VERIFY(i2 < path.size());

        size_t length = 0;
        for(size_t i = i1; i <= i2; i++)
            length += g_.length(path[i]);

        return length <= max_loop_length_;
    }

    bool IsRegionLoop(vector<EdgeId> path, size_t i1, size_t i2){
        return g_.EdgeStart(path[i1]) == g_.EdgeEnd(path[i2]);
    }

    bool IsPathCorrect(vector<EdgeId> path, vector<size_t> pos){
        return (IsLeftTailCorrect(path, pos) && IsRightTailCorrect(path, pos));
    }

    bool IsRegionBulgeSide(vector<EdgeId> path, size_t ind1, size_t ind2){
        return g_.EdgeStart(path[ind1]) != g_.EdgeEnd(path[ind2]);
    }

    bool AreRegionsBulge(vector<EdgeId> path1, size_t i_11, size_t i_12,
            vector<EdgeId> path2, size_t i_21, size_t i_22){
        return IsRegionBulge(g_, CutSubpathByRegion(path1, make_pair(i_11, i_12)),
                CutSubpathByRegion(path2, make_pair(i_21, i_22)));
    }

    bool AreRegionsDiploidBulge(vector<EdgeId> path1, size_t i_11, size_t i_12,
            vector<EdgeId> path2, size_t i_21, size_t i_22){

        TRACE("Bulge: pos1: " << i_11 << " - " << i_12 << ", pos2: " << i_21 << " - " << i_22 );

        if(dsp_cfg::get().cc.align_bulge_sides){
            Bulge bulge(g_, k_value_, path1, make_pair(i_11, i_12), path2, make_pair(i_21, i_22));
            return bulge.IsBulgeDiploid(dsp_cfg::get().pbr.rel_bulge_align,
                    dsp_cfg::get().pbr.rel_bulge_align);
        }
        return true;
    }

    bool IsLCSCorrect(vector<EdgeId> path1, vector<size_t> pos1,
            vector<EdgeId> path2, vector<size_t> pos2){

        VERIFY(pos1.size() == pos2.size());

        size_t pos_len = pos1.size();
        if(pos_len <= 1)
            return false;

        size_t lcs_len = min<size_t>(GetLCSLengthByPath(path1, pos1), GetLCSLengthByPath(path2, pos2));
        size_t path1_len = GetPathLength(g_, path1), path2_len = GetPathLength(g_, path2);

        TRACE("LCS length - " << lcs_len);
        TRACE("Path length1 - " << path1_len << ", path length2 - " << path2_len);

        if(lcs_len <= min_lcs_length_ &&
                min<size_t>(path1_len, path2_len) > min_lcs_length_){
            return false;
        }

        for(size_t i = 0; i < pos_len - 1; i++){

            TRACE("Pos1 - " << pos1[i] << ", " << pos1[i + 1]);
            TRACE("Pos2 - " << pos2[i] << ", " << pos2[i + 1]);
            // if bath are not neighbors
            if(pos1[i] + 1 != pos1[i + 1] || pos2[i] + 1 != pos2[i + 1]){

                TRACE("1st loop checking");
                bool is_1st_loop = false;
                bool is_1st_corr = true;

                size_t i_11, i_12;
                if(pos1[i] + 1 != pos1[i + 1]){
                    TRACE("Positions are not consecutive");
                    // it may be loop
                    i_11 = ConvInd(pos1[i], path1.size());
                    i_12 = ConvInd(pos1[i + 1], path1.size()) - 1;

                    is_1st_loop = IsRegionLoop(path1, i_11, i_12);
                    TRACE("Is loop - " << is_1st_loop);
                    if(is_1st_loop){
                        is_1st_corr = IsLoopCorrect(path1, i_11, i_12);
                    }
                    else{ // then region is bulge
                        VERIFY(IsRegionBulgeSide(path1, i_11, i_12));
                    }
                }
                else{
                    i_11 = pos1[i];
                    i_12 = pos1[i];
                }

                TRACE("2nd loop checking");
                bool is_2nd_loop = false;
                bool is_2nd_corr = true;
                size_t i_21, i_22;
                if(pos2[i] + 1 != pos2[i + 1]){
                    TRACE("Positions are not consecutive");
                    // it may be loop
                    i_21 = ConvInd(pos2[i], path2.size());
                    i_22 = ConvInd(pos2[i + 1], path2.size()) - 1;

                    is_2nd_loop = IsRegionLoop(path2, i_21, i_22);
                    TRACE("Is loop - " << is_2nd_loop );
                    if(is_2nd_loop){
                        is_2nd_corr = IsLoopCorrect(path2, i_21, i_22);
                    }
                    else{
                        VERIFY(IsRegionBulgeSide(path2, i_21, i_22));
                    }
                }
                else{
                    i_21 = pos2[i];
                    i_22 = pos2[i];
                }

                if(!is_1st_loop && !is_2nd_loop){

                    i_12 = (pos1[i + 1] == path1.size()) ? path1.size() - 1 : i_12;
                    i_22 = (pos2[i + 1] == path2.size()) ? path2.size() - 1 : i_22;

                    if(AreRegionsBulge(path1, i_11, i_12, path2, i_21, i_22))
                        if(!AreRegionsDiploidBulge(path1, i_11, i_12, path2, i_21, i_22))
                            return false;
                }
                else
                    if(!is_1st_corr || !is_2nd_corr){
                        return false;
                    }
            }
        }

        return true;
    }

    size_t GetLCSLengthByPath(vector<EdgeId> path, vector<size_t> pos){
        if(pos.size() <= 1)
            return 0;
        size_t pos_len = pos.size();
        size_t last_pos = pos[pos_len - 1];

        size_t ind_start = ConvInd(pos[0], path.size());
        size_t ind_end = ConvInd(last_pos, path.size());
        if(last_pos != path.size())
            ind_end--;

        size_t len = 0;
        for(size_t i = ind_start; i <= ind_end; i++)
            len += g_.length(path[i]);

        return len;
    }

    size_t GetLeftTailLength(vector<EdgeId> path, vector<size_t> pos){
        if(pos.size() <= 1)
            return 0;

        size_t first_pos = ConvInd(pos[0], path.size());

        size_t tail_len = 0;
        for(size_t i = 0; i < first_pos; i++)
            tail_len += g_.length(path[i]);

        return tail_len;
    }

    bool IsLeftTailCorrect(vector<EdgeId> path, vector<size_t> pos){
        if(pos.size() <= 1)
            return false;

        size_t tail_len = GetLeftTailLength(path, pos);
//        TRACE("Left tail length - " << tail_len );
        return (tail_len <= max_tail_length_);
    }

    size_t GetRightTailLength(vector<EdgeId> path, vector<size_t> pos){
        if(pos.size() <= 1)
            return 0;

        size_t last_pos = pos[pos.size() - 1];

        size_t tail_len = 0;
        for(size_t i = last_pos; i < path.size(); i++)
            tail_len += g_.length(path[i]);

        return tail_len;
    }

    bool IsRightTailCorrect(vector<EdgeId> path, vector<size_t> pos){
        if(pos.size() <= 1)
            return false;

        size_t tail_len = GetRightTailLength(path, pos);
//        TRACE("Right tail length - " << tail_len );
        return (tail_len <= max_tail_length_);
    }

    bool AreLeftTailsCorrect(vector<EdgeId> path1, vector<size_t> pos1,
            vector<EdgeId> path2, vector<size_t> pos2){

        VERIFY(pos1.size() == pos2.size());

        if(pos1.size() <= 1)
            return false;

        size_t tail_length1 = GetLeftTailLength(path1, pos1);
        size_t tail_length2 = GetLeftTailLength(path2, pos2);

        if(min<size_t>(tail_length1, tail_length2) > max_tail_length_)
            return false;

        VertexId start1 = g_.EdgeStart(path1[0]); //g_.EdgeStart(path1[first_pos1]);
        VertexId start2 = g_.EdgeStart(path2[0]); //g_.EdgeStart(path2[first_pos2]);

        bool are_tails_correct = false;
        if(g_.IncomingEdgeCount(start1) == 0 &&
                g_.IncomingEdgeCount(start2) == 0){
            are_tails_correct = true;
        }
        else{

            if(dsp_cfg::get().cc.tails_lie_on_bulges){
                // find vertex v, such that paths starts are reachable from v
                auto path_searcher1 = DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(g_,
                        max_tail_length_);
                path_searcher1.Run(start1);
                auto reached_vert1 = path_searcher1.ReachedVertices();

                auto path_searcher2 = DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(g_,
                        max_tail_length_);
                path_searcher2.Run(start2);
                auto reached_vert2 = path_searcher2.ReachedVertices();

                for(size_t i = 0; i < pos1[0]; i++)
                    reached_vert1.push_back(g_.EdgeStart(path1[i]));

                for(size_t i = 0; i < pos2[0]; i++)
                    reached_vert2.push_back(g_.EdgeStart(path2[i]));

                for(auto v1 = reached_vert1.begin(); v1 != reached_vert1.end(); v1++){
                    for(auto v2 = reached_vert2.begin(); v2 != reached_vert2.end(); v2++){
                        if(*v1 == *v2){
                            are_tails_correct = true;
                            break;
                        }
                    }
                    if(are_tails_correct)
                        break;
                }
            }
            else{
                are_tails_correct = true;
            }
        }

        if(!are_tails_correct)
            return false;

        if(!dsp_cfg::get().cc.align_bulge_sides)
            return true;

        Sequence tail_seq1 = GetSequenceOfPathRegion(g_, k_value_, path1,
                pair<size_t, size_t>(0, pos1[0] - 1));

        Sequence tail_seq2 = GetSequenceOfPathRegion(g_, k_value_, path2,
                pair<size_t, size_t>(0, pos2[0] - 1));

        Sequence trim_seq1, trim_seq2;
        if(min<size_t>(tail_seq1.size(), tail_seq2.size()) == tail_seq1.size()){
            trim_seq1 = tail_seq1;
            trim_seq2 = tail_seq2.Subseq(tail_seq2.size() - tail_seq1.size(),
                    tail_seq2.size());
        }
        else{
            trim_seq1 = tail_seq1.Subseq(tail_seq1.size() - tail_seq2.size(),
                    tail_seq1.size());
            trim_seq2 = tail_seq2;
        }

        if(trim_seq1.size() > max_tail_length_)
            return false;

        return RelAlignmentOfSequences(trim_seq1, trim_seq2) <=
                dsp_cfg::get().pbr.rel_bulge_align;

    }

    bool AreRightTailsCorrect(vector<EdgeId> path1, vector<size_t> pos1,
                vector<EdgeId> path2, vector<size_t> pos2){

            VERIFY(pos1.size() == pos2.size());

            if(pos1.size() <= 1)
                return false;

            size_t tail_length1 = GetRightTailLength(path1, pos1);
            size_t tail_length2 = GetRightTailLength(path2, pos2);

            if(min<size_t>(tail_length1, tail_length2) > max_tail_length_)
                return false;

            VertexId end1 = g_.EdgeEnd(path1[path1.size() - 1]);
            VertexId end2 = g_.EdgeEnd(path2[path2.size() - 1]);

            bool are_tails_correct = false;

            if(g_.OutgoingEdgeCount(end1) == 0 && g_.OutgoingEdgeCount(end2) == 0){
                are_tails_correct = true;
            }
            else{

                if(dsp_cfg::get().cc.tails_lie_on_bulges){
                    // find vertex v, such that paths ends are reachable from v
                    auto path_searcher1 = DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(g_,
                            max_tail_length_);
                    path_searcher1.Run(end1);
                    auto reached_vert1 = path_searcher1.ReachedVertices();

                    auto path_searcher2 = DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(g_,
                            max_tail_length_);
                    path_searcher2.Run(end2);
                    auto reached_vert2 = path_searcher2.ReachedVertices();

                    for(size_t i = ConvInd(pos1[pos1.size() - 1], path1.size()); i < path1.size(); i++)
                        reached_vert1.push_back(g_.EdgeEnd(path1[i]));

                    for(size_t i = ConvInd(pos2[pos2.size() - 1], path2.size()); i < path2.size(); i++)
                        reached_vert2.push_back(g_.EdgeEnd(path2[i]));

                    for(auto v1 = reached_vert1.begin(); v1 != reached_vert1.end(); v1++){
                        for(auto v2 = reached_vert2.begin(); v2 != reached_vert2.end(); v2++){
                            if(*v1 == *v2){
                                are_tails_correct = true;
                                break;
                            }
                        }
                        if(are_tails_correct)
                            break;
                    }
                }
                else{
                    // tail lengths comparison?
                    are_tails_correct = true;
                }
            }

            if(!are_tails_correct)
                return false;

            if(!dsp_cfg::get().cc.align_bulge_sides)
                return true;

            Sequence tail_seq1 = GetSequenceOfPathRegion(g_, k_value_, path1,
                    pair<size_t,size_t>(pos1[pos1.size() - 1], path1.size() - 1));

            Sequence tail_seq2 = GetSequenceOfPathRegion(g_, k_value_, path2,
                    pair<size_t,size_t>(pos2[pos2.size() - 1], path2.size() - 1));

            Sequence trim_seq1, trim_seq2;
            if(min<size_t>(tail_seq1.size(), tail_seq2.size()) == tail_seq1.size()){
                trim_seq1 = tail_seq1;
                trim_seq2 = tail_seq2.Subseq(0, tail_seq1.size());
            }
            else{
                trim_seq1 = tail_seq1.Subseq(0, tail_seq2.size());
                trim_seq2 = tail_seq2;
            }

            if(trim_seq1.size() > max_tail_length_)
                return false;

            return (RelAlignmentOfSequences(trim_seq1, trim_seq2) <= dsp_cfg::get().pbr.rel_bulge_align);
        }

    bool IsLeftTailExist(vector<EdgeId>, vector<size_t> pos){
        size_t first_index = pos[0]; //ConvInd(pos[0], path.size());
        return (first_index != 0);
    }

    bool AreBothLeftTailsExist(vector<EdgeId> path1, vector<size_t> pos1,
            vector<EdgeId> path2, vector<size_t> pos2){

        VERIFY(pos1.size() == pos2.size());
        if(pos1.size() == 0)
            return false;

        TRACE("Left: " << IsLeftTailExist(path1, pos2) << " " << IsLeftTailExist(path2, pos2) );
        return (IsLeftTailExist(path1, pos1) && IsLeftTailExist(path2, pos2));
    }

    bool IsRightTailExist(vector<EdgeId> path, vector<size_t> pos){
        size_t last_index = pos[pos.size() - 1]; //ConvInd(pos[pos.size() - 1], path.size());
        return (last_index != path.size());
    }

    bool AreBothRightTailsExist(vector<EdgeId> path1, vector<size_t> pos1,
            vector<EdgeId> path2, vector<size_t> pos2){
        VERIFY(pos1.size() == pos2.size());
        if(pos1.size() == 0)
            return false;

        TRACE("Right: " << IsRightTailExist(path1, pos1) << " " << IsRightTailExist(path2, pos2) );
        return IsRightTailExist(path1, pos1) && IsRightTailExist(path2, pos2);

    }

    // yana todo replace
    vector<VertexId> RearrangementSearch(vector<EdgeId> path1, vector<EdgeId> path2){

        vector<VertexId> common_vertices;
        if(path1.size() == 0 || path2.size() == 0)
            return common_vertices;

        map<VertexId, int> vertex_count;

        set<VertexId> vertices1; vertices1.insert(g_.EdgeStart(path1[0]));
        for(auto e = path1.begin(); e != path1.end(); e++){
            vertices1.insert(g_.EdgeEnd(*e));
        }

        set<VertexId> vertices2; vertices2.insert(g_.EdgeStart(path2[0]));
        for(auto e = path2.begin(); e != path2.end(); e++){
            vertices2.insert(g_.EdgeEnd(*e));
        }

        for(auto v = vertices1.begin(); v != vertices1.end(); v++)
            vertex_count[*v]++;

        for(auto v = vertices2.begin(); v != vertices2.end(); v++)
            vertex_count[*v]++;

        for(auto it = vertex_count.begin(); it != vertex_count.end(); it++)
            if(it->second == 2)
                common_vertices.push_back(it->first);

//        TRACE("Common vertices: " );
//        PrintVectorOfVertices(cout, g_, common_vertices);
        return common_vertices;
    }

    bool ArePathsCorrect(vector<EdgeId> path1, vector<size_t> pos1,
            vector<EdgeId> path2, vector<size_t> pos2){

        VERIFY(pos1.size() == pos2.size());

        if(AreBothLeftTailsExist(path1, pos1, path2, pos2))
        {
            TRACE("Both left tails exist" );
            bool tail_corr = AreLeftTailsCorrect(path1, pos1, path2, pos2);
            if(!tail_corr){
                TRACE("One of left tails is not correct" );
                return false;
            }
        }

        if(AreBothRightTailsExist(path1, pos1, path2, pos2))
        {
            TRACE("Both right tails exist" );
            bool tail_corr = AreRightTailsCorrect(path1, pos1, path2, pos2);
            if(!tail_corr){
                TRACE("One of right tails is not correct" );
                return false;
            }
        }

        bool lcs_corr = IsLCSCorrect(path1, pos1, path2, pos2);

        if(!lcs_corr){
            TRACE("LCS is not correct" );
            auto common_vert = RearrangementSearch(path1, path2);
            if(common_vert.size() > pos1.size())
                TRACE("Possible rearrangement!");
            return false;
        }

        size_t lcs_length1 = GetLCSLengthByPath(path1, pos1),
                lcs_length2 = GetLCSLengthByPath(path2, pos2);

        return (min<size_t>(lcs_length1, lcs_length2) > 0);
    }

    bool IsPathRedundant(vector<EdgeId> path, vector<size_t> pos){
        if(pos.size() <= 1) return true;
        return IsLeftTailCorrect(path, pos) && IsRightTailCorrect(path, pos);
    }

    void CorrectPositionVertor(vector<EdgeId> path, vector<size_t> & pos){
        if(pos.size() <= 1)
            return;

        for(size_t i = 0; i < pos.size() - 1; i++){
            if(pos[i] + 1 != pos[i + 1]){

                size_t i1 = ConvInd(pos[i], path.size()) + 1;
                size_t i2 = ConvInd(pos[i + 1], path.size()) - 1;

                if(IsRegionLoop(path, i1, i2))
                    if(!IsLoopCorrect(path, i1, i2)){
                        VertexId v;
                        if(pos[i + 1] == path.size())
                            v = g_.EdgeEnd(path[path.size() - 1]);
                        else
                            v = g_.EdgeStart(path[pos[i + 1]]);
                        for(size_t j = pos[i] + 1; j < pos[i + 1]; j++)
                            if(g_.EdgeStart(path[j]) == v){
                                pos[i + 1] = j;
                                break;
                            }
                    }
            }
        }
    }

    size_t GetNumberOfErrorsFromLCS(vector<EdgeId> path, vector<size_t> pos){
        if(pos.size() <= 1)
            return 0;

        size_t error_num = 0;

        for(size_t i = 0; i < pos.size() - 1; i++){
            if(pos[i] + 1 != pos[i + 1]){

                size_t i1 = ConvInd(pos[i], path.size()) + 1;
                size_t i2 = ConvInd(pos[i + 1], path.size()) - 1;

                if(IsRegionLoop(path, i1, i2)){
                    if(!IsLoopCorrect(path, i1, i2))
                        error_num++;
                }
            }
        }

        return error_num;
    }

    pair<vector<size_t>, vector<size_t> > GetBestPosVectors(LCSCalculator<VertexId> & calc,
            vector<EdgeId> path1, vector<VertexId> vert_path1,
            vector<EdgeId> path2, vector<VertexId> vert_path2,
            vector<VertexId> lcs){

        // first path processing
        auto pos_right1 = calc.GetPosVector(vert_path1, lcs);
        auto pos_left1 = calc.GetPosVectorFromLeft(vert_path1, lcs);

        bool equal_num_err1 = true;
        vector<size_t> best_vect1;

        {
            size_t err_right1 = GetNumberOfErrorsFromLCS(path1, pos_right1);
            size_t err_left1 = GetNumberOfErrorsFromLCS(path1, pos_left1);
            equal_num_err1 = err_left1 == err_right1;
            best_vect1 = (err_left1 < err_right1) ? pos_left1 : pos_right1;
        }

        size_t lcs_right1 = GetLCSLengthByPath(path1, pos_right1);
        size_t lcs_left1 = GetLCSLengthByPath(path1, pos_left1);

        // second path processing
        auto pos_right2 = calc.GetPosVector(vert_path2, lcs);
        auto pos_left2 = calc.GetPosVectorFromLeft(vert_path2, lcs);

        bool equal_num_err2 = true;
        vector<size_t> best_vect2;

        {
            size_t err_right2 = GetNumberOfErrorsFromLCS(path2, pos_right2);
            size_t err_left2 = GetNumberOfErrorsFromLCS(path2, pos_left2);
            equal_num_err2 = err_left2 == err_right2;
            best_vect2 = (err_left2 < err_right2) ? pos_left2 : pos_right2;
        }

        size_t lcs_right2 = GetLCSLengthByPath(path2, pos_right2);
        size_t lcs_left2 = GetLCSLengthByPath(path2, pos_left2);

        if(equal_num_err1 && !equal_num_err2){

            size_t best_lcs2 = GetLCSLengthByPath(path2, best_vect2);

            if(abs_diff(lcs_right1, best_lcs2) < abs_diff(lcs_left1, best_lcs2))
                return pair<vector<size_t>, vector<size_t> >(pos_right1, best_vect2);
            else
                return pair<vector<size_t>, vector<size_t> >(pos_left1, best_vect2);
        }

        if(!equal_num_err1 && equal_num_err2){
            size_t best_lcs1 = GetLCSLengthByPath(path1, best_vect1);

            if(abs_diff(lcs_right2, best_lcs1) < abs_diff(lcs_left2, best_lcs1))
                return pair<vector<size_t>, vector<size_t> >(best_vect1, pos_right2);
            else
                return pair<vector<size_t>, vector<size_t> >(best_vect1, pos_left2);
        }

        if(equal_num_err1 && equal_num_err2){

            // best pair computing
            size_t left_left = abs_diff(lcs_left1, lcs_left2);
            size_t left_right = abs_diff(lcs_left1, lcs_right2);
            size_t right_left = abs_diff(lcs_right1, lcs_left2);
            size_t right_right = abs_diff(lcs_right1, lcs_right2);

            size_t min_diff = min<size_t>(min<size_t>(left_left, left_right),
                    min<size_t>(right_left, right_right));

            if(min_diff == left_left){
                return pair<vector<size_t>, vector<size_t> >(pos_left1, pos_left2);
            }

            if(min_diff == left_right){
                return pair<vector<size_t>, vector<size_t> >(pos_left1, pos_right2);
            }

            if(min_diff == right_left){
                return pair<vector<size_t>, vector<size_t> >(pos_right1, pos_left2);
            }

            if(min_diff == right_right){
                return pair<vector<size_t>, vector<size_t> >(pos_right1, pos_right2);
            }
        }

        return pair<vector<size_t>, vector<size_t> >(best_vect1, best_vect2);
    }

    vector<size_t> GetBestPosVector(LCSCalculator<VertexId> & calc, vector<EdgeId> path,
            vector<VertexId> vert_path, vector<VertexId> lcs){

        auto pos_right = calc.GetPosVector(vert_path, lcs);
        auto pos_left = calc.GetPosVectorFromLeft(vert_path, lcs);

        size_t err_right = GetNumberOfErrorsFromLCS(path, pos_right);
        size_t err_left = GetNumberOfErrorsFromLCS(path, pos_left);

        if(min<size_t>(err_left, err_right) == err_left)
            return pos_left;
        else
            return pos_right;
    }

    void InitializeMap(ContigStoragePtr contigs){
        for(size_t i = 0; i < contigs->Size(); i++){
            size_t id = (*contigs)[i]->id();
            res.redundancy_map.AddNewKey(id);
        }
    }

    void AddRedundantContig(ContigStoragePtr contigs, size_t index_red, size_t index_main){
        size_t id_main = (*contigs)[index_main]->id(),
                id_rc_main = (*contigs)[index_main]->rc_id();
        size_t id_red = (*contigs)[index_red]->id(),
                id_rc_red = (*contigs)[index_red]->rc_id();
        redundant_contigs.insert(id_red);
        redundant_contigs.insert(id_rc_red);
        // current contig
        res.redundancy_map.AddNewPair(id_main, id_red);
        res.redundancy_map.AddNewPair(id_rc_main, id_rc_red);
    }

    CorrectionResult res;
    set<size_t> redundant_contigs;

public:
    LoopBulgeDeletionCorrector(Graph &g, size_t k_value, size_t max_loop_length,
            size_t max_tail_length, size_t min_lcs_length, VertexPathIndex &path_index,
            ostream &out = cout) : AbstractContigCorrector(g), k_value_(k_value),
            path_index_(path_index), out_(out) {

        max_loop_length_ = max_loop_length;
        max_tail_length_ = max_tail_length;
        min_lcs_length_ = min_lcs_length;
    }

    virtual ContigStoragePtr Correct(ContigStoragePtr contigs)    {

        INFO("Computing redundant contigs starts");

        redundant_contigs.clear();

        InitializeMap(contigs);

        LCSCalculator<VertexId> lcs_calc;

        vector<vector<VertexId> > seqs;
        for(size_t i = 0; i < contigs->Size(); i++){
            vector<VertexId> seq = GetListOfVertices((*contigs)[i]->path_seq());
            seqs.push_back(seq);
        }

        set<size_t> processed_contigs;
        set<size_t> absolutely_redundant;

        size_t contigs_number = seqs.size();
        double processed_perc = 0.1;
        double processed_step = 0.1;

        for(size_t i = 0; i < seqs.size() - 1; i++){

            size_t id_i = (*contigs)[i]->id();
            size_t rc_id_i = (*contigs)[i]->rc_id();

            processed_contigs.insert(id_i);

            if(processed_contigs.find(rc_id_i) == processed_contigs.end() &&
                    absolutely_redundant.find(i) == absolutely_redundant.end()){

                vector<EdgeId> path1 = (*contigs)[i]->path_seq();
                set<int> analyzed_contigs;

                auto contigs_for_analyze = path_index_.GetPathsIntersectedWith(path1);
                for(auto it = contigs_for_analyze.begin(); it != contigs_for_analyze.end(); it++){

                    size_t j = *it;
                    size_t id_j = (*contigs)[j]->id();
                    size_t rc_id_j = (*contigs)[j]->rc_id();

                    bool need_process = !((i % 2 == 0 && i + 1 == j) || j <= i);
                    need_process = need_process &&
                            absolutely_redundant.find(j) == absolutely_redundant.end();
                    if(need_process){

                        vector<EdgeId> path2 = (*contigs)[j]->path_seq();
                        vector<VertexId> lcs_res = lcs_calc.LCS(seqs[i], seqs[j]);
                        vector<size_t> pos1, pos2;

                        auto pos_vectors_pair = GetBestPosVectors(lcs_calc, path1, seqs[i], path2, seqs[j], lcs_res);
                        pos1 = pos_vectors_pair.first;
                        pos2 = pos_vectors_pair.second;

                        {
                            TRACE("--------------------------------");
                            TRACE("Indexes " << i << " " << j);
                            TRACE("IDs " << id_i << " " << id_j);
                            TRACE("RC_Ids " << rc_id_i << " " << rc_id_j);

                            TRACE("Path1. " << SimplePathWithVerticesToString(g_, path1));
                            TRACE("Path2. " << SimplePathWithVerticesToString(g_, path2));

                            TRACE("LCS string: " << VerticesVectorToString(g_, lcs_res));

                            TRACE("Pos1. " << VectorToString<size_t>(pos1));
                            TRACE("Pos2. " << VectorToString<size_t>(pos2));
                        }

                        if(pos1.size() > 1){

                            bool paths_corr = ArePathsCorrect(path1, pos1, path2, pos2);

                            {
                                TRACE("ArePathsCorrect - " << paths_corr);
                            }

                            if(paths_corr){

                                size_t first_tail1 = GetLeftTailLength(path1, pos1);
                                size_t first_tail2 = GetRightTailLength(path1, pos1);
                                size_t first_tails = first_tail1 + first_tail2;

                                size_t second_tail1 = GetLeftTailLength(path2, pos2);
                                size_t second_tail2 = GetRightTailLength(path2, pos2);
                                size_t second_tails = second_tail1 + second_tail2;

                                bool first_path_red = IsPathRedundant(path1, pos1);
                                bool second_path_red = IsPathRedundant(path2, pos2);

                                {
                                    TRACE("\tFirst tails length - " << first_tails);
                                    TRACE("\tFirst path is redundant - " << first_path_red);
                                    TRACE("\tSecond tails length - " << second_tails);
                                    TRACE("\tSecond path is redundant - " << second_path_red);
                                }

                                if(first_path_red && second_path_red){
                                    if(first_tails < second_tails){
                                        TRACE(id_i << " is redundant");
                                        AddRedundantContig(contigs, i, j);

                                        if(first_tails == 0)
                                            absolutely_redundant.insert(i);
                                    }
                                    else{
                                        TRACE(id_j << " is redundant");
                                        AddRedundantContig(contigs, j, i);

                                        if(second_tails == 0)
                                            absolutely_redundant.insert(j);
                                    }
                                }
                                else{
                                    if(first_path_red && !second_path_red){
                                        TRACE(id_i << " is redundant");
                                        AddRedundantContig(contigs, i, j);

                                        if(first_tails == 0)
                                            absolutely_redundant.insert(i);

                                    }
                                    else
                                        if(!first_path_red && second_path_red){
                                            TRACE(id_j << " is redundant");
                                            AddRedundantContig(contigs, j, i);

                                            if(second_tails == 0)
                                                absolutely_redundant.insert(j);
                                        }
                                }

                                if(absolutely_redundant.find(i) != absolutely_redundant.end())
                                    break;
                            }
                    }
                    }
                }
            }

            double cur_process_perc = static_cast<double>(i) / static_cast<double>(contigs_number);
            if(cur_process_perc > processed_perc) {
                while(processed_perc + processed_step<= cur_process_perc)
                    processed_perc += processed_step;
                INFO(ToString(processed_perc * 100.0) << "% contigs were processed");
                processed_perc += processed_step;
            }
        }
        INFO("100% contigs were processed");

        RedundancyMapCondenser<size_t> condenser;
        condenser.Condense(res.redundancy_map);

        INFO(ToString(redundant_contigs.size()) + " contigs from " + ToString(contigs->Size()) + " are redundant");

        contigs->DeleteByIDs(redundant_contigs);

        INFO("Computing redundant contigs ends");

        return contigs;
    }

    MappingContigPtr Correct(MappingContigPtr contig){
        return contig;
    }

    CorrectionResult Results(){
        return res;
    }

    virtual ~LoopBulgeDeletionCorrector(){}

protected:
    DECL_LOGGER("LoopBulgeDeletionCorrector");
};

}
