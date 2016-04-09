//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "path_routines.hpp"
#include "element_printers.hpp"

using namespace debruijn_graph;

namespace dipspades {

bool IsRegionBulge(Graph &g, vector<EdgeId> path1, vector<EdgeId> path2){
    if(path1.size() == 0 || path2.size() == 0)
        return true;
    if ((g.EdgeStart(path1[0]) != g.EdgeStart(path2[0])) ||
            (g.EdgeEnd(path1[path1.size() - 1]) != g.EdgeEnd(path2[path2.size() - 1])))
            return false;
    return !PathsShareEdge(path1, path2);
}

size_t AlignmentOfSequencesByParts(Sequence seq1, Sequence seq2){
    size_t max_length = 10000;
    if(min<size_t>(seq1.size(), seq2.size()) > max_length){
        size_t shrink1 = max_length;
        size_t num_full_iter = seq1.size() / shrink1;

        size_t summary_dist = 0;
        size_t shrink2 = size_t((double(shrink1) / double(seq1.size())) * double(seq2.size()));
        for(size_t i = 0; i < num_full_iter; i++){
            Sequence cur_seq1 = seq1.Subseq(shrink1 * i, shrink1 * (i + 1));
            Sequence cur_seq2 = seq2.Subseq(shrink2 * i, shrink2 * (i + 1));
            summary_dist += EditDistance(cur_seq1, cur_seq2);
        }

        if(shrink1 % seq1.size() != 0){
            Sequence cur_seq1 = seq1.Subseq(shrink1 * num_full_iter, seq1.size());
            Sequence cur_seq2 = seq2.Subseq(shrink2 * num_full_iter, seq2.size());
            summary_dist += EditDistance(cur_seq1, cur_seq2);
        }

        return summary_dist;
    }
    return EditDistance(seq1, seq2);
}

double RelAlignmentOfSequences(Sequence seq1, Sequence seq2){
    return double(AlignmentOfSequencesByParts(seq1, seq2)) / double(min<size_t>(seq1.size(), seq2.size()));
}

double RelativeLengthEquality(size_t len1, size_t len2){
    return double(min<size_t>(len1, len2)) / double(max<size_t>(len1, len2));
}

enum glue_direction { direct_gluing, reverse_gluing, undefined };

class BaseBulge{
protected:
    Graph &graph_;
public:
    BaseBulge(Graph &graph) : graph_(graph) { }
    BaseBulge(const BaseBulge& bulge) : graph_(bulge.graph_) {  }

    virtual double relative_length() = 0;
    virtual double relative_align() = 0;
    virtual bool IsBulgeDiploid(double rel_length_threshold, double rel_seq_threshold) = 0;
    virtual vector<EdgeId> path1() = 0; // todo make it const
    virtual vector<EdgeId> path2() = 0; // todo make it const
    virtual Sequence seq1() = 0;
    virtual Sequence seq2() = 0;
    virtual VertexId start_vertex() = 0;
    virtual VertexId end_vertex() = 0;
    virtual bool IsSimple() = 0;
    virtual bool IsEmpty() = 0;

    virtual size_t BulgeLength() = 0;

    virtual string StrId() = 0;
    virtual string BulgeToString() = 0;

    virtual ~BaseBulge() { }
};

class Bulge : public BaseBulge{
    size_t k_value_;
    vector<EdgeId> path1_;
    vector<EdgeId> path2_;
    Sequence seq1_;
    Sequence seq2_;
    double rel_length_;
    double rel_align_;

    void CalculateRelativeLength(size_t length1, size_t length2){
        rel_length_ = double(min<size_t>(length1, length2)) / double(max<size_t>(length1, length2));
    }

    void CalculateRelativeAlign(Sequence seq1, Sequence seq2){
        rel_align_ = RelAlignmentOfSequences(seq1, seq2);
    }

    string GetPathStr(vector<EdgeId> path) {
        string s1 = "";
        for(auto edge = path.begin(); edge != path.end(); edge++)
            s1 = ToString(graph_.int_id(*edge)) + "-";
        return s1.substr(0, s1.size() - 1);
    }

public:
    Bulge(Graph &graph) : BaseBulge(graph), k_value_(graph.k()), path1_(), path2_(),
    seq1_(), seq2_(), rel_length_(), rel_align_() { }

    Bulge(Graph &g, size_t k_value, vector<EdgeId> path1, pair<size_t,size_t> bulge_region1,
            vector<EdgeId> path2, pair<size_t,size_t> bulge_region2) :
    BaseBulge(g), k_value_(k_value),
    path1_(CutSubpathByRegion(path1, bulge_region1)),
    path2_(CutSubpathByRegion(path2, bulge_region2)),
    seq1_(GetSequenceByPath(graph_, k_value_, path1_).str().c_str()), // todo make it lazy
    seq2_(GetSequenceByPath(graph_, k_value_, path2_).str().c_str()), // todo make it lazy
    rel_length_(0), rel_align_(0) {
        VERIFY(IsRegionBulge(graph_, path1_, path2_));
    }

    Bulge(Graph &g, size_t k_value, vector<EdgeId> path1, vector<EdgeId> path2) :
    BaseBulge(g), k_value_(k_value), path1_(path1), path2_(path2),
    seq1_(GetSequenceByPath(graph_, k_value_, path1_).str().c_str()),
    seq2_(GetSequenceByPath(graph_, k_value_, path2_).str().c_str()),
    rel_length_(0), rel_align_(0){
        VERIFY(IsRegionBulge(graph_, path1_, path2_));
    }

    Bulge(Graph &g, size_t k_value, EdgeId edge1, EdgeId edge2) : BaseBulge(g),
            k_value_(k_value),
            path1_(1, edge1), path2_(1, edge2),
            seq1_(graph_.EdgeNucls(edge1).str().c_str()),
            seq2_(graph_.EdgeNucls(edge2).str().c_str()),
            rel_length_(0), rel_align_(0) {
        VERIFY(IsRegionBulge(graph_, path1_, path2_));
    }

    double relative_length(){
        if(rel_length_ == 0){
            size_t length1 = GetPathLength(graph_, path1_);
            size_t length2 = GetPathLength(graph_, path2_);
            CalculateRelativeLength(length1, length2);
        }
        return rel_length_;
    }

    double relative_align(){
        if(rel_align_ == 0){
            Sequence seq1 = GetSequenceByPath(graph_, k_value_, path1_);
            Sequence seq2 = GetSequenceByPath(graph_, k_value_, path2_);
            CalculateRelativeAlign(seq1, seq2);
        }
        return rel_align_;
    }

    bool IsBulgeDiploid(double rel_length_threshold, double rel_seq_threshold){
        if(relative_length() < rel_length_threshold)
            return false;

        return relative_align() <= rel_seq_threshold;
    }

    vector<EdgeId> path1(){
        return path1_;
    }

    vector<EdgeId> path2(){
        return path2_;
    }

    Sequence seq1() { return seq1_; }

    Sequence seq2() { return seq2_; }

    VertexId start_vertex(){
        return graph_.EdgeStart(path1_[0]);
    }

    VertexId end_vertex(){
            return graph_.EdgeEnd(path1_[path1_.size() - 1]);
    }

    bool IsSimple() { return path1_.size() == 1 && path2_.size() == 1; }

    bool IsEmpty() { return path1_.size() == 0 || path2_.size() == 0; }

    string StrId() {
        string s1 = GetPathStr(path1());
        string s2 = GetPathStr(path2());
        return min<string>(s1,s2) + "_" + max<string>(s1,s2);
    }

    size_t BulgeLength() {
        return max<size_t>(GetPathLength(graph_, path1()), GetPathLength(graph_, path2()));
    }

    string BulgeToString() {
        return "Side1: " + SimplePathWithVerticesToString(graph_, path1()) + "\n" +
                "Side2: " + SimplePathWithVerticesToString(graph_, path2());
    }
};

class DirectedBulge : public BaseBulge {
    shared_ptr<BaseBulge> bulge_;
    bool glue_direct_;
public:
    DirectedBulge(Graph &graph, shared_ptr<BaseBulge> bulge, glue_direction glue_direct = direct_gluing) :
        BaseBulge(graph), bulge_(bulge), glue_direct_(glue_direct) { }

    double relative_length() { return bulge_->relative_length(); }

    double relative_align() { return bulge_->relative_align(); }

    bool IsBulgeDiploid(double rel_length_threshold, double rel_seq_threshold) {
        return bulge_->IsBulgeDiploid(rel_length_threshold, rel_seq_threshold);
    }

    vector<EdgeId> path1() {
        if(glue_direct_ == direct_gluing)
            return bulge_->path1();
        return bulge_->path2();
    }

    vector<EdgeId> path2() {
        if(glue_direct_ == direct_gluing)
            return bulge_->path2();
        return bulge_->path1();
    }

    Sequence seq1() {
        if(glue_direct_ == direct_gluing)
            return bulge_->seq1();
        return bulge_->seq1();
    }

    Sequence seq2() {
        if(glue_direct_ == direct_gluing)
            return bulge_->seq2();
        return bulge_->seq2();
    }

    VertexId start_vertex() {
        return bulge_->start_vertex();
    }

    VertexId end_vertex() {
        return bulge_->end_vertex();
    }

    bool IsSimple() { return bulge_->IsSimple(); }

    bool IsEmpty() { return bulge_-> IsEmpty(); }

    string StrId() { return bulge_->StrId(); }

    size_t BulgeLength() { return bulge_->BulgeLength(); }

    string BulgeToString() { return bulge_->BulgeToString(); }
};

}
