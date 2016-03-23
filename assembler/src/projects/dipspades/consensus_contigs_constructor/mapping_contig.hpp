//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * mapping_contig.hpp
 *
 *  Created on: 13.11.2012
 *      Author: yana
 */

#pragma once
#include "overlap_graph.hpp"
#include "../utils/element_printers.hpp"
#include <map>

using namespace debruijn_graph;

namespace dipspades {

class Loop{
    vector<EdgeId> list_;

public:
    Loop(vector<EdgeId> list) : list_(list) {}
    const vector<EdgeId> edges() { return list_; }
    size_t Size() { return list_.size(); }
};

class MappingContig{
public:
    MappingContig() { }
    virtual Sequence seq() = 0;
    virtual vector<EdgeId> path_seq() = 0;
    virtual MappingPath<EdgeId> mapping_path() = 0;

    virtual string name() = 0;
    virtual string src_file() = 0;
    virtual string full_name() { return src_file() + ":" + name(); }

    virtual size_t size() = 0;
    virtual size_t length() = 0;
    virtual size_t id(){ return -1; }
    virtual size_t rc_id() = 0;
    virtual vector<shared_ptr<MappingContig> > AllMappingContigs() = 0;
    virtual void ChangeMappingRange(size_t, MappingRange){    }
    virtual void ChangeName(string new_name) = 0;

    virtual string ToString(Graph &graph) = 0;

    virtual ~MappingContig(){}
};

typedef shared_ptr<MappingContig> MappingContigPtr;

class SimpleMappingContig : public MappingContig{
    string name_;
    string src_file_;
    Sequence seq_;
    MappingPath<EdgeId> map_path_;
    vector<EdgeId> edge_path_;
    size_t id_, rc_id_;

public:
    SimpleMappingContig(){ }

    SimpleMappingContig(string name, string src_file, Sequence seq,
            MappingPath<EdgeId> map_path, size_t id, size_t rc_id) :
                name_(name),
                src_file_(src_file),
                seq_(seq),
                map_path_(map_path),
                edge_path_(map_path_.simple_path()),
                id_(id),
                rc_id_(rc_id) {    }

    SimpleMappingContig(string name, string src_file, Sequence seq,
            MappingPath<EdgeId> map_path, vector<EdgeId> edge_path,
            size_t id, size_t rc_id) :
                name_(name),
                src_file_(src_file),
                seq_(seq),
                map_path_(map_path),
                edge_path_(edge_path),
                id_(id),
                rc_id_(rc_id) {    }

    Sequence seq() { return seq_; }

    vector<EdgeId> path_seq() { return edge_path_; }

    MappingPath<EdgeId> mapping_path() { return map_path_; }

    string name() { return name_; }

    string src_file() { return src_file_; }

    size_t size() { return edge_path_.size(); }

    size_t length() { return seq_.size(); }

    size_t id(){ return id_; }

    size_t rc_id() { return rc_id_; }

    vector<MappingContigPtr> AllMappingContigs(){
        return vector<MappingContigPtr>();
    }

    void ChangeMappingRange(size_t index, MappingRange new_range){
        VERIFY(index < map_path_.size());
        vector<EdgeId> new_path = map_path_.simple_path();
        vector<MappingRange> new_ranges;
        for(size_t i = 0; i < map_path_.size(); i++)
            if(i != index)
                new_ranges.push_back(map_path_[i].second);
            else
                new_ranges.push_back(new_range);
        MappingPath<EdgeId> new_map_path(new_path, new_ranges);
        map_path_ = new_map_path;
    }

    void ChangeName(string new_name) {
        name_ = new_name;
    }

    string ToString(Graph &graph) {
        stringstream ss;
        ss << "Id: " << id_ << ". Seq size: " << seq_.size() <<
                ". Map path: " << MappingPathToString(graph, map_path_);
        return ss.str();
    }
};

class ReplacedPathMappingContig : public MappingContig{
    MappingContigPtr c_;
    vector<EdgeId> new_path_;
    MappingPath<EdgeId> new_map_path_;

public:
    ReplacedPathMappingContig(MappingContigPtr c, vector<EdgeId> new_path) : c_(c), new_path_(new_path) { }

    ReplacedPathMappingContig(MappingContigPtr c, MappingPath<EdgeId> new_map_path) : c_(c), new_map_path_(new_map_path) {
        new_path_ = new_map_path_.simple_path();
    }

    Sequence seq() { return c_->seq(); }

    vector<EdgeId> path_seq() {
        return new_path_;
    }

    MappingPath<EdgeId> mapping_path(){
        if(new_map_path_.size() != 0)
            return new_map_path_;
        return c_->mapping_path();
    }

    string name() { return c_->name(); }

    string src_file() { return c_->src_file(); }

    size_t size() { return new_path_.size(); }

    size_t length() { return c_->length(); }

    size_t id(){ return c_->id(); }

    size_t rc_id() { return c_->rc_id(); }

    vector<MappingContigPtr> AllMappingContigs(){
        return vector<MappingContigPtr>();
    }

    void ChangeName(string new_name) {
        c_->ChangeName(new_name);
    }

    string ToString(Graph &graph) {
        if(new_map_path_.size() == 0)
            return c_-> ToString(graph);
        stringstream ss;
        ss << "Id: " << id() << ". Seq size: " << seq().size() <<
                ". Map path: " << MappingPathToString(graph, new_map_path_);
        return ss.str();
    }
};

class CompositeMappingContig : public MappingContig{

    Graph &g_;
    size_t k_value_;

    vector<MappingContigPtr> contigs_;
    vector<pair<Range, Range> > overlaps_;

    string contig_name_;

    Sequence composite_seq;
    vector<EdgeId> composite_path;
    size_t composite_size;

    size_t IndexOfEdgeByNumberOfVertex(size_t vertex_index){
        if(vertex_index == 0)
            return 0;
        return vertex_index - 1;
    }

public:
    CompositeMappingContig(Graph &g,
            size_t k_value,
            vector<MappingContigPtr> contigs,
            vector<pair<Range, Range> > overlaps) :
                g_(g),
                k_value_(k_value),
                contigs_(contigs),
                overlaps_(overlaps),
                contig_name_("") {
        VERIFY(contigs.size() > 1);
        VERIFY(contigs.size() == overlaps.size() + 1);
        composite_size = 0;
    }

    Sequence seq(){
        if(composite_seq.size() == 0){
            vector<EdgeId> comp_path = path_seq();
            composite_seq = GetSequenceByPath(g_, k_value_, comp_path);
        }
        return composite_seq;
    }

    vector<EdgeId> path_seq(){
        if(composite_path.size() == 0){
            if(overlaps_.size() == 0){
                if(contigs_.size() == 0)
                    return vector<EdgeId>();
                return contigs_[0]->path_seq();
            }
            else{
                TRACE("New composite contig:");
                TRACE("Path construction of composite contig starts");

                TRACE("Ranges: ");
                for(auto it = overlaps_.begin(); it != overlaps_.end(); it++)
                    TRACE(it->first.start_pos << " - " << it->first.end_pos << ". " <<
                        it->second.start_pos << " - " << it->second.end_pos);

                // first path processing
                {
                    TRACE("First path processing");
                    TRACE("Id - " << contigs_[0]->id());
                    vector<EdgeId> first_path = contigs_[0]->path_seq();
                    size_t end_ind = min<size_t>(IndexOfEdgeByNumberOfVertex(overlaps_[0].first.end_pos),
                            first_path.size() - 1);
                    for(size_t i = 0; i <= end_ind; i++)
                        composite_path.push_back(first_path[i]);
                }

                TRACE("Intermediate paths processing");
                // intermediate paths processing
                for(size_t i = 0; i < overlaps_.size() - 1; i++){
                    auto cur_path = contigs_[i + 1]->path_seq();
                    TRACE("Id: " << contigs_[i + 1]->id());
                    size_t start_ind = min<size_t>(IndexOfEdgeByNumberOfVertex(overlaps_[i].second.end_pos) + 1,
                            cur_path.size() - 1);
                    size_t end_ind = min<size_t>(IndexOfEdgeByNumberOfVertex(overlaps_[i + 1].first.end_pos),
                            cur_path.size() - 1);
                    TRACE("Start - " << start_ind << ", end - " << end_ind);
                    VERIFY(start_ind < cur_path.size() && end_ind < cur_path.size());
                    for(size_t j = start_ind; j <= end_ind; j++)
                        composite_path.push_back(cur_path[j]);
                }

                {
                    // last path processing
                    TRACE("Last path processing");
                    vector<EdgeId> last_path = contigs_[contigs_.size() - 1]->path_seq();
                    TRACE("Id: " << contigs_[contigs_.size() - 1]->id());
                    size_t start_ind = IndexOfEdgeByNumberOfVertex(overlaps_[overlaps_.size() - 1].second.end_pos) + 1;
                    start_ind = min<size_t>(start_ind, last_path.size() - 1);
                    size_t end_ind = last_path.size() - 1;
                    TRACE("Start - " << start_ind << ", end - " << end_ind);
                    VERIFY(start_ind < last_path.size() && end_ind < last_path.size());
                    for(size_t i = start_ind; i <= end_ind; i++)
                        composite_path.push_back(last_path[i]);
                }

                // deletion of repetitive start edge
                TRACE("Deletion of repetitive start-end edge");
                if(composite_path[0] == composite_path[composite_path.size() - 1]){
                    composite_path.erase(composite_path.begin() + composite_path.size() - 1);
                    TRACE("Deletion done");
                }

                TRACE("Path construction of composite contig ends");
            }
        }
        return composite_path;
    }

    MappingPath<EdgeId> mapping_path(){ return MappingPath<EdgeId>(); } // todo refactor

    string name() { return contig_name_; }

    string src_file() { return ""; }

    size_t size(){
        return path_seq().size();
    }

    size_t length() { return seq().size(); }

    size_t id(){ return 0; }

    size_t rc_id() { return 0; }

    void ChangeName(string new_name) {
        contig_name_ = new_name;
    }

    vector<MappingContigPtr> AllMappingContigs(){
        return contigs_;
    }

    string ToString(Graph &){
        return "Composite contig";
    }

private:
    DECL_LOGGER("CompositeMappingContig");
};

class ReplacedNameMappingContig : public MappingContig{
    MappingContigPtr c_;
    string contig_name_;

public:
    ReplacedNameMappingContig(MappingContigPtr c, string contig_name) :
        c_(c),
        contig_name_ (contig_name) { }

    Sequence seq() { return c_->seq(); }

    vector<EdgeId> path_seq() {
        return c_->path_seq();
    }

    MappingPath<EdgeId> mapping_path(){
        return c_->mapping_path();
    }

    string name() { return contig_name_; }

    string src_file() { return c_->src_file(); }

    size_t size() { return c_->size(); }

    size_t length() { return c_->length(); }

    size_t id(){ return c_->id(); }

    size_t rc_id() { return c_->rc_id(); }

    vector<MappingContigPtr> AllMappingContigs(){
        return c_->AllMappingContigs();
    }

    void ChangeName(string new_name) {
        c_->ChangeName(new_name);
    }

    string ToString(Graph &graph) {
        return c_->ToString(graph);
    }
};

}
