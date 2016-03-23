//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "mapping_contig.hpp"

using namespace debruijn_graph;

namespace dipspades {

// interface for contig storage
class ContigStorage{
public:
    virtual void Add(MappingContigPtr new_contig) = 0;
    virtual size_t Size() = 0;
    virtual MappingContigPtr& operator[](size_t index) = 0;
    virtual void ReplaceContig(MappingContigPtr new_contig, size_t index) = 0;
    virtual void DeleteByIDs(set<size_t> ids) = 0;
    virtual MappingContigPtr GetContigById(size_t id) = 0;
    virtual MappingContigPtr GetRCContigById(size_t id) = 0;
    virtual shared_ptr<ContigStorage> Clone() = 0;
    virtual ~ContigStorage(){}

    virtual string ToString(Graph &graph) = 0;
};

typedef shared_ptr<ContigStorage> ContigStoragePtr;

// simple implementation
class SimpleContigStorage : public ContigStorage{
    vector<MappingContigPtr> storage_;

public:
    void Add(MappingContigPtr new_contig) {
        storage_.push_back(new_contig);
    }

    size_t Size() {
        return storage_.size();
    }

    MappingContigPtr& operator[](size_t index){
        VERIFY(index < storage_.size());
        return storage_[index];
    }

    void ReplaceContig(MappingContigPtr new_contig, size_t index) {
        VERIFY(index < storage_.size());
        storage_[index] = new_contig;
    }

    void DeleteByIDs(set<size_t> ids){
        vector<MappingContigPtr> new_storage;
        for(size_t i = 0; i < storage_.size(); i++)
            if(ids.find(storage_[i]->id()) == ids.end())
                new_storage.push_back(storage_[i]);
        storage_ = new_storage;
    }

    MappingContigPtr GetContigById(size_t id){
        for(size_t i = 0; i < storage_.size(); i++)
            if(storage_[i]->id() == id)
                return storage_[i];
        return MappingContigPtr(new SimpleMappingContig());
    }

    MappingContigPtr GetRCContigById(size_t id){
        for(size_t i = 0; i < storage_.size(); i++)
            if(storage_[i]->rc_id() == id)
                return storage_[i];
        return MappingContigPtr(new SimpleMappingContig());
    }

    ContigStoragePtr Clone(){
        ContigStoragePtr clone_storage(new SimpleContigStorage());
        for(size_t i = 0; i < storage_.size(); i++)
            clone_storage->Add(storage_[i]);
        return clone_storage;
    }

    string ToString(Graph &graph) {
        stringstream ss;
        for(auto c = storage_.begin(); c != storage_.end(); c++)
            ss << (*c)->ToString(graph) << endl;
        return ss.str();
    }
};

//-------------------------------------------------------------------------
void save_contig_storage(Graph&g, ContigStoragePtr stor, string fname){

    ofstream save(fname.c_str());
    for(size_t i = 0; i < stor->Size(); i++){
        save << "#" << i << " contig" << endl;
        auto contig = (*stor)[i];
        save << "id " << contig->id() << endl;
        save << "rc_id " << contig->rc_id() << endl;

        auto path = contig->path_seq();
        for(size_t j = 0; j < path.size(); j++){
            save << g.int_id(path[j]) << " ";
        }
        save << endl;
    }
    save.close();
}
//-------------------------------------------------------------------------

}
