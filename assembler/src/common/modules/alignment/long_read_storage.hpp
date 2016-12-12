//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * long_edge_storage.hpp
 *
 *  Created on: Feb 7, 2013
 *      Author: lab42
 */

#pragma once

#include <algorithm>

namespace debruijn_graph {

template<class Graph>
class PathInfo {
public:
    typedef typename Graph::EdgeId EdgeId;
    vector<EdgeId> path;

private:
    mutable size_t w;

public:
    const vector<EdgeId>& getPath() const {
        return path;
    }

    size_t getWeight() const {
        return w;
    }

    void increaseWeight(int addition = 1) const {
        w += addition;
    }

    bool operator<(const PathInfo<Graph> &other) const {
        return path < other.path;
    }

    PathInfo(const vector<EdgeId> &p, size_t weight = 0) :
            path(p), w(weight) {
    }

    PathInfo(const PathInfo<Graph> &other) {
        path = other.path;
        w = other.w;
    }

    string str(const Graph &g_) const {
        stringstream s;
        for(auto iter = path.begin(); iter != path.end(); iter ++ ){
            s << g_.int_id(*iter) << " ";
        }
        return s.str();
    }

};

template<class Graph>
class PathStorage {
    friend class PathInfo<Graph> ;
    typedef typename Graph::EdgeId EdgeId;
    typedef map<EdgeId, set<PathInfo<Graph> > > InnerIndex;

    const Graph &g_;
    InnerIndex inner_index_;
    static const size_t kLongEdgeForStats = 500;

    void HiddenAddPath(const vector<EdgeId> &p, int w){
        if (p.size() == 0 ) return;
        for (typename set<PathInfo<Graph> >::iterator iter = inner_index_[p[0]].begin(); iter != inner_index_[p[0]].end(); ++iter) {

            if (iter->path == p) {
                iter->increaseWeight(w);
                return;
            }
        }
        inner_index_[p[0]].insert(PathInfo<Graph>(p, w));
        size_++;
    }

public:

    PathStorage(const Graph &g)
            : g_(g),
              inner_index_(),
              size_(0) {
    }

    PathStorage(const PathStorage & p)
            : g_(p.g_),
              inner_index_(),
              size_(0) {
        for (auto iter = p.inner_index_.begin(); iter != p.inner_index_.end();
                iter++) {
            for (auto j_iter = iter->second.begin();
                    j_iter != iter->second.end(); j_iter++) {
                this->AddPath(j_iter->path, (int) j_iter->getWeight());
            }
        }
    }

    void ReplaceEdges(map<EdgeId, EdgeId> &old_to_new){
        map<int, EdgeId> tmp_map;
//        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter ){
//            tmp_map[g_.int_id(*iter)] = *iter;
//        }
        InnerIndex new_index;
        for (auto iter = inner_index_.begin(); iter != inner_index_.end(); iter++) {
            auto tmp = iter->second;
            EdgeId new_first;
            if (old_to_new.find(iter->first) == old_to_new.end())
                new_first = iter->first;
            else {
                DEBUG("new first edge: "<< g_.int_id(old_to_new[iter->first]) << " with " << tmp.size() << " edges ");
                new_first = old_to_new[iter->first];
            }
            set<PathInfo<Graph> > new_tmp;
            for (auto j_iter = tmp.begin(); j_iter != tmp.end(); j_iter++) {
                PathInfo<Graph> pi = *(j_iter);
                for (size_t k = 0; k < pi.path.size(); k++)
                    if (old_to_new.find(pi.path[k]) != old_to_new.end()) {
//                        INFO(g_.int_id(old_to_new[pi.path[k]]));
                        pi.path[k] = old_to_new[pi.path[k]];
                    }
                DEBUG(pi.str(g_));
                new_tmp.insert(pi);

            }
            if (new_first != iter->first) {
                TRACE("and mmew_tmp.size: "<< new_tmp.size());
            }
            if (new_index.find(new_first) == new_index.end()) {
                new_index[new_first] = new_tmp;
            } else {
                for (auto j_iter = new_tmp.begin(); j_iter != new_tmp.end(); j_iter++) {
                    new_index[new_first].insert(*j_iter);
                }
            }

        }

        inner_index_ = new_index;
    }

    void AddPath(const vector<EdgeId> &p, int w, bool add_rc = false) {
        HiddenAddPath(p, w);
        if (add_rc) {
            vector<EdgeId> rc_p(p.size());
            for (size_t i = 0; i < p.size(); i++)
                rc_p[i] = g_.conjugate(p[p.size() - 1 - i]);
            HiddenAddPath(rc_p, w);
        }
    }

    void DumpToFile(const string& filename) const{
        map <EdgeId, EdgeId> auxilary;
        DumpToFile(filename, auxilary);
    }

    void DumpToFile(const string& filename, const map<EdgeId, EdgeId>& replacement,
                    size_t stats_weight_cutoff = 1, bool need_log = false) const {
        ofstream filestr(filename);
        set<EdgeId> continued_edges;

        for(auto iter = inner_index_.begin(); iter != inner_index_.end(); ++iter){
            filestr<< iter->second.size() << endl;
            int non1 = 0;
            for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
                filestr << " Weight: " << j_iter->getWeight();
                if (j_iter->getWeight() > stats_weight_cutoff)
                    non1++;

                filestr << " length: " << j_iter->path.size() << " ";
                for (auto p_iter = j_iter->path.begin(); p_iter != j_iter->path.end(); ++p_iter) {
                    if (p_iter != j_iter->path.end() - 1 && j_iter->getWeight() > stats_weight_cutoff) {
                        continued_edges.insert(*p_iter);
                    }

                    filestr << g_.int_id(*p_iter) << "(" << g_.length(*p_iter) << ") ";
                }
                filestr << endl;
            }
            filestr << endl;
        }

        int noncontinued = 0;
        int long_gapped = 0;
        int continued = 0;
        if (need_log) {
            for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
                EdgeId e = *iter;
                if (g_.length(e) > kLongEdgeForStats) {
                    if (!g_.IsDeadEnd(g_.EdgeEnd(e))) {
                        if (continued_edges.find(e) == continued_edges.end()) {
                            auto replacement_it = replacement.find(e);
                            if (replacement_it != replacement.end() &&
                                continued_edges.find(replacement_it->second) != continued_edges.end()) {
                                TRACE("found in teplacement, edges " << g_.int_id(e) << " " <<
                                      g_.int_id(replacement_it->second) << " skipping ");
                                continue;
                            }
                            TRACE("noncontinued end left " << g_.int_id(e));
                            noncontinued++;
                        } else
                            continued++;
                    } else {
                        TRACE("dead end left " << g_.int_id(e));
                        long_gapped++;
                    }
                }
            }
            INFO("After PacBio (long reads) aligning, for edges longer than " << kLongEdgeForStats << ":");
            INFO("No continuation found for " << noncontinued + long_gapped << " edges of " <<
                 noncontinued + continued + long_gapped);
        }
    }

     void SaveAllPaths(vector<PathInfo<Graph>> &res) const {
        for (auto iter = inner_index_.begin(); iter != inner_index_.end(); ++iter) {
            for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
                res.push_back(*j_iter);
            }
        }
    }

    void LoadFromFile(const string s, bool force_exists = true) {
        FILE* file = fopen(s.c_str(), "r");
        if (force_exists) {
            VERIFY(file != NULL);
        } else if (file == NULL) {
            INFO("Long reads not found, skipping");
            return;
        }
        fclose(file);

        INFO("Loading long reads alignment...");
        ifstream filestr(s);
        INFO("loading from " << s);
        map<size_t, EdgeId> tmp_map;
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            tmp_map[g_.int_id(*iter)] = *iter;
        }
        int fl;

        file = fopen((s).c_str(), "r");
        char ss[14];
        while (!feof(file)) {
            int n;

            fl = fscanf(file, "%d\n", &n);
            if (fl != 1)
                break;
            TRACE(n);
            for (int i = 0; i < n; i++) {

                int w = -1, l = -1;
                fl = fscanf(file, "Weight: %d length: %d", &w, &l);
                TRACE(w << " " << l);
                VERIFY(fl == 2);
                vector<EdgeId> p;
                for (int j = 0; j < l; j++) {
                    size_t e;
                    int x;
                    fl = fscanf(file, "%zu(%d)", &e, &x);
                    VERIFY(fl == 2);
                    VERIFY(tmp_map.find(e) != tmp_map.end());
                    p.push_back(tmp_map[e]);
                }
                fl = fscanf(file, "%[^\n]\n", ss);
                TRACE(ss[0]);
                AddPath(p, w);
            }
        }
        fclose(file);
        INFO("Loading finished.");
    }

    void AddStorage(PathStorage<Graph> & to_add) {

        for(auto iter = to_add.inner_index_.begin(); iter != to_add.inner_index_.end(); iter++) {
            for(auto j_iter = iter->second.begin(); j_iter != iter->second.end(); j_iter ++) {
                this->AddPath(j_iter->path, (int) j_iter->getWeight());
            }
        }
    }

    void Clear() {
        inner_index_.clear();
        size_ = 0;
    }

    size_t size() const {
        return size_;
    }

//    typename InnerIndex::iterator begin() const {
//        return inner_index.begin();
//    }
//
//    typename InnerIndex::iterator end() const {
//        return inner_index.end();
//    }
//    typename InnerIndex::iterator operator*(){
//        return this->first;
//    }
private:
    size_t size_;
};

template<class Graph>
class LongReadContainer {
    Graph& g_;
    vector<PathStorage<Graph>> data_;

public:

    LongReadContainer(Graph& g, size_t count = 0): g_(g) {
        for (size_t i = 0; i < count; ++i) {
            data_.emplace_back(g_);
        }
    }

    PathStorage<Graph>& operator[](size_t index) {
        return data_[index];
    }

    const PathStorage<Graph>& operator[](size_t index) const {
        return data_[index];
    }

    size_t size() const {
        return data_.size();
    }

    void Clear() {
        for (auto& storage : data_) {
            storage.Clear();
        }
    }

};


}


