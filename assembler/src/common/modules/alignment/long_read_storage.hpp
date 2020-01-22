//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/binary/binary.hpp"
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

namespace debruijn_graph {

template<class Graph>
class PathInfo {
public:
    typedef typename Graph::EdgeId EdgeId;

private:
    std::vector<EdgeId> path_;
    mutable size_t w_;
    mutable std::string barcodes_;
    mutable size_t cut_from_begin_;
    mutable size_t cut_from_end_;

public:
    const std::vector<EdgeId> &path() const {
        return path_;
    }

    std::vector<EdgeId> &path() {
        return path_;
    }

    size_t weight() const {
        return w_;
    }

    size_t get_cut_from_begin() const {
        return cut_from_begin_;
    }

    size_t get_cut_from_end() const {
        return cut_from_end_;
    }

    std::string get_barcodes() const {
        return barcodes_;
    }

    void add_barcodes(const std::string &new_barcode) const {
        barcodes_ += "," + new_barcode;
    }

    void increase_weight(int addition = 1) const {
        w_ += addition;
    }

    bool operator<(const PathInfo<Graph> &other) const {
        if (path_ < other.path_) {
            return true;
        } else if (path_ > other.path_) {
            return false;
        }
        if (barcodes_ < other.barcodes_) {
            return true;
        } else if (barcodes_ > other.barcodes_) {
            return false;
        }
        return false;
    }

    PathInfo(const std::vector<EdgeId> &p, size_t weight = 0, const std::string &barcodes = "", size_t cut_from_begin = 0, size_t cut_from_end = 0) :
            path_(p), w_(weight), barcodes_(barcodes), cut_from_begin_(cut_from_begin), cut_from_end_(cut_from_end) {
    }

    PathInfo(const PathInfo<Graph> &other) {
        path_ = other.path_;
        w_ = other.w_;
        barcodes_ = other.barcodes_;
        cut_from_begin_ = other.cut_from_begin_;
        cut_from_end_ = other.cut_from_end_;
    }

    std::string str(const Graph &g_) const {
        std::stringstream s;
        for (EdgeId e : path_) {
            s << g_.int_id(e) << " ";
        }
        return s.str();
    }
};

template<class Graph>
class PathStorage {
    friend class PathInfo<Graph> ;
    typedef typename Graph::EdgeId EdgeId;
    typedef std::map<EdgeId, std::set<PathInfo<Graph>>> InnerIndex;
    using outer_iterator = typename InnerIndex::iterator;
    using inner_iterator = typename std::set<PathInfo<Graph>>::iterator;

    const Graph &g_;
    InnerIndex inner_index_;
    static const size_t kLongEdgeForStats = 500;

    bool IsCanonical(const std::vector<EdgeId> &p) {
        for (size_t i = 0; i < p.size(); ++i) {
            if (p[i].int_id() < g_.conjugate(p[p.size() - 1 - i]).int_id()) {
                return true;
            }
            if (p[i].int_id() > g_.conjugate(p[p.size() - 1 - i]).int_id()) {
                return false;
            }
        }
        return true;
    }

    void HiddenAddPath(const std::vector<EdgeId> &p, int w, const std::string &barcodes, size_t forward_gap, size_t backward_gap) {
        if (p.size() == 0 ) return;
        if (!IsCanonical(p))
            return;
        typename std::set<PathInfo<Graph> >::iterator iter = inner_index_[p[0]].begin();
        for (; iter != inner_index_[p[0]].end(); ++iter) {
            if (iter->path() == p && std::abs((int)forward_gap - (int)iter->get_cut_from_begin()) < 50 && std::abs((int)backward_gap - (int)iter->get_cut_from_end()) < 50) {
                iter->increase_weight(w);
                iter->add_barcodes(barcodes);
                return;
            }
        }
        inner_index_[p[0]].insert(PathInfo<Graph>(p, w, barcodes, forward_gap, backward_gap));
        size_++;
    }

    class flattening_iterator {
    public:
        flattening_iterator() { }
        flattening_iterator(outer_iterator it) : outer_it_(it), outer_end_(it) { }

        flattening_iterator(outer_iterator it, outer_iterator end)
                : outer_it_(it),
                  outer_end_(end)
        {
            if (outer_it_ == outer_end_) { return; }

            inner_it_ = outer_it_->second.begin();
            advance_past_empty_inner_containers();
        }

        PathInfo<Graph> operator*() {
            return *inner_it_;
        }

        flattening_iterator& operator++()
        {
            ++inner_it_;
            if (inner_it_ == outer_it_->second.end())
                advance_past_empty_inner_containers();
            return *this;
        }

        flattening_iterator operator++(int)
        {
            flattening_iterator it(*this);
            ++*this;
            return it;
        }

        friend bool operator==(const flattening_iterator& a,
                               const flattening_iterator& b)
        {
            if (a.outer_it_ != b.outer_it_)
                return false;

            if (a.outer_it_ != a.outer_end_ &&
                b.outer_it_ != b.outer_end_ &&
                a.inner_it_ != b.inner_it_)
                return false;

            return true;
        }

        friend bool operator!=(const flattening_iterator& a,
                               const flattening_iterator& b)
        {
            return !(a == b);
        }
    private:
        void advance_past_empty_inner_containers()
        {
            while (outer_it_ != outer_end_ && inner_it_ == outer_it_->second.end())
            {
                ++outer_it_;
                if (outer_it_ != outer_end_)
                    inner_it_ = outer_it_->second.begin();
            }
        }

        outer_iterator outer_it_;
        outer_iterator outer_end_;
        inner_iterator inner_it_;
    };

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
                this->AddPath(j_iter->path(), (int) j_iter->weight());
            }
        }
    }

    flattening_iterator end() {
        return flattening_iterator(inner_index_.end());
    }

    flattening_iterator begin() {
        return flattening_iterator(inner_index_.begin(), inner_index_.end());
    }

    void ReplaceEdges(std::map<EdgeId, EdgeId> &old_to_new){
        std::map<int, EdgeId> tmp_map;
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
            std::set<PathInfo<Graph>> new_tmp;
            for (auto j_iter = tmp.begin(); j_iter != tmp.end(); j_iter++) {
                PathInfo<Graph> pi = *(j_iter);
                for (size_t k = 0; k < pi.path().size(); k++)
                    if (old_to_new.find(pi.path()[k]) != old_to_new.end()) {
//                        INFO(g_.int_id(old_to_new[pi.path[k]]));
                        pi.path()[k] = old_to_new[pi.path()[k]];
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

    void AddPath(const std::vector<EdgeId> &p, int w, bool add_rc = false, const std::string &barcodes = "", size_t cut_from_begin = 0, size_t cut_from_end = 0) {
        HiddenAddPath(p, w, barcodes, cut_from_begin, cut_from_end);
        if (add_rc) {
            std::vector<EdgeId> rc_p(p.size());
            for (size_t i = 0; i < p.size(); i++)
                rc_p[i] = g_.conjugate(p[p.size() - 1 - i]);
            HiddenAddPath(rc_p, w, barcodes, cut_from_end, cut_from_begin);
        }
    }

    void DumpToFile(const std::string &filename) const{
        std::map<EdgeId, EdgeId> auxilary;
        DumpToFile(filename, auxilary);
    }

    void BinWrite(std::ostream &str) const {
        using io::binary::BinWrite;
        BinWrite(str, inner_index_.size());
        for (const auto &i : inner_index_) {
            BinWrite(str, (size_t)i.second.size());
            for (const auto &j : i.second) {
                BinWrite(str, j.weight());
                BinWrite(str, j.path().size());
                for (const auto &p : j.path()) {
                    BinWrite(str, g_.int_id(p));
                }
            }
        }
    }

    void BinRead(std::istream &str) {
        inner_index_.clear();
        using io::binary::BinRead;

        auto size = BinRead<size_t>(str);
        while (size--) {
            auto count = BinRead<size_t>(str);
            while (count--) {
                auto weight = BinRead<size_t>(str);
                auto length = BinRead<size_t>(str);
                std::vector<EdgeId> path;
                path.reserve(length);
                while (length--) {
                    auto eid = BinRead<uint64_t>(str);
                    path.push_back(eid);
                }
                AddPath(path, (int)weight);
            }
        }
    }

    void DumpToFile(const std::string& filename, const std::map<EdgeId, EdgeId>& replacement,
                    size_t stats_weight_cutoff = 1, bool need_log = false) const {
        std::ofstream filestr(filename);
        std::set<EdgeId> continued_edges;

        for(auto iter = inner_index_.begin(); iter != inner_index_.end(); ++iter){
            filestr << iter->second.size() << std::endl;
            int non1 = 0;
            for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
                filestr << " Weight: " << j_iter->weight();
                if (j_iter->weight() > stats_weight_cutoff)
                    non1++;

                filestr << " length: " << j_iter->path().size() << " ";
                for (auto p_iter = j_iter->path().begin(); p_iter != j_iter->path().end(); ++p_iter) {
                    if (p_iter != j_iter->path().end() - 1 && j_iter->weight() > stats_weight_cutoff) {
                        continued_edges.insert(*p_iter);
                    }

                    filestr << g_.int_id(*p_iter) << "(" << g_.length(*p_iter) << ") ";
                }
                filestr << std::endl;
            }
            filestr << std::endl;
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

    void SaveAllPaths(std::vector<PathInfo<Graph>> &res) const {
        for (auto iter = inner_index_.begin(); iter != inner_index_.end(); ++iter) {
            for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
                res.push_back(*j_iter);
            }
        }
    }

    void LoadFromFile(const std::string &s, bool force_exists = true) {
        FILE* file = fopen(s.c_str(), "r");
        if (force_exists) {
            VERIFY(file != NULL);
        } else if (file == NULL) {
            INFO("Long reads not found, skipping");
            return;
        }
        fclose(file);

        INFO("Loading long reads alignment...");
        std::ifstream filestr(s);
        INFO("loading from " << s);
        std::map<size_t, EdgeId> tmp_map;
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
                std::vector<EdgeId> p;
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

    void AddStorage(PathStorage<Graph> &to_add) {
        for (auto iter = to_add.inner_index_.begin(); iter != to_add.inner_index_.end(); iter++) {
            for(auto j_iter = iter->second.begin(); j_iter != iter->second.end(); j_iter ++) {
                this->AddPath(j_iter->path(), (int) j_iter->weight());
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
    std::vector<PathStorage<Graph>> data_;

public:
    typedef PathStorage<Graph> value_type;

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

} // namespace debruijn_graph
