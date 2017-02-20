//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/standard_base.hpp"
#include "AlignmentAnalyserNew.hpp"
#include "consistent_mapping.h"

namespace alignment_analysis {
    using omnigraph::MappingRange;

    bool ConsistentMapping::CheckConnect(EdgeId e, Range r) const {
        return CheckConnect(mapped_path.back(), EdgeRange(e, r));
    }

    bool ConsistentMapping::CheckConnect(const EdgeRange &er) const {
        return CheckConnect(mapped_path.back(), er);
    }

    bool ConsistentMapping::CheckConnect(const vector <EdgeRange> &path) const {
        for (size_t i = 0; i + 1 < path.size(); i++) {
            if (!CheckConnect(path[i], path[i + 1]))
                return false;
        }
        return true;
    }

    bool ConsistentMapping::CheckConnect(EdgeId e, MappingRange r) const {
        return CheckConnect(e, r.mapped_range) && initial_range.end_pos == r.initial_range.start_pos;
    }

    bool ConsistentMapping::CheckConnect(const ConsistentMapping &other) const {
        return this->CheckConnect(other.Front()) && initial_range.end_pos == other.initial_range.start_pos;
    }

    void ConsistentMapping::Join(const ConsistentMapping &other) {
        if (!this->IsEmpty() && !other.IsEmpty()) {
            VERIFY(this->initial_range.end_pos == other.initial_range.start_pos);
            VERIFY(CheckConnect(this->mapped_path.back(), other.mapped_path.front()));
        }
        this->initial_range = this->initial_range.Merge(other.initial_range);
        this->mapped_path.insert(this->mapped_path.end(), other.mapped_path.begin(), other.mapped_path.end());
    }

    void ConsistentMapping::Join(const ConsistentMapping &other, const vector <EdgeRange> &path) {
        VERIFY(!this->IsEmpty());
        VERIFY(!other.IsEmpty());
        VERIFY(path.size() == 0 || CheckConnect(this->Back(), path.front()));
        VERIFY(path.size() == 0 || CheckConnect(path.back(), other.Front()));
        VERIFY(CheckConnect(path));
        this->initial_range = this->initial_range.Merge(other.initial_range);
        this->mapped_path.insert(this->mapped_path.end(), path.begin(), path.end());
        this->mapped_path.insert(this->mapped_path.end(), other.mapped_path.begin(), other.mapped_path.end());
    }

    void ConsistentMapping::ForceJoin(const ConsistentMapping &other, const vector <EdgeId> &path) {
        VERIFY(!this->IsEmpty());
        VERIFY(!other.IsEmpty());
        auto pos = other.mapped_path.begin();
        if (path.empty()) {
            VERIFY(graph_.EdgeEnd(this->mapped_path.back().first) == graph_.EdgeEnd(other.mapped_path.front().first));
        } else {
            CutToVertex(graph_.EdgeStart(path.front()));
            while (pos != other.mapped_path.end() && graph_.EdgeEnd(path.back()) != graph_.EdgeStart(pos->first)) {
                ++pos;
            }
            VERIFY(pos != other.mapped_path.end());
        }
        this->mapped_path.back().second.end_pos = graph_.length(this->mapped_path.back().first);
        for (auto it = path.begin(); it != path.end(); ++it) {
            this->mapped_path.push_back(EdgeRange(*it, Range(0, graph_.length(*it))));
        }
        this->initial_range = this->initial_range.Merge(other.initial_range);
        EdgeRange er = *pos;
        er.second.start_pos = 0;
        this->mapped_path.push_back(er);
        this->mapped_path.insert(this->mapped_path.end(), pos + 1, other.mapped_path.end());
    }

    void ConsistentMapping::CutToVertex(VertexId path_start) {
        while (mapped_path.size() > 0 &&
               graph_.EdgeEnd(mapped_path.back().first) != path_start) {
            initial_range.end_pos -= mapped_path.back().second.size();
            mapped_path.pop_back();
        }
        VERIFY(this->mapped_path.size() > 0);
    }

    bool ConsistentMapping::IsEmpty() const {
        return initial_range.empty();
    }

    bool ConsistentMapping::CheckConnect(const EdgeRange &r1, const EdgeRange &r2) const {
        bool result = true;
        if (r1.first != r2.first) {
            result &= graph_.EdgeEnd(r1.first) == graph_.EdgeStart(r2.first);
            result &= r1.second.end_pos == graph_.length(r1.first);
            result &= r2.second.start_pos == 0;
        } else {
            result &= r1.second.end_pos == r2.second.start_pos;
        }
        return result;
    }

    ConsistentMapping::ConsistentMapping(const Graph &graph, const omnigraph::MappingPath<EdgeId> &path) : graph_(graph) {
        VERIFY(path.size() > 0);
        this->initial_range = Range(path.start_pos(), path.end_pos());
        for (size_t i = 0; i < path.size(); i++) {
            VERIFY(i == 0 ||graph.EdgeEnd(path[i - 1].first) == graph.EdgeStart(path[i].first));
            EdgeRange p(path[i].first, path[i].second.mapped_range);
            mapped_path.push_back(p);
        }
        VERIFY(CheckConnect(mapped_path));
    }

    ConsistentMapping::ConsistentMapping(const Graph &graph) : graph_(graph) {
    }

    ConsistentMapping::ConsistentMapping(const Graph &graph, EdgeId e, MappingRange m)
            : graph_(graph), initial_range(m.initial_range), mapped_path{EdgeRange(e, m.mapped_range)} {
    }

    const Range &ConsistentMapping::GetInitialRange() const {
        return initial_range;
    }

    const EdgeRange &ConsistentMapping::Back() const {
        return this->mapped_path.back();
    }

    const EdgeRange &ConsistentMapping::Front() const {
        return this->mapped_path.front();
    }

    Sequence ConsistentMapping::CorrectSequence() const {
        SequenceBuilder sb;
        if(mapped_path.size() == 1) {
            return graph_.EdgeNucls(Front().first).Subseq(Front().second.start_pos, Front().second.end_pos + graph_.k());
        }
        sb.append(graph_.EdgeNucls(Front().first).Subseq(Front().second.start_pos));
        for(auto it = mapped_path.begin(); it != mapped_path.end(); ++it) {
            EdgeId e = it->first;
            Range r = it->second;
            VERIFY(it == mapped_path.begin() || r.start_pos == 0);
            sb.append(graph_.EdgeNucls(e).Subseq(graph_.k(), graph_.k() + r.end_pos));
        }
        return sb.BuildSequence();

    }

    size_t ConsistentMapping::size() const {
        size_t result = 0;
        for(auto it = mapped_path.begin(); it != mapped_path.end(); ++it) {
            result += it->second.size();
        }
        return result;
    }

    VertexId ConsistentMapping::StartVertex() const {
        return graph_.EdgeStart(Front().first);
    }

    VertexId ConsistentMapping::EndVertex() const {
        return graph_.EdgeEnd(Back().first);
    }

    EdgeId ConsistentMapping::StartEdge() const {
        return Front().first;
    }

    EdgeId ConsistentMapping::EndEdge() const {
        return Back().first;
    }

    ConsistentMapping::ConsistentMapping(Graph const &graph, Range r, const vector<EdgeRange> &path) : graph_(graph), initial_range(r), mapped_path(path){

    }

    vector <EdgeRange> ConsistentMapping::GenerateMappingPath(const vector <EdgeId> &path) const {
        vector <EdgeRange> result;
        for(auto it = path.begin(); it != path.end(); ++it) {
            result.push_back(EdgeRange(*it, Range(0, graph_.length(*it))));
        }
        return result;
    }

    const vector <EdgeRange> &ConsistentMapping::GetMappedPath() const {
        return mapped_path;
    }

    string ConsistentMapping::CompareToReference(string const &reference) const {
        string reference_part = reference.substr(initial_range.start_pos, initial_range.size() + graph_.k());
        string correct = CorrectSequence().str();
        if (correct == reference_part) {
            return "Match";
        }
        size_t l = 0;
        size_t r = 0;
        while(l < min(reference_part.size(), correct.size()) && reference_part[l] == correct[l])
            l++;
        while(l + r <= min(reference_part.size(), correct.size()) && reference_part[reference_part.size() - 1 - r] == correct[correct.size() - 1 - r])
            r++;
        stringstream ss;
        if(l + r == reference_part.size()) {
            ss << "Insertion (" << initial_range.start_pos + l << "): Length: " << correct.substr(l, correct.size() - l - r).size();
        } else if(l + r == correct.size()) {
            ss << "Deletion (" << initial_range.start_pos + l + 1 << ", " << initial_range.end_pos - r + graph_.k() << "): Length: " << reference_part.substr(l, reference_part.size() - l - r).size();
        } else {
            ss << "Substitution (" << initial_range.start_pos + l + 1 << ", " <<
                initial_range.end_pos - r + graph_.k() << "): Lengths: " <<
                reference_part.substr(l, reference_part.size() - l - r).size() << " -> " <<
                correct.substr(l, correct.size() - l - r).size();
        }
        return ss.str();
    }

    ostream &operator<<(ostream& os, const EdgeRange& er) {
        os << "EdgeRange(" << er.first.int_id() << " : " << er.second << ")";
        return os;
    }

//    void ConsistentMapping::CloseEnd() {
//        EdgeRange & er = this->mapped_path.back();
//        er.second.end_pos = graph_.length(er.first);
//    }
//
//    void ConsistentMapping::CloseStart() {
//        EdgeRange & er = this->mapped_path.front();
//        er.second.start_pos = 0;
//    }

    ostream &operator<<(ostream& os, const ConsistentMapping& cm) {
        os << cm.GetInitialRange() << " -> ( ";
        for(auto it = cm.GetMappedPath().begin(); it != cm.GetMappedPath().end(); ++it) {
            EdgeRange er = *it;
            os << er << " ";
        }
        os << ")";
        return os;
    }
}
