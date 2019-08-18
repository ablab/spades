//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/handlers/id_track_handler.hpp"

#include <vector>

namespace debruijn_graph {

namespace coverage_profiles {

//TODO always working with doubles seems easier and more correct
class EdgeProfileStorage : public omnigraph::GraphActionHandler<Graph> {
    typedef Graph::EdgeId EdgeId;
    typedef Graph::VertexId VertexId;
    typedef std::vector<size_t> RawAbundanceVector;
    typedef std::vector<double> AbundanceVector;

    size_t sample_cnt_;
    std::unordered_map<EdgeId, RawAbundanceVector> profiles_;

    // FIXME self-conjugate edge coverage?!
    template<class SingleStream, class Mapper>
    void Fill(SingleStream &reader, size_t stream_id, const Mapper &mapper) {
        typename SingleStream::ReadT read;
        while (!reader.eof()) {
            reader >> read;
            //TRACE("Aligning " << read.name());

            for (const auto &e_mr: mapper.MapSequence(read.sequence())) {
                // FIXME: should initial_range be used instead?
                profiles_[e_mr.first][stream_id] += e_mr.second.mapped_range.size();
            }
        }
    };

    AbundanceVector Normalize(const RawAbundanceVector &p, size_t length) const {
        AbundanceVector answer(sample_cnt_);
        for (size_t i = 0; i < sample_cnt_; ++i) {
            answer[i] = double(p[i]) / double(length);
        }
        return answer;
    }

    void Add(RawAbundanceVector &p, const RawAbundanceVector &to_add) const {
        for (size_t i = 0; i < sample_cnt_; ++i) {
            p[i] += to_add[i];
        }
    }

    RawAbundanceVector MultiplyEscapeZero(const AbundanceVector &p, size_t factor) const {
        RawAbundanceVector answer(sample_cnt_);
        for (size_t i = 0; i < sample_cnt_; ++i) {
            answer[i] = size_t(math::round(p[i] * double(factor)));
            if (answer[i] == 0 && math::gr(p[i], 0.))
                answer[i] = 1;
        }
        return answer;
    }

    AbundanceVector LoadAbundanceVector(std::istream &is) const {
        AbundanceVector total(sample_cnt_, 0.);
        for (size_t i = 0; i < sample_cnt_; ++i) {
            is >> total[i];
            VERIFY_MSG(is, "Unexpected end of stream while reading " << sample_cnt_ << " values");
        }
        return total;
    }

public:
    EdgeProfileStorage(const Graph &g, size_t sample_cnt) :
            omnigraph::GraphActionHandler<Graph>(g, "EdgeProfileStorage"),
            sample_cnt_(sample_cnt) {}

    void HandleDelete(EdgeId e) override {
        profiles_.erase(e);
    }

    size_t sample_cnt() const {
        return sample_cnt_;
    }

    void HandleMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) override {
        RawAbundanceVector total(sample_cnt_, 0);
        for (EdgeId e : old_edges) {
            Add(total, utils::get(profiles_, e));
        }
        profiles_[new_edge] = std::move(total);
    }

    void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) override {
        RawAbundanceVector total(utils::get(profiles_, edge1));
        Add(total, utils::get(profiles_, edge2));
        profiles_[new_edge] = std::move(total);
    }

    void HandleSplit(EdgeId old_edge, EdgeId new_edge1, EdgeId new_edge2) override {
        AbundanceVector abund = profile(old_edge);
        if (old_edge == g().conjugate(old_edge)) {
            RawAbundanceVector raw1 = MultiplyEscapeZero(abund, g().length(new_edge1));
            profiles_[new_edge1] = raw1;
            profiles_[g().conjugate(new_edge1)] = std::move(raw1);
            profiles_[new_edge2] = MultiplyEscapeZero(abund, g().length(new_edge2));
        } else {
            profiles_[new_edge1] = MultiplyEscapeZero(abund, g().length(new_edge1));
            profiles_[new_edge2] = MultiplyEscapeZero(abund, g().length(new_edge2));
        }
    }

    template<class SingleStreamList, class Mapper>
    void Fill(SingleStreamList &streams, const Mapper &mapper) {
        for (auto it = g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            profiles_[*it] = RawAbundanceVector(sample_cnt_, 0);
        }

#       pragma omp parallel for
        for (size_t i = 0; i < sample_cnt_; ++i) {
            Fill(streams[i], i, mapper);
        }
    }

    AbundanceVector profile(EdgeId e) const {
        return Normalize(utils::get(profiles_, e), g().length(e));
    }

    //FIXME support namer (in particular IdMapper)
    void Save(std::ostream &os) const {
        for (auto it = g().ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            os << g().int_id(e) << '\t';
            auto prof = profile(e);
            std::copy(prof.begin(), prof.end(), std::ostream_iterator<double>(os, "\t"));
            os << '\n';
        }
    }

    void Load(std::istream &is, const omnigraph::GraphElementFinder<Graph> &element_finder,
            std::unique_ptr<std::unordered_map<std::string, size_t>> label2int_id = nullptr) {
        std::string s;
        while (std::getline(is, s)) {
            std::istringstream ss(s);
            std::string label;
            ss >> label;
            size_t int_id;
            if (!label2int_id) {
                int_id = std::stoi(label);
            } else {
                int_id = utils::get(*label2int_id, label);
            }
            EdgeId e = element_finder.ReturnEdgeId(int_id);
            auto p = MultiplyEscapeZero(LoadAbundanceVector(ss), g().length(e));
            profiles_[e] = p;
            profiles_[g().conjugate(e)] = std::move(p);
        }

        //TODO make safeguard check that all the profiles loaded successfully optional?
        //VERIFY(label2int_id);
        //std::map<size_t, std::string> i2l;
        //for (auto k_v : (*label2int_id)) {
        //    i2l[k_v.second] = k_v.first;
        //}
        for (auto it = g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            VERIFY_MSG(profiles_.count(e) > 0, "Failed to load profile for one of the edges");
            //label " << i2l[g().int_id(e)] << ".");
        }
    }

};

//void SaveStorage(const EdgeProfileStorage &storage, std::ostream &os) {
//    os << storage.sample_cnt() << 'n';
//    storage.Save(os);
//}
//
//EdgeProfileStorage LoadStorage(const Graph &g,
//                               const omnigraph::GraphElementFinder<Graph> &element_finder,
//                               std::istream &is) {
//    size_t sample_cnt;
//    is >> sample_cnt;
//    EdgeProfileStorage storage(g, sample_cnt);
//    storage.Load(is, element_finder);
//    return storage;
//}

}

}
