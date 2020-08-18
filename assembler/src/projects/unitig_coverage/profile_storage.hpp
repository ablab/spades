//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/handlers/id_track_handler.hpp"
#include "toolchain/edge_label_helper.hpp"

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
            CHECK_FATAL_ERROR(is, "Unexpected end of stream while reading " << sample_cnt_ << " values");
        }
        return total;
    }

    template<class SingleStream, class Mapper>
    void Fill(SingleStream &reader, size_t stream_id, const Mapper &mapper) {
        typename SingleStream::ReadT read;
        while (!reader.eof()) {
            reader >> read;

            for (const auto &e_mr: mapper.MapSequence(read.sequence())) {
                profiles_[e_mr.first][stream_id] += e_mr.second.mapped_range.size();
            }
        }
    };

public:
    EdgeProfileStorage(const Graph &g, size_t sample_cnt) :
            omnigraph::GraphActionHandler<Graph>(g, "EdgeProfileStorage"),
            sample_cnt_(sample_cnt) {}

    template<class SingleStreamList, class Mapper>
    void Fill(SingleStreamList &streams, const Mapper &mapper) {
        //Initialize profiles
        for (auto it = g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            profiles_[*it] = RawAbundanceVector(sample_cnt_, 0);
        }

#       pragma omp parallel for
        for (size_t i = 0; i < sample_cnt_; ++i) {
            Fill(streams[i], i, mapper);
        }
    }

    size_t sample_cnt() const {
        return sample_cnt_;
    }

    AbundanceVector profile(EdgeId e) const {
        return Normalize(utils::get(profiles_, e), g().length(e));
    }

    void HandleDelete(EdgeId e) override;

    void HandleMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) override;

    void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) override;

    void HandleSplit(EdgeId old_edge, EdgeId new_edge1, EdgeId new_edge2) override;

    void Save(std::ostream &os,
              const io::EdgeNamingF<Graph> &edge_namer = io::IdNamingF<Graph>()) const;

    //TODO maybe pass EdgeDereferenceF?
    void Load(std::istream &is,
              const io::EdgeLabelHelper<Graph> &label_helper,
              bool check_consistency = false);

};

}

}
