//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/graphio.hpp"
#include "utils/stl_utils.hpp"
#include "assembly_graph/core/graph.hpp"
#include "coordinates_handler.hpp"
#include "math/xmath.h"
#include <iostream>
#include <vector>
#include "utils/logger/logger.hpp"
#include "io/reads/multifile_reader.hpp"
#include "io/reads/splitting_wrapper.hpp"
#include "io/reads/modifying_reader_wrapper.hpp"
#include "io/reads/vector_reader.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace cap {
using namespace debruijn_graph;

template <class Graph>
MappingRange TrivialRange(const Graph& g, typename Graph::EdgeId e, size_t& offset) {
    size_t l = g.length(e);
    offset += l;
    return MappingRange(Range(offset - l, offset), Range(0, l));
}

template <class Graph>
MappingPath<EdgeId> TrivialMappingPath(const Graph& g
        , const vector<typename Graph::EdgeId>& edges) {
  INFO("start tripath");
    vector<MappingRange> ranges;
    size_t offset = 0;
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        ranges.push_back(TrivialRange(g, *it, offset));
    }
  INFO("end tripath");
    return MappingPath<EdgeId>(edges, ranges);
}

inline Sequence ReadSequence(ContigStream& reader) {
    VERIFY(!reader.eof());
    io::SingleRead read;
    reader >> read;
    return read.sequence();
}

template<class Graph, class Index>
void ConstructGraph(Graph& g, Index& index,
        ContigStream& stream) {
    vector<ContigStream*> streams = { &stream };
    ConstructGraph<Graph>(streams, g, index);
}

/*
template<class Graph, class Index>
void ConstructGraph(Graph& g, Index& index,
        ContigStream& stream1,
        ContigStream& stream2) {
    io::MultifileStream<io::SingleRead> composite_reader(stream1, stream2);
    ConstructGraph<Graph, Index>(g, index, composite_reader);
}
*/

inline Sequence ReadGenome(const string& filename) {
    fs::CheckFileExistenceFATAL(filename);
    io::FileReadStream genome_stream(filename);
    return ReadSequence(genome_stream);
}

void WriteGenome(const Sequence& genome, const string& filename) {
  io::osequencestream stream(filename);
  io::SingleRead read("genome", genome.str());
  stream << read;
}

inline vector<io::SingleRead> MakeReads(const vector<Sequence>& ss) {
    vector<io::SingleRead> ans;
    for (size_t i = 0; i < ss.size(); ++i) {
        ans.push_back(io::SingleRead("read_" + std::to_string(i), ss[i].str()));
    }
    return ans;
}

inline Sequence FirstSequence(ContigStream& stream) {
    stream.reset();
    io::SingleRead r;
    VERIFY(!stream.eof());
    stream >> r;
    return r.sequence();
}

inline vector<Sequence> AllSequences(ContigStream& stream) {
    vector<Sequence> answer;
    stream.reset();
    io::SingleRead r;
    while (!stream.eof()) {
        stream >> r;
        answer.push_back(r.sequence());
    }
    return answer;
}

inline vector<Sequence> ReadContigs(const string& filename) {
    fs::CheckFileExistenceFATAL(filename);
    io::FileReadStream genome_stream(filename);
    return AllSequences(genome_stream);
}

//Prints only basic graph structure!!!
//todo rewrite with normal splitter usage instead of filtering
inline void PrintGraphComponentContainingEdge(const string& file_name, const Graph& g,
        size_t split_edge_length, const omnigraph::GraphElementFinder<Graph>& element_finder,
        int int_edge_id) {
    shared_ptr<GraphSplitter<Graph>> inner_splitter = ReliableSplitter<Graph>(g, split_edge_length);

//    VERIFY_MSG(element_finder.ReturnEdgeId(int_edge_id) != NULL,
//            "Couldn't find edge with id = " << int_edge_id);

    shared_ptr<GraphComponentFilter<Graph>> filter = make_shared<AnyEdgeContainFilter<Graph>>(g, element_finder.ReturnEdgeId(int_edge_id));
    FilteringSplitterWrapper<Graph> splitter(inner_splitter, filter);
    vector<vector<VertexId>> components;
    while (splitter.HasNext()) {
        auto component = splitter.Next();
        components.push_back(vector<VertexId>(component.vertices().begin(), component.vertices().end()));
    }
    VERIFY(components.size() == 1);
    debruijn_graph::graphio::ConjugateDataPrinter<Graph> printer(g, components.front().begin(), components.front().end());
    debruijn_graph::graphio::PrintBasicGraph<Graph>(file_name, printer);
}

template<class Graph>
class EdgeCoordinatesGraphLabeler: public visualization::graph_labeler::AbstractGraphLabeler<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
public:
    const CoordinatesHandler<Graph>& edge_pos_;
    const std::vector<std::string> genome_names_;

    EdgeCoordinatesGraphLabeler(const Graph& g,
                              const CoordinatesHandler<Graph>& edge_pos,
                              const std::vector<std::string> &genome_names)
      : AbstractGraphLabeler<Graph>(g),
        edge_pos_(edge_pos),
        genome_names_(genome_names) {
    }

    virtual std::string label(EdgeId edge) const {
    auto ranges = edge_pos_.GetRanges(edge);
    std::sort(ranges.begin(), ranges.end());

    std::stringstream ss;
    for (const auto &entry : ranges) {
      Range genome_range = CoordinatesHandler<Graph>::GetPrintableRange(
              entry.second.first);
      Range seq_range = CoordinatesHandler<Graph>::GetPrintableRange(
              entry.second.second);
      // Make inclusive
      genome_range.end_pos--;
      seq_range.end_pos--;

      ss << genome_names_[size_t(entry.first)] << ": " <<
        "G" << genome_range << ", Seq" << seq_range << "\\n";
    }

        return ss.str();
    }
};

template <class Graph>
class BulgeRemoverCallbackToCoordinatesHandlerAdapter {
 public:
  typedef typename CoordinatesHandler<Graph>::EdgeId EdgeId;

  BulgeRemoverCallbackToCoordinatesHandlerAdapter(
      CoordinatesHandler<Graph> &coordinates_handler)
      : coordinates_handler_(coordinates_handler) {
  }

  void Project(const EdgeId edge_from, const std::vector<EdgeId> &to) {
    std::vector<EdgeId> from;
    from.push_back(edge_from);

    coordinates_handler_.ProjectPath(from, to);

    // Do the same for conjugate sequences as bulge reomver does not provide
    // such functionality :)
    from[0] = coordinates_handler_.GetGraph()->conjugate(from[0]);
    std::vector<EdgeId> to_conj = to;
    std::reverse(to_conj.begin(), to_conj.end());
    for (auto &edge : to_conj)
      edge = coordinates_handler_.GetGraph()->conjugate(edge);

    coordinates_handler_.ProjectPath(from, to_conj);
  }

 private:
  CoordinatesHandler<Graph> &coordinates_handler_;
};

}
