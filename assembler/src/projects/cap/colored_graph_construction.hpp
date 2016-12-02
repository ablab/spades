//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "sequence/runtime_k.hpp"
#include "compare_standard.hpp"
#include "cap_kmer_index.hpp"
#include "modules/graph_construction.hpp"

namespace cap {

template<class Graph, class Mapper>
class CoveredRangesFinder {
    const Graph& g_;
    const Mapper mapper_;

    typedef typename Graph::EdgeId EdgeId;
    typedef restricted::map<EdgeId, vector<Range>> CoveredRanges;

    vector<Range> ProcessRange(Range new_range,
            const vector<Range>& curr_ranges) const {
        vector < Range > answer;
        size_t i = 0;
        while (i < curr_ranges.size()
                && curr_ranges[i].end_pos < new_range.start_pos) {
            answer.push_back(curr_ranges[i]);
            ++i;
        }

        size_t merge_start =
                (i != curr_ranges.size()) ?
                        std::min(curr_ranges[i].start_pos,
                                new_range.start_pos) :
                        new_range.start_pos;

        size_t merge_end = new_range.end_pos;
        while (i < curr_ranges.size()
                && curr_ranges[i].start_pos <= new_range.end_pos) {
            if (curr_ranges[i].end_pos > merge_end)
                merge_end = curr_ranges[i].end_pos;
            ++i;
        }
        answer.push_back(Range(merge_start, merge_end));
        while (i < curr_ranges.size()) {
            answer.push_back(curr_ranges[i]);
            ++i;
        }
        return answer;
    }

    void ProcessPath(const MappingPath<EdgeId>& path,
            CoveredRanges& crs) const {
        for (size_t i = 0; i < path.size(); ++i) {
            auto mapping = path[i];
            EdgeId edge = mapping.first;
            const vector<Range>& curr_ranges = crs[edge];
            Range mapping_range = mapping.second.mapped_range;
            VERIFY(g_.length(edge) >= mapping_range.end_pos);
            crs[edge] = ProcessRange(mapping_range, curr_ranges);
            VERIFY(g_.length(edge) >= crs[edge].back().end_pos);
        }
    }

public:

    CoveredRangesFinder(const Graph& g, const Mapper& mapper) :
            g_(g), mapper_(mapper) {

    }

    void FindCoveredRanges(CoveredRanges& crs, ContigStream& stream) const {
        io::SingleRead read;
        stream.reset();
//        BasicSequenceMapper<gp_t::k_value + 1, Graph> mapper(gp_.g,
//                gp_.index, gp_.kmer_mapper);
        while (!stream.eof()) {
            stream >> read;
            ProcessPath(mapper_.MapSequence(read.sequence()), crs);
        }
    }

};

template<class Graph, class Mapper>
class ColoredGraphConstructor {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef restricted::map<EdgeId, vector<Range>> CoveredRanges;
    typedef restricted::map<EdgeId, vector<size_t>> BreakPoints;

    Graph& g_;
    ColorHandler<Graph>& coloring_;
    const Mapper mapper_;

    void AddBreaks(set<size_t>& breaks, const vector<Range>& ranges) const {
        for (auto it = ranges.begin(); it != ranges.end(); ++it) {
            breaks.insert(it->start_pos);
            breaks.insert(it->end_pos);
        }
    }

    vector<size_t> PostProcessBreaks(const set<size_t>& tmp_breaks) const {
        vector<size_t> breaks(tmp_breaks.begin(), tmp_breaks.end());
        //breaks contain 0 and edge_length here!
        VERIFY(breaks.size() >= 2);
        //cleaning breaks from 0 and edge_length
        vector<size_t> final_breaks;
        for (size_t i = 1; i < breaks.size() - 1; ++i) {
            final_breaks.push_back(breaks[i]);
        }
        return final_breaks;
    }

    void FindBreakPoints(BreakPoints& bps,
            const vector<CoveredRanges>& crss) const {
        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            set<size_t> tmp_breaks;
            for (size_t i = 0; i < crss.size(); ++i) {
                auto crs_it = crss[i].find(e);
                if (crs_it != crss[i].end()) {
                    AddBreaks(tmp_breaks, crs_it->second);
                }
            }
            bps[e] = PostProcessBreaks(tmp_breaks);
            VERIFY(bps[e].empty() || bps[e].back() < g_.length(e));
        }
    }

    void SplitEdge(const vector<size_t>& breaks, EdgeId e) {
        vector<size_t> shifts(breaks.size());
        if (!breaks.empty()) {
            shifts[0] = breaks[0];
            for (size_t i = 1; i < breaks.size(); ++i) {
                shifts[i] = breaks[i] - breaks[i - 1];
            }
        }
        EdgeId curr_e = e;
        for (size_t i = 0; i < breaks.size(); ++i) {
            auto split_result = g_.SplitEdge(curr_e, shifts[i]);
            curr_e = split_result.second;
        }
    }

//    void PrintCRS(const CoveredRanges& crs) {
//        for (auto it = crs.begin(); it != crs.end(); ++it) {
//            DEBUG(
//                    "For edge " << gp_.g.str(it->first) << " ranges "
//                            << it->second);
//        }
//    }

    void SplitGraph(ContigStreams& streams) {
        INFO("Determining covered ranges");
        CoveredRangesFinder<Graph, Mapper> crs_finder(g_, mapper_);
        vector<CoveredRanges> crss(streams.size());
        for (size_t i = 0; i < streams.size(); ++i) {
            crs_finder.FindCoveredRanges(crss[i], streams[i]);
            //        DEBUG("Printing covered ranges for stream i");
            //        PrintCRS(crss[i]);
        }
        BreakPoints bps;
        INFO("Determining breakpoints");
        FindBreakPoints(bps, crss);

        INFO("Splitting graph");
        SplitGraph(bps);
    }

    void SplitGraph(/*const */BreakPoints& bps) {
        vector<EdgeId> initial_edges;
        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            initial_edges.push_back(*it);
        }
        for (auto it = SmartSetIterator<Graph, EdgeId>(g_,
                initial_edges.begin(), initial_edges.end()); !it.IsEnd();
                ++it) {
            EdgeId e = *it;
            VERIFY(bps.find(e) != bps.end());
            VERIFY(bps[e].empty() || bps[e].back() < g_.length(e));
            //todo temporary fix!!!
            if (e == g_.conjugate(e))
                continue;
            SplitEdge(bps[e], e);
        }
    }

    void PaintGraph(ContigStream& stream, TColorSet color) {
        io::SingleRead read;
        stream.reset();
        while (!stream.eof()) {
            stream >> read;
            PaintPath(mapper_.MapSequence(read.sequence()).path(),
                    color);
        }
    }

    void PaintGraph(ContigStreams& streams, const vector<TColorSet>& stream_colors) {
        VERIFY(streams.size() == stream_colors.size());
        for (size_t i = 0; i < streams.size(); ++i) {
            PaintGraph(streams[i], stream_colors[i]);
        }
    }

    void PaintPath(const Path<EdgeId>& path, TColorSet color) {
        for (size_t i = 0; i < path.size(); ++i) {
            coloring_.PaintEdge(path[i], color);
            coloring_.PaintVertex(g_.EdgeStart(path[i]), color);
            coloring_.PaintVertex(g_.EdgeEnd(path[i]), color);
        }
    }

    void CompressGraph(Graph& g, ColorHandler<Graph>& coloring) {
        for (auto it = g.SmartVertexBegin(); !it.IsEnd(); ++it) {
            VertexId v = *it;
            if (g.CanCompressVertex(v)
                    && coloring.Color(g.GetUniqueOutgoingEdge(v))
                            == coloring.Color(g.GetUniqueIncomingEdge(v))) {
                g.CompressVertex(v);
            }
        }
    }

public:
    ColoredGraphConstructor(Graph& g, ColorHandler<Graph>& coloring,
            const Mapper& mapper) :
            g_(g), coloring_(coloring), mapper_(mapper) {

    }

    void ConstructGraph(ContigStreams& streams) {
// It is not truth anymore?
//        VERIFY(streams.size() == 2);

//        if (detailed_output) {
//            //saving for debug and investigation
//            SaveOldGraph(output_folder + "saves/init_graph");
//        }

        SplitGraph(streams);

//        if (detailed_output) {
//            //saving for debug and investigation
//            SaveOldGraph(output_folder + "saves/split_graph");
//        }

        vector<TColorSet> stream_colors;

        // Obsolete two-coloring
//        stream_mapping.push_back(make_pair(streams[0], kRedColorSet));
//        stream_mapping.push_back(make_pair(streams[1], kBlueColorSet));

        TColor color_number = 0;
        for (auto it = streams.begin(); it != streams.end(); ++it) {
            stream_colors.push_back(TColorSet::SingleColor(color_number));
            ++color_number;
        }

        INFO("Coloring graph");
        PaintGraph(streams, stream_colors);
        INFO("Coloring done.");

        //situation in example 6 =)
        INFO("Compressing graph");
        CompressGraph(g_, coloring_);
        INFO("Compressing done.");
    }
};

template<class Graph>
void SimplifyGraph(Graph& g, size_t br_delta) {
    //outdated
    //debruijn_config::simplification::bulge_remover br_config;
    //br_config.max_bulge_length_coefficient = 20;
    //br_config.max_coverage = 1000.;
    //br_config.max_relative_coverage = 1.2;
    //br_config.max_delta = br_delta;
    //br_config.max_relative_delta = 0.1;
    INFO("Removing bulges");
    RemoveBulges(g, br_config);

//        debruijn_config::simplification::tip_clipper tc;
//        tc.max_coverage = 1000;
//        tc.max_relative_coverage = 1000;
//        tc.max_tip_length_coefficient = 6;
//        ClipTips(gp.g, tc, 10 * gp.g.k());
}

template<class gp_t>
void SplitAndColorGraph(gp_t& gp,
        ColorHandler<typename gp_t::graph_t>& coloring,
        ContigStreams& streams) {

    typedef typename gp_t::graph_t Graph;
    typedef typename gp_t::index_t Index;
    typedef NewExtendedSequenceMapper<Graph, Index> Mapper;

    ColoredGraphConstructor<Graph, Mapper> colored_graph_constructor(gp.g,
            coloring, *MapperInstance<gp_t>(gp));


    colored_graph_constructor.ConstructGraph(streams);
}

template<class Graph, class Index, class Streams>
size_t CapConstructGraph(Streams& streams, Graph& g,
        Index& index) {
    return ConstructGraphUsingOldIndex(streams, g, index);
}

template<class gp_t>
void FillPositions(const gp_t &gp, ContigStreams &streams,
    CoordinatesHandler<typename gp_t::graph_t>& coordinates_handler) {
    typedef NewExtendedSequenceMapper<typename gp_t::graph_t,
                                    typename gp_t::index_t> Mapper;

  VERIFY(coordinates_handler.GetGraph() == NULL);
  coordinates_handler.SetGraph(&(gp.g));

  unsigned contig_id = 0;
  std::shared_ptr<const Mapper> mapper = MapperInstance<gp_t>(gp);

  for (auto it = streams.begin(); it != streams.end(); ++it) {
    //cap::RCWrapper stream(**it);
    ContigStream &stream = *it;
    stream.reset();

    io::SingleRead contig;
    // for forward and reverse directions
    while (!stream.eof()) {
      stream >> contig;

      MappingPath<EdgeId> mapping_path = mapper->MapRead(contig);
      const std::vector<EdgeId> edge_path =
          mapping_path.simple_path();
      coordinates_handler.AddGenomePath(contig_id, edge_path);
      contig_id++;
    }

    stream.reset();
  }
}

template<class gp_t>
void ConstructColoredGraph(gp_t& gp,
        ColorHandler<typename gp_t::graph_t>& coloring,
    CoordinatesHandler<typename gp_t::graph_t>& coordinates_handler,
        ContigStreams& streams) {

  INFO("Constructing de Bruijn graph for k=" << gp.k_value);

    CapConstructGraph(streams,
            gp.g, gp.index);
    SplitAndColorGraph(gp, coloring, streams);
  FillPositions(gp, streams, coordinates_handler);
}

//template<class gp_t>
//void ConstructColoredGraph(gp_t& gp,
//        ColorHandler<typename gp_t::graph_t>& coloring,
//        vector<ContigStream*>& streams, const string& reference, bool fill_pos = true, int br_delta = -1) {
//    typedef typename gp_t::graph_t Graph;
//    const size_t k = gp_t::k_value;
//    typedef BasicSequenceMapper<k + 1, Graph> Mapper;
//
//    INFO("Constructing de Bruijn graph for k=" << k);
//
//    //dirty hack because parallel construction uses cfg::get!!!
//    vector<ContigStream*>& tmp_streams(streams.begin(), streams.end());
//    tmp_streams.push_back(EasyC)
//    io::MultifileReader<Contig> stream(tmp_streams);
//    ConstructGraph<k, Graph>(gp.g, gp.index, stream);
//
//    //TODO do we still need it?
//    if (br_delta > 0)
//        SimplifyGraph(gp.g, br_delta);
//
//    ColoredGraphConstructor<Graph, Mapper> colored_graph_constructor(gp.g,
//            coloring, *MapperInstance < gp_t > (gp));
//    colored_graph_constructor.ConstructGraph(streams);
//
//    if (fill_pos) {
//        INFO("Filling contig positions");
//        for (auto it = streams.begin(); it != streams.end(); ++it) {
//            ContigStream& stream = **it;
//            stream.reset();
//            visualization::position_filler::FillPos(gp, stream);
//        }
//    }
//}

}
