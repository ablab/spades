//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include "compare_standard.hpp"

namespace cap {

struct EdgeData {
    std::string sequence;
    TColorSet color;

    EdgeData() {}

    EdgeData(const string& sequence_, const TColorSet& color_)
            : sequence(sequence_),
              color(color_) {}

    bool operator < (const EdgeData &other) const {
        return sequence < other.sequence;
    }
    bool operator == (const EdgeData &other) const {
        return sequence == other.sequence && color == other.color;
    }
};

inline double uniform_01() {
    static boost::mt19937 rng(43);
    static boost::uniform_01<boost::mt19937> zeroone(rng);
    return zeroone();
}

inline bool event_happened(double rate) {
    return math::ls(uniform_01(), rate);
}

inline size_t rand_int(size_t min, size_t max) {
    static boost::mt19937 rng(43);
    boost::uniform_int<size_t> un_int(min, max);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<size_t> > die(rng,
            un_int);
    return die();
}

inline char switch_nucl(char n) {
    VERIFY(is_nucl(n));
    return nucl((dignucl(n) + rand_int(1, 3)) % 4);
}

inline Sequence IntroduceReversal(const Sequence& s, size_t min_len, size_t max_len) {
    VERIFY(s.size() > min_len);
    //inclusive
    size_t start = rand_int(0, s.size() - min_len);
    size_t len = rand_int(min_len, std::min(max_len, s.size() - start));
    //exclusive
    size_t end = start + len;
    INFO(
            "Reversing fragment of length " << len << " from " << start << " to " << end);
    return s.Subseq(0, start) + !s.Subseq(start, end) + s.Subseq(end);
}

inline Sequence IntroduceReversals(const Sequence& s, size_t rev_count, size_t min_len,
        size_t max_len) {
    Sequence res = s;
    for (size_t i = 0; i < rev_count; ++i) {
        res = IntroduceReversal(res, min_len, max_len);
    }
    return res;
}

inline Sequence IntroduceMutations(const Sequence& s, double rate) {
    VERIFY(math::ge(rate, 0.) && math::ls(rate, 1.0));
    string as_str = s.str();
    for (size_t i = 0; i < s.size(); ++i) {
        if (event_happened(rate)) {
            as_str[i] = switch_nucl(as_str[i]);
        }
    }
    return Sequence(as_str);
}

template<class gp_t>
void ConstructRepeatGraph(gp_t& gp) {
    io::RCWrapper<io::SingleRead> rc_stream(make_shared<io::VectorReadStream<io::SingleRead>>(
                                            io::SingleRead("genome", gp.genome.str())));
    ConstructGraph<gp_t::k_value, typename gp_t::graph_t>(gp.g, gp.index,
            rc_stream);
}

template<class Graph>
vector<Sequence> EdgesSequences(const Graph& g) {
    vector<Sequence> res;
    set<EdgeId> edges;
    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        if (edges.find(*it) == edges.end()) {
            res.push_back(g.EdgeNucls(*it));
            edges.insert(g.conjugate(*it));
        }
    }
    return res;
}

template<class gp_t>
vector<Sequence> RepeatGraphEdges(const Sequence& genome) {
    typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;

    gp_t gp(genome);
    ConstructRepeatGraph(gp);
    return EdgesSequences(gp.g);
}

bool CheckFileDiff(const string& file1, const string& file2) {
    INFO("Checking differences between " << file1 << " and " << file2);
    path::CheckFileExistenceFATAL(file1);
    path::CheckFileExistenceFATAL(file2);
    ifstream f1(file1.c_str());
    ifstream f2(file2.c_str());
    while (!f1.eof() && !f2.eof()) {
        string s1;
        f1 >> s1;
        string s2;
        f2 >> s2;
    TRACE("COMPARING " + s1 + " " + s2);
        if (s1 != s2)
            return false;
    }
    if (!f1.eof() || !f2.eof())
        return false;
    return true;
}

void ReadLabelsInMap(ifstream &ifs, std::map <int, EdgeData> &edges) {
    char temp_char;

    bool read_flag = 0;
    int edge_id = 0;
    string sequence_string;

    while (!ifs.eof()) {
        if (read_flag == 0) {
            ifs >> std::ws >> temp_char;
            VERIFY(temp_char == '>');

            ifs >> edge_id;
        } else {
            ifs >> sequence_string;

            EdgeData &node = edges[edge_id];
            node.sequence = sequence_string;
        }

        read_flag ^= 1;
    }
}

void ReadColorsInMap(ifstream &ifs, std::map <int, EdgeData> &edges) {
    size_t nv, ne, temp;
    int edge_id;
    string coloring_string;

    // First, read colors of vertices (we dont need em)
    ifs >> nv;
    for (size_t i = 0; i < 2 * nv; ++i) {
        ifs >> temp;
    }
    // Now that's for edges
    ifs >> ne;
    for (size_t i = 0; i < ne; ++i) {
        ifs >> edge_id >> coloring_string;

        EdgeData &node = edges[edge_id];
        // TColorSet constructor knows about 'uint' legacy data
        node.color = TColorSet(coloring_string);
    }
}

vector <EdgeData> EdgeDataMapSortedValueSet(std::map <int, EdgeData> m) {
    vector <EdgeData> result;
    result.reserve(m.size());
    for (auto it = m.begin(); it != m.end(); ++it) {
        result.push_back(it->second);
    }
    sort(result.begin(), result.end());

    return result;
}

bool MapsValueSetEquals(std::map <int, EdgeData> m1, std::map <int, EdgeData> m2) {
    vector <EdgeData> v1 = EdgeDataMapSortedValueSet(m1),
                      v2 = EdgeDataMapSortedValueSet(m2);

    bool has_errors = false;
    if (v1.size() != v2.size()) {
        INFO("ERROR: Number of edges differ!");
        return false;
    }

    for (size_t i = 0; i < v1.size() && i < v2.size(); ++i) {
        if (v1[i] == v2[i]) continue;

        has_errors = true;
        TRACE(v1[i].sequence << ", " << v1[i].color.ToString() << " --- " << v2[i].sequence << ", " << v2[i].color.ToString());
    }

    if (has_errors) {
        INFO("ERROR: vectors differ!");
        return false;
    }

    return true;
}

template<class gp_t>
inline void LoadWithColoring(gp_t& gp, ColorHandler<typename gp_t::graph_t>& coloring, const string& path) {
    typedef typename debruijn_graph::graphio::ScannerTraits<typename gp_t::graph_t>::Scanner Scanner;
    Scanner scanner(gp.g);
    scanner.loadGraph(path);
    LoadColoring(gp.g, gp.element_finder, coloring, path);
}

template<class gp_t>
class ColoredGraphIsomorphismChecker {

    typedef typename gp_t::graph_t Graph;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef ColorHandler<Graph> Coloring;
    typedef Sequence Kmer;
    typedef map<Kmer, size_t> Mapping;

    struct Pack {
        gp_t gp;
        Coloring col;
        Mapping map;

        Pack(size_t k, const string& work_dir)
                : gp(k, work_dir, 0),
                  col(gp.g) {

        }

        void FillMapping() {
            for (auto it = gp.g.SmartVertexBegin(); !it.IsEnd(); ++it) {
                map[gp.g.VertexNucls(*it)] = (*it).int_id();
            }
        }
    };

    Pack pack1_;
    Pack pack2_;

    void LoadPack(Pack& pack, const string& path) {
        typedef typename debruijn_graph::graphio::ScannerTraits<typename gp_t::graph_t>::Scanner Scanner;
        Scanner scanner(pack.gp.g);
        pack.gp.index.Detach();
        scanner.LoadGraph(path);
        pack.gp.index.Refill();
        pack.gp.index.Attach();
        LoadColoring(pack.gp.g, pack.gp.element_finder, pack.col, path);
        pack.FillMapping();
    }

    VertexId GetVertexId(const Pack& pack, const Kmer& kmer) const {
        return pack.gp.element_finder.ReturnVertexId(get(pack.map, kmer));
    }

    bool CheckEdgeIsomorphism(vector<EdgeData> edges1, vector<EdgeData> edges2) const {
        std::sort(edges1.begin(), edges1.end());
        std::sort(edges2.begin(), edges2.end());
        return edges1 == edges2;
    }

    template<class EdgeContainer>
    vector<EdgeData> ConvertData(const EdgeContainer& edges, const Pack& pack) const {
        vector<EdgeData> ans;
        for (EdgeId e : edges) {
            EdgeData data(pack.gp.g.EdgeNucls(e).str(), pack.col.Color(e));
            ans.push_back(data);
        }
        return ans;
    }

    bool CheckVertexIsomorphism(VertexId v1, VertexId v2) const {
        return CheckEdgeIsomorphism(
                ConvertData(pack1_.gp.g.OutgoingEdges(v1), pack1_),
                ConvertData(pack2_.gp.g.OutgoingEdges(v2), pack2_))
                && CheckEdgeIsomorphism(
                        ConvertData(pack1_.gp.g.IncomingEdges(v1), pack1_),
                        ConvertData(pack2_.gp.g.IncomingEdges(v2), pack2_))
                && pack1_.col.Color(v1) == pack2_.col.Color(v2);
    }

    bool CheckColoredGraphIsomorphism() const {
        typedef typename gp_t::graph_t Graph;
        typedef typename Graph::VertexId VertexId;

        if (pack1_.map.size() != pack2_.map.size())
            return false;

        FOREACH (Kmer kmer, key_set(pack1_.map)) {
            if (pack2_.map.count(kmer) == 0)
                return false;
            if (!CheckVertexIsomorphism(GetVertexId(pack1_, kmer), GetVertexId(pack2_, kmer)))
                return false;
        }
        return true;
    }

public:
    ColoredGraphIsomorphismChecker(size_t k, const string& work_dir)
    : pack1_(k, work_dir),
      pack2_(k, work_dir) {
    }

    bool Check(const string &prefix1, const string &prefix2) {
        INFO("Checking colored graphs in files " + prefix1 + ".* and " + prefix2 + ".* for isomorphic equality");

        LoadPack(pack1_, prefix1);
        LoadPack(pack2_, prefix2);

        return CheckColoredGraphIsomorphism();
    }
};

const io::SingleRead MakeRead(const string& read) {
    //todo fill with good quality
    std::string qual;
    qual.resize(read.size());
    return io::SingleRead("", read, qual);
}

const vector<io::SingleRead> MakeReads(const vector<string>& reads) {
    vector<io::SingleRead> ans;
    for (size_t i = 0; i < reads.size(); ++i) {
        ans.push_back(MakeRead(reads[i]));
    }
    return ans;
}

}
