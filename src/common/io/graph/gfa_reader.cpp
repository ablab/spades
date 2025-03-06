//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "gfa_reader.hpp"
#include "gfa.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/construction_helper.hpp"

#include "io/utils/id_mapper.hpp"
#include "utils/logger/logger.hpp"
#include "utils/verify.hpp"

#include <string>
#include <memory>
#include <numeric>
#include <string_view>
#include <vector>
#include <unordered_set>

#include <zlib.h>

using namespace debruijn_graph;

namespace gfa {

static ssize_t gzgetdelim(char **buf, size_t *bufsiz, int delimiter, gzFile fp) {
    char *ptr, *eptr;

    if (*buf == NULL || *bufsiz == 0) {
        *bufsiz = BUFSIZ;
        if ((*buf = (char*)malloc(*bufsiz)) == NULL)
            return -1;
    }

    for (ptr = *buf, eptr = *buf + *bufsiz;;) {
        char c = (char)gzgetc(fp);
        if (c == -1) {
            if (gzeof(fp)) {
                ssize_t diff = (ssize_t) (ptr - *buf);
                if (diff != 0) {
                    *ptr = '\0';
                    return diff;
                }
            }
            return -1;
        }
        *ptr++ = c;
        if (c == delimiter) {
            *ptr = '\0';
            return ptr - *buf;
        }
        if (ptr + 2 >= eptr) {
            char *nbuf;
            size_t nbufsiz = *bufsiz * 2;
            ssize_t d = ptr - *buf;
            if ((nbuf = (char*)realloc(*buf, nbufsiz)) == NULL)
                return -1;

            *buf = nbuf;
            *bufsiz = nbufsiz;
            eptr = nbuf + nbufsiz;
            ptr = nbuf + d;
        }
    }
}

static ssize_t gzgetline(char **buf, size_t *bufsiz, gzFile fp) {
    return gzgetdelim(buf, bufsiz, '\n', fp);
}

static void HandleSegment(const gfa::segment &record,
                          io::IdMapper<std::string> &mapper,
                          ConjugateDeBruijnGraph &g,
                          ConjugateDeBruijnGraph::HelperT &helper) {
    unsigned cov = 0;
    if (auto cval = getTag<int64_t>("KC", record.tags))
        cov = unsigned(*cval);

    EdgeId e = helper.AddEdge(DeBruijnEdgeData(Sequence{record.seq}));

    g.coverage_index().SetRawCoverage(e, cov);
    g.coverage_index().SetRawCoverage(g.conjugate(e), cov);

    std::string name{record.name};
    DEBUG("Map ids: " << e.int_id() << ":" << name);
    mapper.map(name, e.int_id());

    EdgeId ce = g.conjugate(e);
    if (e != ce) {
        DEBUG("Map ids: " << ce.int_id() << ":" << name << "'");
        mapper.map(name + '\'', ce.int_id());
    }
}

typedef std::vector<std::tuple<EdgeId, EdgeId, gfa::cigar_string>> Links;

static void HandleLink(Links &links,
                       const gfa::link &record,
                       const io::IdMapper<std::string> &mapper,
                       const ConjugateDeBruijnGraph &g) {
    EdgeId e1 = mapper[std::string(record.lhs)];
    if (record.lhs_revcomp)
        e1 = g.conjugate(e1);

    EdgeId e2 = mapper[std::string(record.rhs)];
    if (record.rhs_revcomp)
        e2 = g.conjugate(e2);

    links.emplace_back(e1, e2, record.overlap);
    if (e1 != g.conjugate(e2))
        links.emplace_back(g.conjugate(e2), g.conjugate(e1), record.overlap);
}

static void HandleGapLink(GFAReader::GapLinks &links,
                          const gfa::gaplink &record,
                          const io::IdMapper<std::string> &mapper,
                          const ConjugateDeBruijnGraph &g) {
    EdgeId e1 = mapper[std::string(record.lhs)];
    if (record.lhs_revcomp)
        e1 = g.conjugate(e1);

    EdgeId e2 = mapper[std::string(record.rhs)];
    if (record.rhs_revcomp)
        e2 = g.conjugate(e2);

    links.emplace_back(e1, e2);
    if (e1 != g.conjugate(e2))
        links.emplace_back(g.conjugate(e2), g.conjugate(e1));
}


static void HandlePath(std::vector<GFAReader::GFAPath> &paths,
                       const gfa::path &record,
                       const io::IdMapper<std::string> &mapper,
                       const ConjugateDeBruijnGraph &g) {
    paths.emplace_back(std::string{record.name});
    auto &cpath = paths.back();
    for (const std::string_view &oriented_segment : record.segments) {
        bool rc = oriented_segment.back() == '-';
        std::string_view segment(oriented_segment.data(), oriented_segment.size() - 1);
        EdgeId e = mapper[std::string(segment)];
        if (rc)
            e = g.conjugate(e);
        cpath.edges.push_back(e);
    }
}

static std::pair<unsigned, bool> ProcessLinks(DeBruijnGraph &g, const Links &links) {
    // First, determine if all overlaps are simple and same (de Bruijn graph, etc.)
    unsigned k = -1U; bool simple = true;
    for (const auto &link: links) {
        const auto &overlap = std::get<2>(link);
        unsigned ovl = -1U;
        if (overlap.size() == 1 && overlap.front().op == 'M') {
            ovl = overlap.front().count;
        } else {
            k = 0; simple = false;
            break;
        }

        if (k == -1U)
            k = ovl;
        else if (k && k != ovl) {
            k = 0; simple = false;
            break;
        }
    }

    INFO("Overlaps: " << (simple ? "simple" : "complex"));

    auto helper = g.GetConstructionHelper();
    for (const auto &link: links) {
        EdgeId e1 = std::get<0>(link), e2 = std::get<1>(link);
        const auto &overlap = std::get<2>(link);

        unsigned ovl = std::accumulate(overlap.begin(), overlap.end(), 0,
                                       [](size_t sum, const auto &cigar) {
                                           return sum + cigar.count;
                                       });

        if (!g.EdgeEnd(e1)) {
            std::vector<LinkId> empty_links;
            helper.LinkIncomingEdge(helper.CreateVertex(DeBruijnVertexData(empty_links)), e1);
        }
        if (!g.EdgeEnd(g.conjugate(e2))) {
            std::vector<LinkId> empty_links;
            helper.LinkOutgoingEdge(helper.CreateVertex(DeBruijnVertexData(empty_links)), e2);
        }

        VertexId v1 = g.EdgeEnd(e1);
        VertexId v2 = g.EdgeStart(e2);
        if (simple) {
            g.set_overlap(v1, ovl);
            g.set_overlap(v2, ovl);
        } else {
            LinkId link_idx = g.add_link(e1, e2, ovl);
            g.add_link(v1, link_idx);
            if (v1 != v2) {
                g.add_links(v1, g.links(v2));
                g.clear_links(v2);
            }
        }

        helper.LinkEdges(e1, e2);
    }

    auto result = std::pair<unsigned, bool>(k, simple);
    return result;
}

unsigned GFAReader::to_graph(ConjugateDeBruijnGraph &g,
                             io::IdMapper<std::string> *id_mapper) {
    num_links_ = num_edges_ = 0; paths_.clear();

    auto helper = g.GetConstructionHelper();
    std::unique_ptr<io::IdMapper<std::string>> local_mapper;
    if (id_mapper == nullptr) {
        local_mapper.reset(new io::IdMapper<std::string>);
        id_mapper = local_mapper.get();
    }

    std::unique_ptr<std::remove_pointer<gzFile>::type, decltype(&gzclose)>
        fp(gzopen(filename_.c_str(), "r"), gzclose);
    if (!fp)
        FATAL_ERROR("Failed to open file: " << filename_);

    char *line = nullptr;
    size_t len = 0;
    ssize_t read;

    Links links;
    while ((read = gzgetline(&line, &len, fp.get())) != -1) {
        if (read <= 1)
            continue; // skip empty lines

        auto result = gfa::parseRecord(line, read - 1);
        if (!result)
            continue;

        std::visit([&](const auto &record) {
            using T = std::decay_t<decltype(record)>;
            if constexpr (std::is_same_v<T, gfa::segment>) {
                num_edges_ += 1;
                HandleSegment(record, *id_mapper, g, helper);
            } else if constexpr (std::is_same_v<T, gfa::link>) {
                num_links_ += 1;
                HandleLink(links, record, *id_mapper, g);
            } else if constexpr (std::is_same_v<T, gfa::path>) {
                HandlePath(paths_, record, *id_mapper, g);
            } else if constexpr (std::is_same_v<T, gfa::gaplink>) {
                HandleGapLink(gap_links_, record, *id_mapper, g);
            }
        },
            *result);
    }

    auto k_and_type = ProcessLinks(g, links);
    unsigned k = k_and_type.first;
    bool is_simple = k_and_type.second;

    // Add "point tips" of edges
    for (EdgeId e: g.edges()) {
        if (g.EdgeEnd(e))
            continue;
        if (is_simple) {
            helper.LinkIncomingEdge(helper.CreateVertex(DeBruijnVertexData(k)), e);
        } else {
            std::vector<LinkId> empty_links;
            helper.LinkIncomingEdge(helper.CreateVertex(DeBruijnVertexData(empty_links)), e);
        }
    }

    free(line);

    // INFO("Filtering dangling vertices");
    for (VertexId v : g.vertices()) {
        if (g.OutgoingEdgeCount(v) > 0 || g.IncomingEdgeCount(v) > 0)
            continue;
        g.DeleteVertex(v);
    }

    if (k != -1U)
        helper.master().set_k(k);

    return k;
}

}
