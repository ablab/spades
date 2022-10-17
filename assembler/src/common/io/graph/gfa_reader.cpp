//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
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

    DeBruijnEdgeData edata(Sequence{record.seq});
    EdgeId e = helper.AddEdge(edata);

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
    std::vector<DeBruijnDataMaster::LinkPtr> empty_links;
    VertexId v1 = helper.CreateVertex(DeBruijnVertexData(empty_links));
    helper.LinkIncomingEdge(v1, e);

    if (e != ce) {
        auto next_storage = new DeBruijnDataMaster::OverlapStorage(empty_links);
        VertexId v2 = helper.CreateVertex(DeBruijnVertexData(next_storage));
        helper.LinkIncomingEdge(v2, ce);
    }
}

void GFAReader::HandleLink(const link &record) {
    link_storage_.push_back(std::make_shared<link>(record));
}

void GFAReader::HandlePath(const gfa::path &record,
                           const io::IdMapper<std::string> &mapper,
                           const ConjugateDeBruijnGraph &g) {
    paths_.emplace_back(std::string{record.name});
    GFAPath &cpath = paths_.back();
    for (const std::string_view &oriented_segment : record.segments) {
        bool rc = oriented_segment.back() == '-';
        std::string_view segment(oriented_segment.data(), oriented_segment.size() - 1);
        EdgeId e = mapper[std::string(segment)];
        if (rc)
            e = g.conjugate(e);
        cpath.edges.push_back(e);
    }
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

    while ((read = gzgetline(&line, &len, fp.get())) != -1) {
        if (read <= 1)
            continue; // skip empty lines

        auto result = gfa::parse_record(line, read - 1);
        if (!result)
            continue;

        std::visit([&](const auto &record) {
            using T = std::decay_t<decltype(record)>;
            if constexpr (std::is_same_v<T, gfa::segment>) {
                num_edges_ += 1;
                HandleSegment(record, *id_mapper, g, helper);
            }
            else if constexpr (std::is_same_v<T, gfa::link>) {
                num_links_ += 1;
                HandleLink(record);
            } else if constexpr (std::is_same_v<T, gfa::path>) {
                HandlePath(record, *id_mapper, g);
            }
        },
            *result);
    }

    unsigned k = ProcessLinks(g, *id_mapper);

    // Add "point tips" of edges
    for (EdgeId e : g.edges()) {
        if (g.EdgeEnd(e))
            continue;

        helper.LinkIncomingEdge(helper.CreateVertex(DeBruijnVertexData()),
                                e);
    }

    free(line);

    if (k != -1U)
        helper.master().set_k(k);

    return k;
}
unsigned GFAReader::ProcessLinks(DeBruijnGraph &g, const io::IdMapper<std::string> &mapper) {
    auto helper = g.GetConstructionHelper();
    unsigned k = -1U;

    for (const auto &link: link_storage_) {
        const auto &record = *link;
        EdgeId e1 = mapper[std::string(record.lhs)], ce1 = g.conjugate(e1);
        if (record.lhs_revcomp)
            std::swap(e1, ce1);

        EdgeId e2 = mapper[std::string(record.rhs)], ce2 = g.conjugate(e2);
        if (record.rhs_revcomp)
            std::swap(e2, ce2);

        const auto &overlap = record.overlap;
        unsigned ovl = -1U;
        if (overlap.size() > 1 ||
            (overlap.size() == 1 && overlap.front().op != 'M')) {
            ovl = std::accumulate(overlap.begin(), overlap.end(), 0, [](size_t sum, const auto &cigar) {
              return sum + cigar.count;
            });
        } else if (overlap.size() == 1) {
            ovl = overlap.front().count;
        }
        if (k == -1U)
            k = ovl;
        else if (k && k != ovl)
            k = 0;

        VertexId v1 = g.EdgeEnd(e1);
        VertexId v2 = g.EdgeStart(e2);
        std::pair<EdgeId, EdgeId> edge_pair(e1, e2);
        auto current_link = std::make_shared<Link>(edge_pair, ovl);
        g.add_link(current_link);
        g.data(v1).add_link(current_link);
        g.data(v1).add_links(g.data(v2).move_links());
        helper.LinkEdges(e1, e2);
    }
    return k;
}

}
