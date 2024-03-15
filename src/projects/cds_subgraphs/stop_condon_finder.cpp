//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "stop_condon_finder.hpp"
#include "assembly_graph/paths/path_utils.hpp"

namespace cds_subgraphs {

static CodonSet RC(const CodonSet &codons) {
    CodonSet rc(codons.size());
    std::transform(codons.begin(), codons.end(), rc.begin(),
                   [](const Sequence &s) {return !s;});
    return rc;
}

CodonSet STOP_CODONS = {Sequence("TAG"), Sequence("TAA"), Sequence("TGA")};
CodonSet RC_STOP_CODONS = RC(STOP_CODONS);

void CodonFinder::NextToQueue(FramedPos fpos) {
    if (fpos.last(g_)) {
        for (EdgeId e : g_.OutgoingEdges(g_.EdgeEnd(fpos.e))) {
            queue_.push(std::make_pair(fpos.next(e), fpos));
        }
    } else {
        queue_.push(std::make_pair(fpos.next(), fpos));
    }
}

std::unordered_map<GraphPos, size_t> CodonFinder::Terminates(const std::vector<EdgePath> &paths) const {
    std::unordered_map<GraphPos, size_t> pos_len;
    for (const auto &p : paths) {
        pos_len[GraphPos(p.sequence().back(), p.end_pos())] = PathLength(g_, p);
    }
    return pos_len;
};

std::set<EdgeId> CodonFinder::Edges(const std::vector<EdgePath> &paths) const {
    std::set<EdgeId> answer;
    for (const auto &p : paths) {
        utils::insert_all(answer, p.sequence());
    }
    return answer;
}

std::vector<EdgePath> CodonFinder::Go() {
    FramedPos init(init_pos_.first, init_pos_.second);
    queue_.push(std::make_pair(init, FramedPos()));

    std::vector<FramedPos> terminated;

    while (!queue_.empty()) {
        FramedPos fpos;
        FramedPos prev;
        std::tie(fpos, prev) = queue_.front();
        queue_.pop();
        if (!prev_.count(fpos)) {
            prev_[fpos] = prev;
            if (Terminate(fpos)) {
                DEBUG("Terminate graph pos: " << g_.str(fpos.e) << " " << fpos.offset);
                DEBUG("Codon start coord " << fpos.offset + g_.k() - 1);
                DEBUG("Codon " << Codon(GraphPos(fpos.e, fpos.offset)));
                terminated.push_back(fpos);
            } else {
                NextToQueue(fpos);
            }
        }
    }

    std::vector<EdgePath> paths;

    for (FramedPos t : terminated) {
        std::vector<EdgeId> reverse_path;
        reverse_path.push_back(t.e);
        if (t == init) {
            paths.push_back(EdgePath(std::vector<EdgeId>(reverse_path.rbegin(), reverse_path.rend()),
                                     init.offset + 1, init.offset + 1));
        } else {
            FramedPos fpos = t;
            while (true) {
                VERIFY(prev_.count(fpos));
                auto prev = prev_[fpos];
                if (prev == init) {
                    break;
                }
                if (prev.e != fpos.e || fpos.offset != prev.offset + 1) {
                    reverse_path.push_back(prev.e);
                }
                fpos = prev;
            }
            paths.push_back(EdgePath(std::vector<EdgeId>(reverse_path.rbegin(), reverse_path.rend()),
                                     fpos.offset, t.offset + 1));
        }
    }

    return paths;
}
}
