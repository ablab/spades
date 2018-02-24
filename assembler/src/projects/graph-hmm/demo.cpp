#undef NDEBUG
#include <cassert>
#include "graph.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <queue>
#include <deque>
#include <algorithm>


extern "C" {
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
}

#include <limits>
#include "hmmfile.hpp"


// to compile: gcc this_prog.c -lz
#include <zlib.h>
#include <stdio.h>

#include "utils.hpp"

#include <kseq/kseq.h>
KSEQ_INIT(gzFile, gzread)


#include "demo.hpp"

#include "pathtree.hpp"

// ASL example
// #include "kseq.h"
// KSEQ_INIT(FILE*, fread)
// {
// std::unique_ptr<FILE, void(*)(FILE*)> fp(fopen(""), fclose);
// std::unique_ptr<kseq_t, void(*)(kseq_t*)> seq(kseq_init(fp.get()), kseq_destroy);
// int l;
// size_t reads = 0;
// while ((l = kseq_read(seq.get())) >= 0) {
//   if (++reads % 1000000 == 0)
//     fprintf(stderr, "Processed %zu reads so far\n", reads);
//   std::string header(seq->name.s);
//   std::string sequence(seq->seq.s);
// }
// }

std::string rev_comp(const std::string &s) {
    std::string result;
    std::unordered_map<char, char> rc = { {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'} };
    for (auto it = s.rbegin(), end = s.rend(); it != end; ++it) {
        result += rc[*it];
    }
    return result;
}


Fees levenshtein_fees(const std::string &s, double mismatch, double gap_open, double gap_ext) {
    size_t M = s.size();
    Fees fees;
    fees.M = M;
    DigitalCodind encode;
    fees.code = encode;

    fees.t.resize(M + 1);
    fees.mat.resize(M + 1);
    fees.ins.resize(M + 1);
    for (size_t i = 0; i <= M; ++i) {
        fees.mat[i].resize(4, mismatch);
        fees.ins[i].resize(4, 0);
        fees.t[i].resize(p7H_NTRANSITIONS);
    }

    const double inf = std::numeric_limits<double>::infinity();
    fees.mat[0] = {inf, inf, inf, inf};
    for (size_t i = 1; i <= M; ++i) {
        fees.mat[i][encode(s[i - 1])] = 0;
    }

    // Insertions
    for (size_t i = 0; i <= M; ++i) {
        fees.t[i][p7H_IM] = 0;
        fees.t[i][p7H_II] = gap_ext;
        fees.t[i][p7H_MI] = gap_open;
    }

    // Deletions
    fees.t[0][p7H_MD] = gap_open;
    fees.t[0][p7H_DD] = inf;
    fees.t[0][p7H_DM] = inf;
    for (size_t i = 1; i < M; ++i) {
        fees.t[i][p7H_MD] = gap_open;
        fees.t[i][p7H_DD] = gap_ext;
        fees.t[i][p7H_DM] = 0;
    }
    fees.t[M][p7H_MD] = inf;
    fees.t[M][p7H_DD] = inf;
    fees.t[M][p7H_DM] = 0;

    // Matches
    for (size_t i = 0; i <= M; ++i) {
        fees.t[i][p7H_MM] = 0;
    }

    return fees;
}

Fees hmm_fees(const P7_HMM *hmm, const double lambda = 0) {
    size_t M = hmm->M;
    Fees fees;
    fees.M = M;

    fees.t.resize(M + 1);
    fees.mat.resize(M + 1);
    fees.ins.resize(M + 1);
    for (size_t i = 0; i <= M; ++i) {
        fees.mat[i].resize(4);
        fees.ins[i].resize(4);
        fees.t[i].resize(p7H_NTRANSITIONS);
    }

    for (size_t i = 0; i <= M; ++i) {
        for (size_t j = 0; j < p7H_NTRANSITIONS; ++j) {
            fees.t[i][j] = -log(hmm->t[i][j]);
        }

        for (size_t j = 0; j < 4; ++j) {
            fees.mat[i][j] = -log(hmm->mat[i][j]) + log(0.25) + lambda;
            fees.ins[i][j] = -log(hmm->ins[i][j]) + log(0.25) + lambda;
        }
    }

    return fees;
}


std::vector<std::string> read_fasta_edges(const std::string &filename, bool add_rc) {
    std::vector<std::string> edges;

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(filename.c_str(), "r");
    assert(fp);
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        // printf("name: %s\n", seq->name.s);
        // if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
        // printf("seq: %s\n", seq->seq.s);
        // if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        std::string read = seq->seq.s;

        // U -> T
        for (char &ch : read) {
            if (ch == 'U') {
                ch = 'T';
            }
        }

        edges.push_back(read);
        if (add_rc) {
            edges.push_back(rev_comp(read));
        }
    }

    kseq_destroy(seq);
    gzclose(fp);

    return edges;
}

Fees fees_from_hmm(const P7_HMM *hmm, const ESL_ALPHABET *abc) {
    auto code = DigitalCodind(abc);
    auto fees = hmm_fees(hmm);
    fees.code = code;

    return fees;
}

Fees read_hmm_file(const std::string &filename) {
    hmmer::HMMFile hmm_file(filename);
    auto hmm = hmm_file.read();
    return fees_from_hmm(hmm->get(), hmm->abc());
}

namespace impl {

struct OldTrajectory {
    OldTrajectory() = default;
    explicit OldTrajectory(size_t i) : id{i} {}
    size_t id;
    std::vector<char> turns;
    std::string info;

    OldTrajectory next(char turn = char(-1)) const{
        OldTrajectory traj = *this;
        if (turn != char(-1)) {
            traj.turns.push_back(turn);
        }
        return traj;
    };
};


struct FakeTrajectory {
    template <typename... Args>
    FakeTrajectory(Args...) {}
    auto next(char turn = char(-1)) const {
        return *this;
    }
};

using pathtree::PathLink;

template <typename GraphPointer>
class StateSet : public std::unordered_map<GraphPointer, std::shared_ptr<PathLink<GraphPointer>>> {
public:
    const std::shared_ptr<PathLink<GraphPointer>> &get_or_create(const GraphPointer &key) {
        return this->insert({key, std::make_shared<PathLink<GraphPointer>>()}).first->second;
    }

    StateSet clone() const {
        StateSet copy{*this};
        for (auto &kv : copy) {
            kv.second = kv.second->clone();
        }

        return copy;
    }

    bool update(const GraphPointer& key, double score, GraphPointer from,
                const std::shared_ptr<PathLink<GraphPointer>> &traj) {
        auto it_fl = this->insert({key, std::make_shared<PathLink<GraphPointer>>()});
        double prev = it_fl.second ? std::numeric_limits<double>::infinity() : it_fl.first->second->score();
        // double prev = it_fl.first->second->score();
        it_fl.first->second->update(from, score, traj);
        return prev > score;
    }

    // TODO implement method detecting upd of any kind
};

template <typename GraphPointer>
StateSet<GraphPointer> top_filter(const StateSet<GraphPointer> &S, size_t top, double threshold) {
    using StateSet = StateSet<GraphPointer>;
    top = std::min(top, S.size());
    std::vector<std::pair<typename StateSet::key_type, typename StateSet::mapped_type>> v(S.cbegin(), S.cend());

    // TODO Use std::nth_element
    std::sort(v.begin(), v.end(),
              [](const auto &e1, const auto &e2) { return e1.second->score() < e2.second->score(); });
    if (v.size() < top) {
        top = v.size();
    }
    size_t size = 0;
    while (size < top && v[size].second->score() < threshold) {
        ++size;
    }
    v.resize(size);

    StateSet result;
    result.insert(CONST_ALL(v));
    return result;
}


template <typename GraphPointer>
std::vector<std::pair<std::string, double>> find_best_path(const Fees &fees,
                                                           const std::vector<GraphPointer> &initial) {
    using StateSet = StateSet<GraphPointer>;
    const auto &code = fees.code;

    INFO("pHMM size: " << fees.M);
    if (!fees.check_i_loop(0)) {
        WARN("Negative-cost insertion at the beginning");
    }
    if (!fees.check_i_loop(fees.M)) {
        WARN("Negative-cost insertion at the end");
    }

    for (size_t i = 0; i <= fees.M; ++i) {
        if (!fees.check_i_loop(i)) {
            WARN("Negative-cost insertion at position " << i);
        }
    }

    if (!fees.check_i_negative_loops()) {
        WARN("MODEL CONTAINS NEGATIVE I-LOOPS");
    }

    auto transfer = [&code, &initial](StateSet &to, const StateSet &from,
                                      double transfer_fee,
                                      const std::vector<double> &emission_fees,
                                      const std::string &info = "") {
        assert(&to != &from);
        for (const auto &kv : from) {
            const auto &cur = kv.first;
            const auto &fee = kv.second->score();
            const auto &id = kv.second;
            if (cur.is_empty()) {
                // This branch is used only during BEGIN->M, BEGIN->I and D->M transfers
                for (size_t i = 0; i < initial.size(); ++i) {
                    const auto &next = initial[i];
                    double cost = fee + transfer_fee + emission_fees[code(next.letter())];
                    to.update(next, cost, cur, id);
                }
            } else {
                auto next_pairs = cur.next_pairs();
                for (size_t i = 0; i < next_pairs.size(); ++i) {
                    const auto &next = next_pairs[i].first;
                    char letter = next_pairs[i].second;
                    double cost = fee + transfer_fee + emission_fees[code(letter)];
                    to.update(next, cost, cur, id);
                }
            }
        }
    };
    auto transfer_upd = [&code, &initial](StateSet &to, const StateSet &from,
                                          double transfer_fee,
                                          const std::vector<double> &emission_fees,
                                          const std::string &info = "",
                                          const auto &keys) {
        assert(&to != &from);
        std::unordered_set<GraphPointer> updated;
        for (const auto &cur : keys) {
            auto it = from.find(cur);
            assert(it != from.cend());
            const auto &fee = it->second->score();
            const auto &id = it->second;
            auto next_pairs = cur.next_pairs();
            for (size_t i = 0; i < next_pairs.size(); ++i) {
                const auto &next = next_pairs[i].first;
                char letter = next_pairs[i].second;
                double cost = fee + transfer_fee + emission_fees[code(letter)];
                if (to.update(next, cost, cur, id)) {
                    updated.insert(next);
                }
            }
        }
        return updated;
    };

    auto i_loop_processing = [&transfer_upd, &fees](StateSet &I, size_t m) {
        const size_t max_insertions = 30;

        std::unordered_set<GraphPointer> updated;
        for (const auto &kv : I) {
            updated.insert(kv.first);
        }
        StateSet Inew = I.clone();
        for (size_t i = 0; i < max_insertions; ++i) {
            updated = transfer_upd(Inew, I, fees.t[m][p7H_II], fees.ins[m], "o", updated);
            TRACE(updated.size() << " items updated");
            for (const auto &cur : updated) {
                I[cur] = Inew[cur]->clone();  // TODO Implement minor updation detection
            }
        }
        I = std::move(Inew);  // It is necessary to copy minorly updated states
    };

    auto merge_state_set = [](StateSet &target, const StateSet &source, double transfer_fee = 0) {
        for (const auto &kv : source) {
            const auto &cur = kv.first;
            const auto &id = kv.second;
            target.get_or_create(cur)->merge_update(id.get(), transfer_fee);
        }
    };

    auto dm_new = [&](StateSet &D, StateSet &M, const StateSet &I, size_t m) {
        StateSet Dnew;
        merge_state_set(Dnew, M, fees.t[m - 1][p7H_MD]);
        merge_state_set(Dnew, D, fees.t[m - 1][p7H_DD]);

        StateSet Mnew;
        transfer(Mnew, M, fees.t[m - 1][p7H_MM], fees.mat[m], "m");
        transfer(Mnew, D, fees.t[m - 1][p7H_DM], fees.mat[m], "m");
        transfer(Mnew, I, fees.t[m - 1][p7H_IM], fees.mat[m], "m");

        M = std::move(Mnew);
        D = std::move(Dnew);
    };

    INFO("Initial set size: " << initial.size());

    StateSet I, M, D;
    const auto empty = GraphPointer();
    auto base = PathLink<GraphPointer>::master_source();
    M[empty] = base;  // TODO Implement and use empty Trajectory() instead of Trajectory(0)

    INFO("The number of links (M): " << fees.M);

    transfer(I, M, fees.t[0][p7H_MI], fees.ins[0], "i");
    i_loop_processing(I, 0); // Do we really need I at the beginning???
    for (size_t m = 1; m <= fees.M; ++m) {
        INFO("Step #: " << m);

        dm_new(D, M, I, m);
        I.clear();
        transfer(I, M, fees.t[m][p7H_MI], fees.ins[m], "i");
        i_loop_processing(I, m);

        size_t top = std::max({D.size(), I.size(), M.size()});
        if (m > 10) {
            top = std::min<size_t>(1000000, top);
        }
        if (m > 50) {
            top = std::min<size_t>(10000, top);
        }
        if (m > 500) {
            top = std::min<size_t>(10000, top);
        }

        I = top_filter(I, top, 10);
        M = top_filter(M, top, 10);
        D = top_filter(D, top, 10);
    }

    PathLink<GraphPointer> terminal;
    auto upd_terminal = [&](const StateSet &S, double fee) {
        for (const auto &kv : S) {
            terminal.update(kv.first, kv.second->score() + fee, kv.second);
        }
    };

    upd_terminal(D, fees.t[fees.M][p7H_DM]);
    upd_terminal(I, fees.t[fees.M][p7H_DM]);  // Do we really need I at the end?
    upd_terminal(M, fees.t[fees.M][p7H_MM]);

    INFO("Best score: " << terminal.score());
    INFO("Best of the best");
    INFO(terminal.best_path_string());

    auto result = terminal.top_k_string(10000);
    return result;
}

}  // namespace impl

std::vector<std::pair<std::string, double>> find_best_path_rev(const Fees &fees,
                                                               const std::vector<ReversalGraphPointer<Graph::GraphPointer>> &initial) {
    return impl::find_best_path(fees, initial);
}

std::vector<std::pair<std::string, double>> find_best_path_rev(const Fees &fees,
                                                               const std::vector<ReversalGraphPointer<DBGraph::GraphPointer>> &initial) {
    return impl::find_best_path(fees, initial);
}

std::vector<std::pair<std::string, double>> find_best_path(const Fees &fees, const std::vector<DBGraph::GraphPointer> &initial) {
    return impl::find_best_path(fees, initial);
}

std::vector<std::pair<std::string, double>> find_best_path(const Fees &fees, const std::vector<Graph::GraphPointer> &initial) {
    return impl::find_best_path(fees, initial);
}
