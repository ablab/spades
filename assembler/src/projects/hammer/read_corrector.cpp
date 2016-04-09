//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "read_corrector.hpp"

#include "kmer_data.hpp"
#include "kmer_stat.hpp"
#include "valid_kmer_generator.hpp"

#include <string>
#include <vector>
#include <queue>

using namespace hammer;

using positions_t = std::array<uint16_t, 4>;

struct state {
    state(size_t p, std::string s, double pen, KMer l, positions_t c)
            : pos(p), str(s), penalty(pen), last(l), cpos(c) {}

    size_t pos; std::string str; double penalty; KMer last; positions_t cpos;
};

std::ostream& operator<<(std::ostream &os, const state &state) {
    os << "[pos: " << state.pos << ", last: " << state.last << " penalty: " << state.penalty << "]";
    return os;
}

std::ostream& operator<<(std::ostream &os, const positions_t &pos) {
    os << "[" << pos[0] << ", " << pos[1] << ", " << pos[2] << ", " << pos[3] << "]";
    return os;
}

namespace std {
template<>
struct less<state> {
    bool operator()(const state &lhs, const state &rhs) const {
        return lhs.penalty < rhs.penalty ||
               (lhs.penalty == rhs.penalty && lhs.pos < rhs.pos);
    }
};
};

static void FlushCandidates(std::priority_queue<state>& corrections, std::priority_queue<state> &candidates,
                            size_t size_limit) {
    if (candidates.empty())
        return;

    if (corrections.size() > size_limit) {
        corrections.emplace(candidates.top());
    } else {
        while (!candidates.empty()) {
            corrections.emplace(candidates.top());
            candidates.pop();
        }
    }

    std::priority_queue<state>().swap(candidates);
}

std::string ReadCorrector::CorrectReadRight(const std::string &seq, const std::string &qual,
                                            size_t right_pos) {
    const size_t read_size = seq.size();
    std::priority_queue<state> corrections, candidates;
    positions_t cpos{{(uint16_t)-1, (uint16_t)-1U, (uint16_t)-1U, (uint16_t)-1U}};

    const size_t size_thr = size_t(100 * log2(read_size - right_pos)) + 1;
    const double penalty_thr = -(double)(read_size - right_pos) * 15.0 / 100;
    const size_t pos_thr = 8;

    corrections.emplace(right_pos, seq,
                        0.0, KMer(seq, right_pos - K + 1, K, /* raw */ true),
                        cpos);
    while (!corrections.empty()) {
        state correction = corrections.top(); corrections.pop();
        size_t pos = correction.pos + 1;
        if (pos == read_size)
            return correction.str;

        char c = correction.str[pos];

        // See, whether it's enough to perform single nucl extension
        bool extended = false;
        if (is_nucl(c)) {
            KMer last = correction.last << dignucl(c);
            size_t idx = data_.checking_seq_idx(last);
            if (idx != -1ULL) {
                const KMerStat &kmer_data = data_[idx];
                candidates.emplace(pos, correction.str,
                                   correction.penalty - (kmer_data.good() ?
                                                         0.0 :
                                                         (qual[pos] >= 20 ? 1.0 : 2.0)),
                                    last, cpos);
                if (kmer_data.good() && qual[pos] >= 20)
                    extended = true;
            } else {
                candidates.emplace(pos, correction.str,
                                   correction.penalty - (qual[pos] >= 20 ? 2.0 : 3.0),
                                   last, cpos);
            }
        }

        // Ok, it's possible to extend using solely solid k-mer, do not try any other corrections.
        if (extended) {
            FlushCandidates(corrections, candidates, size_thr);
            continue;
        }

        // Do not allow too many corrections
        if (correction.penalty < penalty_thr) {
            FlushCandidates(corrections, candidates, size_thr);
            continue;
        }

        // Do not allow clustered corrections
        if (pos - correction.cpos.front() < pos_thr) {
            // INFO("Cluster " << pos << "," << correction.cpos);
            FlushCandidates(corrections, candidates, size_thr);
            continue;
        }

        // Try corrections
        positions_t cpos = correction.cpos;
        std::copy(cpos.begin() + 1, cpos.end(), cpos.begin());
        cpos.back() = (uint16_t)pos;
        for (char cc = 0; cc < 4; ++cc) {
            char ncc = nucl(cc);
            if (c == ncc)
                continue;

            KMer last = correction.last << cc;
            size_t idx = data_.checking_seq_idx(last);
            if (idx == -1ULL)
                continue;

            const KMerStat &kmer_data = data_[idx];
            if (kmer_data.good()) {
                std::string corrected = correction.str; corrected[pos] = ncc;
                double penalty = correction.penalty - (is_nucl(c) ?
                                                       (qual[pos] >= 20 ? 5.0 : 1.0) :
                                                       0.0);
                candidates.emplace(pos, corrected, penalty, last, cpos);
            }
        }

        FlushCandidates(corrections, candidates, size_thr);
    }

#   pragma omp atomic
    uncorrected_nucleotides_ += read_size - right_pos;

    return seq;
}

bool ReadCorrector::CorrectOneRead(Read & r,
                                   bool, bool, bool) {
    std::string seq = r.getSequenceString();
    const std::string &qual = r.getQualityString();

    size_t read_size = seq.size();

    // Find the longest "solid island"
    size_t lleft_pos = -1ULL, lright_pos = -1ULL, solid_len = 0;

    ValidKMerGenerator<K> gen(seq.data(), qual.data(), read_size);
    size_t left_pos = 0, right_pos = 0;
    while (gen.HasMore()) {
        size_t read_pos = gen.pos() - 1;
        hammer::KMer kmer = gen.kmer();
        size_t idx = data_.checking_seq_idx(kmer);
        if (idx != -1ULL) {
            const KMerStat &kmer_data = data_[idx];
            if (kmer_data.good()) {
                if (read_pos != right_pos - K + 2) {
                    left_pos = read_pos;
                    right_pos = left_pos + K - 1;
                } else
                    right_pos += 1;

                if (right_pos - left_pos + 1 > solid_len) {
                    lleft_pos = left_pos;
                    lright_pos = right_pos;
                    solid_len = right_pos - left_pos + 1;
                }
            }
        }

        // INFO("" << left_pos << ":" << right_pos << ":" << read_pos << ", " << lleft_pos << ":" << lright_pos << "(" << solid_len << "), " << (kmer_data.good() ? "solid" : "non-solid"));

        gen.Next();
    }

#   pragma omp atomic
    total_nucleotides_ += read_size;

    // Now iterate over all the k-mers of a read trying to make all the stuff solid and good.
    if (solid_len && solid_len != read_size) {
        //std::string seq2 = seq;
        //INFO(seq2.insert(lleft_pos, "[").insert(lright_pos + 2, "]"));

        std::string newseq = CorrectReadRight(seq, qual, lright_pos);
        newseq = ReverseComplement(CorrectReadRight(ReverseComplement(newseq), Reverse(qual),
                                                    read_size - 1 - lleft_pos));

        unsigned corrected = 0;
        for (size_t i = 0; i < read_size; ++i)
            corrected += seq[i] != newseq[i];

        if (corrected) {
#           pragma omp atomic
            changed_reads_ += 1;

#           pragma omp atomic
            changed_nucleotides_ += corrected;
            if (correct_stats_) {
                std::string name = r.getName();
                name += " BH:changed:" + std::to_string(corrected);
                r.setName(name.data());
            }
        } else if (correct_stats_) {
            std::string name = r.getName();
            name += " BH:failed";
            r.setName(name.data());
        }

        if (seq.size() != read_size) {
            INFO("Jere");
            return false;
        }

        r.setSequence(newseq.data(), /* preserve_trimming */ true);
        return true;
    } else if (solid_len == read_size && correct_stats_) {
        std::string name = r.getName();
        name += " BH:ok";
        r.setName(name.data());
    }

    return solid_len == read_size;
}
