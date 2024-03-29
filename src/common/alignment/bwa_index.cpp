//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2016-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "bwa_index.hpp"

#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/rle.h"
#include "bwa/rope.h"
#include "bwa/utils.h"

#include <string>
#include <memory>

#define MEM_F_SOFTCLIP  0x200

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= uint8_t((c)<<((~(l)&3)<<1)))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
extern "C" {
int is_bwt(uint8_t *T, bwtint_t n);
int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a);
};

namespace alignment {

BWAIndex::BWAIndex(const debruijn_graph::Graph& g, AlignmentMode mode,
                   RetainAlignments retain)
        : g_(g),
          memopt_(mem_opt_init(), free),
          idx_(nullptr, bwa_idx_destroy),
          mode_(mode),
          retain_(retain) {
    memopt_->flag |= MEM_F_SOFTCLIP;
    switch (mode_) {
        default:
        case AlignmentMode::Illumina:
        case AlignmentMode::Default:
            INFO("Setting BWA alignment mode to 'default'");
            break;
        case AlignmentMode::IntraCtg:
            INFO("Setting BWA alignment mode to 'intractg'");
            memopt_->o_del = 16; memopt_->o_ins = 16;
            memopt_->b = 9;
            memopt_->pen_clip5 = 5; memopt_->pen_clip3 = 5;
            break;
        case AlignmentMode::PacBio:
        case AlignmentMode::Ont2D:
            INFO("Setting BWA alignment mode to 'pacbio'");
            memopt_->o_del = 1; memopt_->e_del = 1;
            memopt_->o_ins = 1; memopt_->e_ins = 1;
            memopt_->b = 1;
            memopt_->split_factor = 10.;
            memopt_->pen_clip5 = 0; memopt_->pen_clip3 = 0;
            memopt_->min_seed_len = 14;
            memopt_->mask_level = 20;
            memopt_->drop_ratio = 20;
            memopt_->min_chain_weight = 40;
            if (retain_ == RetainAlignments::Default)
                retain_ = RetainAlignments::All;
            break;
        case AlignmentMode::HiC:
            INFO("Setting BWA alignment mode to 'hic'");
            // -M5
            // FIXME: Do we realy need -M here?
            memopt_->flag |= MEM_F_NO_MULTI;  // -M
            memopt_->flag |= MEM_F_PRIMARY5;  // -5
            memopt_->pen_clip5 = 0; memopt_->pen_clip3 = 0;
            break;
    };

    if (retain_ == RetainAlignments::Default)
        retain_ = RetainAlignments::OnlyPrimary;

    bwa_fill_scmat(memopt_->a, memopt_->b, memopt_->mat);
    Init();
}

BWAIndex::~BWAIndex() {}

static uint8_t* seqlib_add1(const Sequence &seq,
                            uint8_t *pac, size_t &l_pac, size_t &m_pac) {
    for (size_t i = 0; i < seq.size(); ++i) {
        unsigned c = seq[i];
        // fill buffer
        if (l_pac == m_pac) { // double the pac size
            m_pac <<= 1;
            pac = (uint8_t*)realloc(pac, m_pac/4);
            memset(pac + l_pac/4, 0, (m_pac - l_pac)/4);
        }
        _set_pac(pac, l_pac, c);
        ++l_pac;
    }

    return pac;
}

static uint8_t* seqlib_make_pac(const debruijn_graph::Graph &g,
                                const std::vector<debruijn_graph::EdgeId> &ids,
                                bool for_only) {
    uint8_t *pac = 0;
    size_t m_pac, l_pac;
    ssize_t l;

    m_pac = 0x10000; l_pac = 0;
    pac = (uint8_t*)calloc(m_pac/4, 1);

    // Move through the sequences
    for (auto e : ids) {
        // make the forward only pac
        pac = seqlib_add1(g.EdgeNucls(e), pac, l_pac, m_pac);
    }

    if (!for_only) {
        // add the reverse complemented sequence
        m_pac = (l_pac * 2 + 3) / 4 * 4;
        pac = (uint8_t*)realloc(pac, m_pac/4);
        memset(pac + (l_pac+3)/4, 0, (m_pac - (l_pac+3)/4*4) / 4);
        for (l = l_pac - 1; l >= 0; --l, ++l_pac)
            _set_pac(pac, l_pac, 3 - _get_pac(pac, l));
    }

    return pac;
}

static bwt_t *seqlib_bwt_pac2bwt(const uint8_t *pac, size_t bwt_seq_lenr) {
    bwt_t *bwt;
    ubyte_t *buf;

    // Initialization
    bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
    bwt->seq_len = bwt_seq_lenr;
    bwt->bwt_size = (bwt->seq_len + 15) >> 4;

    // Prepare sequence
    memset(bwt->L2, 0, 5 * 4);
    buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
    for (bwtint_t i = 0; i < bwt->seq_len; ++i) {
        buf[i] = pac[i>>2] >> ((3 - (i&3)) << 1) & 3;
        ++bwt->L2[1+buf[i]];
    }
    for (bwtint_t i = 2; i <= 4; ++i)
        bwt->L2[i] += bwt->L2[i-1];

    // Burrows-Wheeler Transform
    if (bwt_seq_lenr < 50000000) {
        INFO("Using BWA IS algorithm");
        bwt->primary = is_bwt(buf, bwt->seq_len);
    } else {
        INFO("Using BWA RopeBWT algorithm");
        rope_t *r;
        int64_t x, i;
        rpitr_t itr;
        const uint8_t *blk;

        r = rope_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN);
        for (i = bwt->seq_len - 1, x = 0; i >= 0; --i) {
            int c = buf[i] + 1;
            x = rope_insert_run(r, x, c, 1, 0) + 1;
            while (--c >= 0) x += r->c[c];
        }
        bwt->primary = x;
        rope_itr_first(r, &itr);
        x = 0;
        while ((blk = rope_itr_next_block(&itr)) != 0) {
            const uint8_t *q = blk + 2, *end = blk + 2 + *rle_nptr(blk);
            while (q < end) {
                int c = 0;
                int64_t l;
                rle_dec1(q, c, l);
                for (i = 0; i < l; ++i)
                    buf[x++] = ubyte_t(c - 1);
            }
        }
        rope_destroy(r);
    }
    bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
    for (bwtint_t  i = 0; i < bwt->seq_len; ++i)
        bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
    free(buf);
    return bwt;
}

void BWAIndex::Init() {
    idx_.reset((bwaidx_t*)calloc(1, sizeof(bwaidx_t)));
    ids_.clear();

    for (debruijn_graph::EdgeId e : g_.canonical_edges()) {
        ids_.push_back(e);
    }

    // construct the forward-only pac
    uint8_t* fwd_pac = seqlib_make_pac(g_, ids_, true); // true->for_only

    // construct the forward-reverse pac ("packed" 2 bit sequence)
    uint8_t* pac = seqlib_make_pac(g_, ids_, false); // don't write, because only used to make BWT

    size_t tlen = 0;
    for (auto e : ids_)
        tlen += g_.EdgeNucls(e).size();

    // make the bwt
    bwt_t *bwt;
    bwt = seqlib_bwt_pac2bwt(pac, tlen*2); // *2 for fwd and rev
    bwt_bwtupdate_core(bwt);
    free(pac); // done with fwd-rev pac

    // construct sa from bwt and occ. adds it to bwt struct
    bwt_cal_sa(bwt, 32);
    bwt_gen_cnt_table(bwt);

    // make the bns
    bntseq_t * bns = (bntseq_t*) calloc(1, sizeof(bntseq_t));
    bns->l_pac = tlen;
    bns->n_seqs = int(ids_.size());
    bns->seed = 11;
    bns->n_holes = 0;

    // make the anns
    bns->anns = (bntann1_t*)calloc(ids_.size(), sizeof(bntann1_t));
    size_t offset = 0, k = 0;
    for (auto e: ids_) {
        bntann1_t *ann = &bns->anns[k++];
        int len = int(g_.EdgeNucls(e).size());

        ann->offset = offset;
        ann->name = ann->anno = nullptr;
        ann->len = len;
        ann->n_ambs = 0; // number of "holes"
        ann->gi = 0; // gi?
        ann->is_alt = 0;

        offset += len;
    }

    // ambs is "holes", like N bases
    bns->ambs = nullptr;

    // Make the in-memory idx struct
    idx_->bwt = bwt;
    idx_->bns = bns;
    idx_->pac = fwd_pac;
}

#if 0
        fprintf(stderr, "%zu: [%lld, %lld)\t[%d, %d) %c %d %s %ld %d\n",
                i,
                pos, pos + a.re - a.rb, a.qb, a.qe,
                "+-"[is_rev], a.rid,
                idx_->bns->anns[a.rid].name, g_.int_id(ids_[a.rid]), a.secondary);
#endif

#if 0
        mem_aln_t aln = mem_reg2aln(memopt_.get(), idx_->bns, idx_->pac, seq.length(), seq.c_str(), &a);

        // print alignment
        fprintf(stderr, "\t%c\t%s\t%ld %ld %ld\t%d\t", "+-"[aln.is_rev], idx_->bns->anns[aln.rid].name, aln.rid, g_.int_id(ids_[aln.rid]), (long)aln.pos, aln.mapq);
        for (int k = 0; k < aln.n_cigar; ++k) // print CIGAR
            fprintf(stderr, "%d%c", aln.cigar[k]>>4, "MIDSH"[aln.cigar[k]&0xf]);
        fprintf(stderr, "\t%d\n", aln.NM); // print edit distance
        free(aln.cigar);
#endif

inline std::ostream& operator<<(std::ostream& os, const mem_alnreg_s& a) {
    os << a.qb << " - " << a.qe << " ---> " << a.rb << " - " << a.re << "query->ref\n";
    os << a.seedcov << " - seedcov; " << a.score << " - score";
    return os;
}

static bool MostlyInVertex(size_t rb, size_t re, size_t edge_len, size_t k) {
//  k-rb > re - k
    if (rb < k && 2 * k  > re + rb)
        return true;
//  re - edge_len > edge_len - rb
    if (re > edge_len && re + rb > 2 * edge_len)
        return true;
    return false;
}


omnigraph::MappingPath<debruijn_graph::EdgeId> BWAIndex::GetMappingPath(const mem_alnreg_v &ar, const Sequence &seq,
                                                                        bool only_simple) const {
    omnigraph::MappingPath<debruijn_graph::EdgeId> res;

    size_t seq_len = seq.size();
    bool is_short = seq_len <= g_.k();

    for (size_t i = 0; i < ar.n; ++i) {
        const mem_alnreg_t &a = ar.a[i];

        if (a.score < memopt_->T) // skip low-scored alignments
            continue;

        if (retain_ != RetainAlignments::All && a.secondary >= 0) continue; // skip secondary alignments

        if (is_short) {
            // Skipping alignments shorter than half of read length
            if (size_t(a.qe - a.qb) * 2 <= seq_len) continue;
            if (size_t(a.re - a.rb) * 2 <= seq_len) continue;
        } else {
            if (size_t(a.qe - a.qb) <= g_.k()) continue; // skip short alignments
            if (size_t(a.re - a.rb) <= g_.k()) continue;
        }
        int is_rev = 0;
        bntann1_t *ans = idx_->bns->anns + a.rid;
        size_t pos = bns_depos(idx_->bns, a.rb < idx_->bns->l_pac? a.rb : a.re - 1, &is_rev) - ans->offset;
        size_t initial_range_end;
        size_t mapping_range_end;

        // Reduce the range to kmer-based
        if (is_short) {
            initial_range_end = a.qb + 1;
            mapping_range_end = pos + 1;
            if (mapping_range_end > g_.length(ids_[a.rid]))
                continue;
        } else {
            initial_range_end = a.qe - g_.k();
            mapping_range_end = pos + a.re - a.rb - g_.k();
        }
        // Important for alignments shorter than K
        if (MostlyInVertex(pos, pos + a.re - a.rb, g_.length(ids_[a.rid]), g_.k()))
            continue;

        if (retain_ == RetainAlignments::QualityPrimary &&
            mem_approx_mapq_se(memopt_.get(), &a) < 20) // Low-quality alignment: multi-mapped, etc.
            continue;

        DEBUG_EXPR({
            std::string sseq = seq.str();
            mem_aln_t aln = mem_reg2aln(memopt_.get(), idx_->bns, idx_->pac, int(sseq.length()), sseq.c_str(), &a);

            fprintf(stderr, "%zu: [%zu, %lld)\t[%d, %d) %c %ld %s %d %d\n",
                    i,
                    pos, pos + a.re - a.rb, a.qb, a.qe,
                    "+-"[is_rev], g_.int_id(ids_[a.rid]), a.secondary >= 0 ? "secondary" : "primary", aln.mapq, aln.NM);

            // print alignment
            for (int k = 0; k < aln.n_cigar; ++k) // print CIGAR
                fprintf(stderr, "%d%c", aln.cigar[k]>>4, "MIDSH"[aln.cigar[k]&0xf]);
            free(aln.cigar);
        });

        DEBUG(a);
        // FIXME: what about other scoring systems?
        double qual = double(a.score)/double(a.qe - a.qb);
        DEBUG("Edge: "<< ids_[a.rid] << " quality from score: " << qual);

        if (!is_rev) {
            res.push_back(ids_[a.rid],
                          { { (size_t)a.qb, initial_range_end },
                            { pos, mapping_range_end}, qual});
        } else {
            res.push_back(g_.conjugate(ids_[a.rid]),
                          { { (size_t)a.qb, initial_range_end }, //.Invert(read_length),
                            Range(pos,  mapping_range_end).Invert(g_.length(ids_[a.rid])) , qual});

        }
        if (only_simple && res.size() > 1)
            return omnigraph::MappingPath<debruijn_graph::EdgeId>();
    }
    return res;
}


omnigraph::MappingPath<debruijn_graph::EdgeId> BWAIndex::AlignSequence(const Sequence &sequence,
                                                                       bool only_simple) const {
    omnigraph::MappingPath<debruijn_graph::EdgeId> res;
    VERIFY(idx_);

    char *seq = new char[sequence.size()];
    for (size_t i = 0; i < sequence.size(); ++i)
        seq[i] = sequence[i];
    // WARNING: This function modifies seq (for purpose)!
    mem_alnreg_v ar = mem_align1_bin(memopt_.get(), idx_->bwt, idx_->bns, idx_->pac,
                                     int(sequence.size()), seq);
    res = GetMappingPath(ar, sequence, only_simple);

    free(ar.a);
    delete[] seq;

    return res;
}

}
