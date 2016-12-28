//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "bwa_index.hpp"

#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/utils.h"
#include "kseq/kseq.h"

#include <string>
#include <memory>

// all of the bwa and kseq stuff is in unaligned sequence
// best way I had to keep from clashes with klib macros

#define MEM_F_SOFTCLIP  0x200

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
extern "C" {
int is_bwt(uint8_t *T, int n);
};

namespace alignment {

BWAIndex::BWAIndex(const debruijn_graph::Graph& g)
        : g_(g),
          memopt_(mem_opt_init(), free),
          idx_(nullptr, bwa_idx_destroy) {
    memopt_->flag |= MEM_F_SOFTCLIP;
    Init();
}

BWAIndex::~BWAIndex() {}

// modified from bwa (heng li)
static uint8_t* seqlib_add1(const kstring_t *seq, const kstring_t *name,
                            bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q) {
    bntann1_t *p;
    int lasts;
    if (bns->n_seqs == *m_seqs) {
        *m_seqs <<= 1;
        bns->anns = (bntann1_t*)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
    }
    p = bns->anns + bns->n_seqs;
    p->name = strdup((char*)name->s);
    p->anno = strdup("(null");
    p->gi = 0; p->len = seq->l;
    p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
    p->n_ambs = 0;
    for (size_t i = lasts = 0; i < seq->l; ++i) {
        int c = nst_nt4_table[(int)seq->s[i]];
        if (c >= 4) { // N
            if (lasts == seq->s[i]) { // contiguous N
                ++(*q)->len;
            } else {
                if (bns->n_holes == *m_holes) {
                    (*m_holes) <<= 1;
                    bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
                }
                *q = bns->ambs + bns->n_holes;
                (*q)->len = 1;
                (*q)->offset = p->offset + i;
                (*q)->amb = seq->s[i];
                ++p->n_ambs;
                ++bns->n_holes;
            }
        }
        lasts = seq->s[i];
        { // fill buffer
            if (c >= 4) c = lrand48()&3;
            if (bns->l_pac == *m_pac) { // double the pac size
                *m_pac <<= 1;
                pac = (uint8_t*)realloc(pac, *m_pac/4);
                memset(pac + bns->l_pac/4, 0, (*m_pac - bns->l_pac)/4);
            }
            _set_pac(pac, bns->l_pac, c);
            ++bns->l_pac;
        }
    }
    ++bns->n_seqs;

    return pac;
}

static uint8_t* seqlib_make_pac(const debruijn_graph::Graph &g,
                                const std::vector<debruijn_graph::EdgeId> &ids,
                                bool for_only) {
    bntseq_t * bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
    uint8_t *pac = 0;
    int32_t m_seqs, m_holes;
    int64_t m_pac, l;
    bntamb1_t *q;

    bns->seed = 11; // fixed seed for random generator
    m_seqs = m_holes = 8; m_pac = 0x10000;
    bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
    bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
    pac = (uint8_t*) calloc(m_pac/4, 1);
    q = bns->ambs;

    // move through the sequences
    // FIXME: not kstring is required
    for (auto e : ids) {
        std::string ref = std::to_string(g.int_id(e));
        std::string seq = g.EdgeNucls(e).str();

        // make the ref name kstring
        kstring_t * name = (kstring_t*)malloc(1 * sizeof(kstring_t));
        name->l = ref.length() + 1;
        name->m = ref.length() + 3;
        name->s = (char*)calloc(name->m, sizeof(char));
        memcpy(name->s, ref.c_str(), ref.length()+1);

        // make the sequence kstring
        kstring_t * t = (kstring_t*)malloc(sizeof(kstring_t));
        t->l = seq.length();
        t->m = seq.length() + 2;
        //t->s = (char*)calloc(v[k].Seq.length(), sizeof(char));
        t->s = (char*)malloc(t->m);
        memcpy(t->s, seq.c_str(), seq.length());

        // make the forward only pac
        pac = seqlib_add1(t, name, bns, pac, &m_pac, &m_seqs, &m_holes, &q);

        // clear it out
        free(name->s);
        free(name);
        free(t->s);
        free(t);
    }

    if (!for_only) {
        // add the reverse complemented sequence
        m_pac = (bns->l_pac * 2 + 3) / 4 * 4;
        pac = (uint8_t*)realloc(pac, m_pac/4);
        memset(pac + (bns->l_pac+3)/4, 0, (m_pac - (bns->l_pac+3)/4*4) / 4);
        for (l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac)
            _set_pac(pac, bns->l_pac, 3-_get_pac(pac, l));
    }

    bns_destroy(bns);

    return pac;
}

static bwt_t *seqlib_bwt_pac2bwt(const uint8_t *pac, size_t bwt_seq_lenr) {
    bwt_t *bwt;
    ubyte_t *buf;
    int i;

    // initialization
    bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
    bwt->seq_len = bwt_seq_lenr; //bwa_seq_len(fn_pac); //dummy
    bwt->bwt_size = (bwt->seq_len + 15) >> 4;

    // prepare sequence
    //pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
    //buf2 = (ubyte_t*)calloc(pac_size, 1);
    //err_fread_noeof(buf2, 1, pac_size, fp);
    //err_fclose(fp);
    memset(bwt->L2, 0, 5 * 4);
    buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
    for (i = 0; i < (int)bwt->seq_len; ++i) {
        buf[i] = pac[i>>2] >> ((3 - (i&3)) << 1) & 3;
        ++bwt->L2[1+buf[i]];
    }
    for (i = 2; i <= 4; ++i)
        bwt->L2[i] += bwt->L2[i-1];
    //free(buf2);

    // Burrows-Wheeler Transform
    bwt->primary = is_bwt(buf, bwt->seq_len);
    bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
    for (i = 0; i < (int)bwt->seq_len; ++i)
        bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
    free(buf);
    return bwt;
}

static bntann1_t* seqlib_add_to_anns(const std::string& name, const std::string& seq, bntann1_t* ann, size_t offset) {
    ann->offset = offset;
    ann->name = (char*)malloc(name.length()+1); // +1 for \0
    strncpy(ann->name, name.c_str(), name.length()+1);
    ann->anno = (char*)malloc(7);
    strcpy(ann->anno, "(null)\0");
    ann->len = seq.length();
    ann->n_ambs = 0; // number of "holes"
    ann->gi = 0; // gi?
    ann->is_alt = 0;

    return ann;
}

void BWAIndex::Init() {
    idx_.reset((bwaidx_t*)calloc(1, sizeof(bwaidx_t)));
    ids_.clear();

    for (auto it = g_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        ids_.push_back(*it);
    }

    // construct the forward-only pac
    uint8_t* fwd_pac = seqlib_make_pac(g_, ids_, true); //true->for_only

    // construct the forward-reverse pac ("packed" 2 bit sequence)
    uint8_t* pac = seqlib_make_pac(g_, ids_, false); // don't write, becasue only used to make BWT

    size_t tlen = 0;
    for (auto e : ids_)
        tlen += g_.EdgeNucls(e).size();

#ifdef DEBUG_BWATOOLS
    std::cerr << "ref seq length: " << tlen << std::endl;
#endif

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
    bns->n_seqs = ids_.size();
    bns->seed = 11;
    bns->n_holes = 0;

    // make the anns
    // FIXME: Do we really need this?
    bns->anns = (bntann1_t*)calloc(ids_.size(), sizeof(bntann1_t));
    size_t offset = 0, k = 0;
    for (auto e: ids_) {
        std::string name = std::to_string(g_.int_id(e));
        std::string seq = g_.EdgeNucls(e).str();
        seqlib_add_to_anns(name, seq, &bns->anns[k++], offset);
        offset += seq.length();
    }

    // ambs is "holes", like N bases
    bns->ambs = 0; //(bntamb1_t*)calloc(1, sizeof(bntamb1_t));

    // make the in-memory idx struct
    idx_->bwt = bwt;
    idx_->bns = bns;
    idx_->pac = fwd_pac;
}

omnigraph::MappingPath<debruijn_graph::EdgeId> BWAIndex::AlignSequence(const Sequence &sequence) const {
    omnigraph::MappingPath<debruijn_graph::EdgeId> res;

    if (!idx_) return res;

    std::string seq = sequence.str();
    mem_alnreg_v ar = mem_align1(memopt_.get(), idx_->bwt, idx_->bns, idx_->pac,
                                 seq.length(), seq.data());
    for (size_t i = 0; i < ar.n; ++i) {
        const mem_alnreg_t &a = ar.a[i];
        if (a.secondary >= 0) continue; // skip secondary alignments
//        if (a.qe - a.qb < g_.k()) continue; // skip short alignments
//        if (a.re - a.rb < g_.k()) continue;
        int is_rev = 0;
        size_t pos = bns_depos(idx_->bns, a.rb < idx_->bns->l_pac? a.rb : a.re - 1, &is_rev) - idx_->bns->anns[a.rid].offset;
/*        fprintf(stderr, "%zu: [%lld, %lld]\t[%d, %d] %c %d %s %ld %zu\n",
                i,
                a.rb, a.re, a.qb, a.qe,
                "+-"[is_rev], a.rid,
                idx_->bns->anns[a.rid].name, g_.int_id(ids_[a.rid]), pos);
*/
        size_t initial_range_end = a.qe;
        size_t mapping_range_end = pos + a.re - a.rb;
        size_t read_length = seq.length() ;
        //we had to reduce the range to kmer-based
        if (pos + (a.re - a.rb) >= g_.length(ids_[a.rid]) ){
            if (a.qe > g_.k() + a.qb)
                initial_range_end -= g_.k();
            else continue;
            if (a.re > g_.k() + a.rb)
                mapping_range_end -= g_.k();
            else continue;
            if (read_length >= g_.k())
                read_length -= g_.k();
            else continue;
        }
        // FIXME: Check this!
        if (!is_rev) {
            res.push_back(ids_[a.rid],
                          { { (size_t)a.qb, initial_range_end },
                            { pos, mapping_range_end}});
        } else {
//          fprintf (stderr,"%d %d %d\n", a.qb, a.qe  - g_.k(), seq.length() - g_.k());

//            fprintf (stderr,"%d %d %d\n", pos, pos + a.re - a.rb , g_.length(ids_[a.rid]) );

            res.push_back(g_.conjugate(ids_[a.rid]),
                          { omnigraph::Range(a.qb, initial_range_end).Invert(read_length),
                            omnigraph::Range(pos, mapping_range_end ).Invert(g_.length(ids_[a.rid])) });

        }

#if 0
        mem_aln_t aln = mem_reg2aln(memopt_.get(), idx_->bns, idx_->pac, seq.length(), seq.c_str(), &a);

        // print alignment
        printf("\t%c\t%s\t%ld %ld %ld\t%d\t", "+-"[aln.is_rev], idx_->bns->anns[aln.rid].name, aln.rid, g_.int_id(ids_[aln.rid]), (long)aln.pos, aln.mapq);
        for (int k = 0; k < aln.n_cigar; ++k) // print CIGAR
            printf("%d%c", aln.cigar[k]>>4, "MIDSH"[aln.cigar[k]&0xf]);
        printf("\t%d\n", aln.NM); // print edit distance
        free(aln.cigar);
#endif

    }
    free(ar.a);

    return res;
}

}
