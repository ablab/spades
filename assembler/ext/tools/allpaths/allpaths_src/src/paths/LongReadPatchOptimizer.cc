///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS


#include <omp.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/LongReadTools.h"



// This class describes 2 bases at a specific point in an alignment
// At least one of the bases must be present, a or b

class BasesAligned
{
  uint8_t _base_a : 2;  // 2 bits for base in read A
  uint8_t _base_b : 2;  // 2 bits for base in read B
  uint8_t _gap_a  : 1;  // 1 bit  for gap in read A (non-existent base in read A) 
  uint8_t _gap_b  : 1;  // 1 bit  for gap in read B (non-existent base in read B) 

public:
  BasesAligned(const unsigned ba, const unsigned bb,
               const bool ga, const bool gb)
    : _base_a(ba),
      _base_b(bb),
      _gap_a(ga),
      _gap_b(gb)
  {}
  
  unsigned base_a () const { return _base_a; }
  unsigned base_b () const { return _base_b; }
  bool     gap_a  () const { return _gap_a; }
  bool     gap_b  () const { return _gap_b; }

  bool     gapped       () const { return gap_a() || gap_b(); }
  bool     supported    () const { return !gapped(); }
  bool     matched      () const { return supported() && base_a() == base_b(); }
  bool     substitution () const { return supported() && base_a() != base_b(); }
};





// A class that stores the correspondence between the global alignment index 'j'
// and the indices in read A 'ia' and read b 'ib'

class AlignmentIndexes
{
  vec<int> _j_of_ia;   // j as a function of ia
  vec<int> _j_of_ib;   // j as a function of ib

  // below is optimizable but not now. Only one is enough.

  vec<int> _ia_le_j;   // ia as a function of j (if gap choose lower ia) 
  vec<int> _ia_ge_j;   // ia as a function of j (if gap choose higher ia) 
  vec<int> _ib_le_j;   // ib as a function of j (if gap choose lower ib) 
  vec<int> _ib_ge_j;   // ib as a function of j (if gap choose higher ib) 

  vec<BasesAligned> _bases; // stores all the bases and gaps as a function of j
    
public:  
  AlignmentIndexes(const BaseVec & bv_a, 
                   const BaseVec & bv_b,
                   const alignment & al) 
  {
    int ia;
    int ib;
    int errors;
    avector<int> gaps;
    avector<int> lengths;
    al.Unpack(ia, ib, errors, gaps, lengths);
    //cout << "ia = " << ia << endl;
    //cout << "ib = " << ib << endl;
    //cout << "ngaps = " << gaps.length << endl;
    //cout << "nlens = " << lengths.length << endl;
    
    for (int i = 0; i != ia; i++) _j_of_ia.push_back(i - ia); 
    for (int i = 0; i != ib; i++) _j_of_ib.push_back(i - ib);
    
    int j = 0;

    const int n_blocks = gaps.length + lengths.length;

    for (int i_block = 1; i_block != n_blocks; i_block++) { // first gap always = 0;
      if (i_block % 2 == 1) { // a length
        const int len = lengths(i_block / 2);
        //cout << " l=" << len;

        for (int k = 0; k != len; k++) {
          _bases.push_back(BasesAligned(bv_a[ia], bv_b[ib], false, false));
          _j_of_ia.push_back(j);
          _j_of_ib.push_back(j);
          _ia_le_j.push_back(ia);
          _ia_ge_j.push_back(ia);
          _ib_le_j.push_back(ib);
          _ib_ge_j.push_back(ib);
          ia++;
          ib++;
          j++;
        }
      }
      else { // a gap
        const int gap = gaps(i_block / 2);
        //cout << " g=" << gap;

        if (gap < 0) { // gap in read b
          for (int k = 0; k != -gap; k++) {
            _bases.push_back(BasesAligned(bv_a[ia], 0, false, true));
            _j_of_ia.push_back(j);
            _ia_le_j.push_back(ia);
            _ia_ge_j.push_back(ia);
            _ib_le_j.push_back(ib);
            _ib_ge_j.push_back(ib + 1);
            ia++;
            j++;
          }
        }
        else { // gap in read a
          for (int k = 0; k != gap; k++) {
            _bases.push_back(BasesAligned(0, bv_b[ib], true, false));
            _j_of_ib.push_back(j);
            _ia_le_j.push_back(ia);
            _ia_ge_j.push_back(ia + 1);
            _ib_le_j.push_back(ib);
            _ib_ge_j.push_back(ib);
            ib++;
            j++;
          }
        }
      }
    }
    //cout << endl;
  }
  
  int index_B_le_index_A(const int ia) const { return _ib_le_j[_j_of_ia[ia]]; }
  int index_B_ge_index_A(const int ia) const { return _ib_ge_j[_j_of_ia[ia]]; }
  int index_A_le_index_B(const int ib) const { return _ia_le_j[_j_of_ib[ib]]; }
  int index_A_ge_index_B(const int ib) const { return _ia_ge_j[_j_of_ib[ib]]; }


  unsigned base_A_at_A(const int ia) const { return _bases[_j_of_ia[ia]].base_a(); }
  unsigned base_B_at_B(const int ib) const { return _bases[_j_of_ib[ib]].base_b(); }

  bool supported_B_at_A(const int ia) const { return _bases[_j_of_ia[ia]].supported(); }
  bool supported_A_at_B(const int ib) const { return _bases[_j_of_ib[ib]].supported(); }
 
  bool matched_B_at_A(const int ia) const { return _bases[_j_of_ia[ia]].matched(); }
  bool matched_A_at_B(const int ib) const { return _bases[_j_of_ib[ib]].matched(); }
  
  bool deletion_B_at_A(const int ia) const { return _bases[_j_of_ia[ia]].gap_b(); }
  bool deletion_A_at_B(const int ib) const { return _bases[_j_of_ib[ib]].gap_a(); }

  bool substitution_B_at_A(const int ia) const { return _bases[_j_of_ia[ia]].substitution(); }
  bool substitution_A_at_B(const int ib) const { return _bases[_j_of_ib[ib]].substitution(); }

  unsigned base_B_at_A(const int ia) const { return _bases[_j_of_ia[ia]].base_b(); }
  unsigned base_A_at_B(const int ib) const { return _bases[_j_of_ib[ib]].base_a(); }


  bool deletion_B_before_A(const int ia) const { return _bases[_j_of_ia[ia] - 1].gap_b(); }
  bool deletion_A_before_B(const int ib) const { return _bases[_j_of_ib[ib] - 1].gap_a(); }

  bool insertion_B_before_A(const int ia) const { return _bases[_j_of_ia[ia] - 1].gap_a(); }
  bool insertion_A_before_B(const int ib) const { return _bases[_j_of_ib[ib] - 1].gap_b(); }

  unsigned base_B_before_A(const int ia) const { return _bases[_j_of_ia[ia] - 1].base_b(); }
  unsigned base_A_before_B(const int ib) const { return _bases[_j_of_ib[ib] - 1].base_a(); }


  bool deletion_B_after_A(const int ia) const { return _bases[_j_of_ia[ia] + 1].gap_b(); }
  bool deletion_A_after_B(const int ib) const { return _bases[_j_of_ib[ib] + 1].gap_a(); }

  bool insertion_B_after_A(const int ia) const { return _bases[_j_of_ia[ia] + 1].gap_a(); }
  bool insertion_A_after_B(const int ib) const { return _bases[_j_of_ib[ib] + 1].gap_b(); }

  unsigned base_B_after_A(const int ia) const { return _bases[_j_of_ia[ia] + 1].base_b(); }
  unsigned base_A_after_B(const int ib) const { return _bases[_j_of_ib[ib] + 1].base_a(); }
          

  void print_simple(ostream & out = cout) const 
  {
    vec<String> b(4);
    b[0] = "a";
    b[1] = "c";
    b[2] = "g";
    b[3] = "t";

    const int n = _bases.size();
    for (int j = 0; j != n; j++) out << (j % 10);
    out << endl;
    for (int j = 0; j != n; j++) 
      out << (_bases[j].gap_a() ? " " : 
              (_bases[j].base_a() == _bases[j].base_b() || _bases[j].gap_b() ?
               "-" : b[_bases[j].base_a()]));
    out << endl;
    for (int j = 0; j != n; j++) 
      out << (_bases[j].gap_b() ? " " : 
              (_bases[j].base_a() == _bases[j].base_b() || _bases[j].gap_a() ?
               "-" : b[_bases[j].base_b()]));
    out << endl;
  }

  void print(ostream & out = cout) const 
  {
    vec<String> b(4);
    b[0] = "a";
    b[1] = "c";
    b[2] = "g";
    b[3] = "t";

    const int n = _bases.size();
    for (int j = 0; j != n; j++) out << (j % 10);
    out << endl;
    for (int j = 0; j != n; j++) 
      out << (_bases[j].gap_a() ? " " : b[_bases[j].base_a()]);
    out << endl;
    for (int j = 0; j != n; j++) 
      out << (_bases[j].gap_b() ? " " : b[_bases[j].base_b()]);
    out << endl;
  }
 

};






// Very simple classes to describe indels and substitutions and the distance difference
//
// d_dist = dist_mod - dist_orig
//
// d_dist < 0  corresponds to an improvement

class Deletion
{
public:
  int ib;
  int d_dist;
  Deletion(const int ib, const int d_dist) 
    : ib(ib), d_dist(d_dist) {}
};


class Insertion
{
public:
  int ib;
  int d_dist;
  unsigned base;
  Insertion(const int ib, const int d_dist, const unsigned base) 
    : ib(ib), d_dist(d_dist), base(base) {}
};

typedef Insertion Substitution;












// aligns bv1 to all of bvs2 with Smith-Waterman
// and returns a vector of AligmentIndexes

int dist_base_vecs_SW(const BaseVec & bv1, 
                      const vec<BaseVec> & bvs2, 
                      vec<AlignmentIndexes> * p_alis = 0,
                      const unsigned verbosity = 0)
{
  if (p_alis) p_alis->clear();
  const size_t n2 = bvs2.size();
  int dist = 0;
  //cout << endl;
  for (size_t i2 = 0; i2 < n2; i2++) {
    const BaseVec & bv2 = bvs2[i2];
   
    int best_loc;
    int mismatch_penalty = 1;
    int gap_penalty = 1;
    alignment al;
    if (bv2.size() == 0) {
      dist = bv1.size() * gap_penalty;
    }
    else {
      if (bv1.size() <= bv2.size()) {
        dist += SmithWatFree(bv1, bv2, best_loc, al, true, true, mismatch_penalty, gap_penalty);
      }
      else {
        dist += SmithWatFree(bv2, bv1, best_loc, al, true, true, mismatch_penalty, gap_penalty);
        al.Flip();
      }
    }
    
    if (verbosity >= 2) {
      #pragma omp critical 
      {
        cout << "-------------------- i2= " << i2 << " dist= " << dist << endl;
        PrintVisualAlignment(true, cout, bv1, bv2, al);
      }
    }
    if (p_alis) p_alis->push_back(AlignmentIndexes(bv1, bv2, al));
  }

  //cout << "SW: " << bv1.size() << " x " << bv2.size() << " score= " << score << endl;
  
  //AlignmentIndexes an(bv1, bv2, *p_al);
  //an.print();
  //exit(0);
  return dist;
}













void deletions_find(const BaseVec & bv, 
                    const int & ib0, 
                    const int & ib1,
                    const vec<BaseVec> & bvs,
                    const vec<AlignmentIndexes> & alis,
                    const int nb_rad,
                    vec<Deletion> * p_dels,
                    double * time = 0,
                    const unsigned PATCH_VERBOSITY = 0)
{
  p_dels->clear();

  const bool verbose = (PATCH_VERBOSITY >= 4);

  const size_t nbv = bvs.size();

  if (verbose) {
    #pragma omp critical
    for (size_t ibv = 0; ibv != nbv; ibv++) {
      cout << "ibv= " << ibv << endl;
      alis[ibv].print_simple(cout);
    }
  }

  const int nb = bv.size();
  for (int ib = ib0; ib <= ib1; ib++) {

    const int ib0_sub = ib - nb_rad;
    const int ib1_sub = ib + nb_rad;
    const int nb_sub = ib1_sub - ib0_sub + 1;

    // ---- build bv_sub_orig
    BaseVec bv_sub_orig(bv, ib0_sub, nb_sub);

    // ---- build bvs_sub
    vec<BaseVec> bvs_sub(nbv);
    size_t size_sum = 0;
    for (size_t ibv = 0; ibv != nbv; ibv++) {

      const int jb0 = alis[ibv].index_B_ge_index_A(ib0_sub);
      const int jb1 = alis[ibv].index_B_le_index_A(ib1_sub);
      if (jb0 >= bvs[ibv].isize() || 
          jb1 > bvs[ibv].isize()) {
        cout << "OMFG!!! asking for sub bv out of range of bv" << endl
             << " ib0= " << ib0 
             << " ib= " << ib 
             << " ib1= " << ib1 
             << " nb= " << nb << endl
             << "nbv= " << nbv
             << " ibv= " << ibv 
             << " bvs[ibv].size= " << bvs[ibv].size() << endl
             << "jb0= " << jb0 << endl
             << "jb1= " << jb1 << endl
             << "bvs[ibv]= " << bvs[ibv] << endl;
      }
      if (jb1 - jb0 + 1 > 0) {
        bvs_sub[ibv].SetToSubOf(bvs[ibv], jb0, jb1 - jb0 + 1);
      }
      else if (jb1 - jb0 + 1 < 0) {
        cout << "OOOOPPPPSSSS!!!!! asking for negative sub bv" << endl
             << " ib0= " << ib0 
             << " ib= " << ib 
             << " ib1= " << ib1 
             << " nb= " << nb << endl
             << "nbv= " << nbv 
             << " ibv= " << ibv
             << " bvs[ibv].size= " << bvs[ibv].size() << endl
             << "jb0= " << jb0 << endl
             << "jb1= " << jb1 << endl
             << "jb1 - jb0 + 1 = " << jb1 - jb0 + 1 << endl
             << "bv= " << bvs[ibv] << endl;
      }
      size_sum += jb1 - jb0 + 1;
    }

    if (time) *time -= WallClockTime();
    // ---- align bv_sub_orig to bvs_sub
    const int dist_orig = dist_base_vecs_SW(bv_sub_orig, bvs_sub);
    if (time) *time += WallClockTime();
      
    // ---- build bv_sub_mod
    BaseVec bv_sub_mod(nb_sub - 1);
    for (int jb = 0; jb < nb_rad; jb++) {
      bv_sub_mod.set(         jb, bv_sub_orig[             jb]);
      bv_sub_mod.set(nb_rad + jb, bv_sub_orig[1 + nb_rad + jb]);
    }

    // ---- align bv_sub_mod to bvs_sub
    if (time) *time -= WallClockTime();
    const int d_dist = dist_base_vecs_SW(bv_sub_mod, bvs_sub) - dist_orig;
    if (time) *time += WallClockTime();

    // ---- decide whether to keep modification
    if (d_dist < 0) {
      p_dels->push_back(Deletion(ib, d_dist));
      ib += nb_rad; // so that no changes overlap
    }
  }
 
}




void insertions_find(const BaseVec & bv, 
                     const int & ib0, 
                     const int & ib1,
                     const vec<BaseVec> & bvs,
                     const vec<AlignmentIndexes> & alis,
                     const int nb_rad,
                     vec<Insertion> * p_inss,
                     double * time = 0)
{
  p_inss->clear();

  const size_t nbv = bvs.size();
  const int nb = bv.size();
  for (int ib = ib0; ib <= ib1; ib++) {

    const int ib0_sub = ib - nb_rad;
    const int ib1_sub = ib + nb_rad - 1;
    const int nb_sub = ib1_sub - ib0_sub + 1;

    // ---- build bv_sub_orig
    BaseVec bv_sub_orig(bv, ib0_sub, nb_sub);

    // ---- build bvs_sub
    vec<BaseVec> bvs_sub(nbv);
    for (size_t ibv = 0; ibv != nbv; ibv++) {
      const int jb0 = alis[ibv].index_B_ge_index_A(ib0_sub);
      const int jb1 = alis[ibv].index_B_le_index_A(ib1_sub);
      bvs_sub[ibv].SetToSubOf(bvs[ibv], jb0, jb1 - jb0 + 1);
    }

    // ---- align bv_sub_orig to bvs_sub
    if (time) *time -= WallClockTime();
    const int dist_orig = dist_base_vecs_SW(bv_sub_orig, bvs_sub);
    if (time) *time += WallClockTime();
    
    // ---- build bv_sub_mod and aligning it for each of 4 bases
    BaseVec bv_sub_mod(nb_sub + 1);
    for (int jb = 0; jb < nb_rad; jb++) {
      bv_sub_mod.set(             jb, bv_sub_orig[         jb]);
      bv_sub_mod.set(1 + nb_rad + jb, bv_sub_orig[nb_rad + jb]);
    }
    
    Insertion ins(ib, 0, 0);
    if (time) *time -= WallClockTime();
    for (unsigned b_new = 0; b_new != 4; b_new++) {
      bv_sub_mod.set(nb_rad, b_new);
      const int d_dist = dist_base_vecs_SW(bv_sub_mod, bvs_sub) - dist_orig;
      if (d_dist < ins.d_dist) {
        ins.d_dist = d_dist;
        ins.base = b_new;
      } 
    }     
    if (time) *time += WallClockTime();

    // ---- decide whether to keep modification
    if (ins.d_dist < 0) {
      p_inss->push_back(ins);
      ib += nb_rad; // so that no changes overlap
    }
  }
}












void substitutions_find(const BaseVec & bv, 
                        const int & ib0, 
                        const int & ib1,
                        const vec<BaseVec> & bvs,
                        const vec<AlignmentIndexes> & alis,
                        const int nb_rad,
                        vec<Substitution> * p_subs,
                        double * time = 0)
{
  p_subs->clear();
  const size_t nbv = bvs.size();
  const int nb = bv.size();
  for (int ib = ib0; ib <= ib1; ib++) {

    const int ib0_sub = ib - nb_rad;
    const int ib1_sub = ib + nb_rad;
    const int nb_sub = ib1_sub - ib0_sub + 1;

    // ---- build bv_sub_orig
    BaseVec bv_sub_orig(bv, ib0_sub, nb_sub);

    // ---- build bvs_sub
    vec<BaseVec> bvs_sub(nbv);
    for (size_t ibv = 0; ibv != nbv; ibv++) {
      int jb0 = alis[ibv].index_B_ge_index_A(ib0_sub);
      int jb1 = alis[ibv].index_B_le_index_A(ib1_sub);
      bvs_sub[ibv].SetToSubOf(bvs[ibv], jb0, jb1 - jb0 + 1);
    }

    // ---- align bv_sub_orig to bvs_sub
    if (time) *time -= WallClockTime();
    int dist_orig = dist_base_vecs_SW(bv_sub_orig, bvs_sub);
    
    // ---- build bv_sub_mod and aligning it for each of 3 bases
    BaseVec bv_sub_mod = bv_sub_orig;

    Substitution sub(ib, 0, 0);
    for (unsigned db = 0; db != 3; db++) {
      const unsigned b_new = (bv_sub_orig[nb_rad] + db) % 4;
      bv_sub_mod.set(nb_rad, b_new);
      const int d_dist = dist_base_vecs_SW(bv_sub_mod, bvs_sub) - dist_orig;
      if (d_dist < sub.d_dist) {
        sub.d_dist = d_dist;
        sub.base = b_new;
      }      
    }
    if (time) *time += WallClockTime();

    // ---- decide whether to keep modification
    if (sub.d_dist < 0) {
      p_subs->push_back(sub);
      ib += nb_rad; // so that no changes overlap
    }
  }
}







void deletions_apply(const vec<Deletion> & dels, 
                     BaseVec * p_bv, 
                     int * p_ib1,
                     const vec<BaseVec> bvs)
{
  const size_t n_dels = dels.size();
  *p_ib1 -= n_dels;
  
  const int nb = p_bv->size();
  BaseVec bv_new(nb - n_dels);
  size_t i_del = 0;
  for (int ib = 0; ib != nb; ib++) {
    if (i_del != n_dels && 
        ib == dels[i_del].ib) 
      i_del++;
    else 
      bv_new.set(ib - i_del, (*p_bv)[ib]);
  }


  /*
  {
    vec<AlignmentIndexes> alis;
    for (size_t n_del = 0; n_del <= n_dels; n_del++) {
      BaseVec bv_tmp(nb - n_del);
      size_t i_del = 0;
      for (int ib = 0; ib != nb; ib++) {
        if (i_del != n_del && ib == dels[i_del].ib) 
          i_del++;
        else 
          bv_tmp.set(ib - i_del, (*p_bv)[ib]);
      }
      int dist = dist_base_vecs_SW(bv_tmp, bvs, & alis);
      cout << " n_del= " << n_del
           << " ib= " << (n_del ? dels[n_del - 1].ib : 0)
           << " dist= " << dist
           << endl;
    }
  }

  
  vec<int> ibs_del2(9);
  ibs_del2[0] = 441 + 12;
  ibs_del2[1] = 456 + 12;
  ibs_del2[2] = 457 + 12;
  ibs_del2[3] = 465 + 12;
  ibs_del2[4] = 476 + 12;
  ibs_del2[5] = 491 + 12;
  ibs_del2[6] = 495 + 12;
  ibs_del2[7] = 507 + 12;
  ibs_del2[8] = 543 + 12;
  
  {
    vec<AlignmentIndexes> alis;
    for (size_t n_del = 0; n_del <= ibs_del2.size(); n_del++) {
      BaseVec bv_tmp(nb - n_del);
      size_t i_del = 0;
      for (int ib = 0; ib != nb; ib++) {
        if (i_del != n_del && ib == ibs_del2[i_del]) 
          i_del++;
        else 
          bv_tmp.set(ib - i_del, (*p_bv)[ib]);
      }
      int dist = dist_base_vecs_SW(bv_tmp, bvs, & alis);
      cout << " n_del= " << n_del
           << " ib= " << (n_del ? ibs_del2[n_del - 1] : 0)
           << " dist= " << dist
           << endl;
    }
  }
  */

  *p_bv = bv_new;
}



void insertions_apply(const vec<Insertion> & inss, 
                      BaseVec * p_bv, 
                      int * p_ib1)
{
  const size_t n_inss = inss.size();
  *p_ib1 += n_inss;
  
  const int nb = p_bv->size();
  BaseVec bv_new(nb + n_inss);
  size_t i_ins = 0;
  for (int ib = 0; ib != nb; ib++) {
    if (i_ins != n_inss && ib == inss[i_ins].ib) {
      bv_new.set(ib + i_ins, inss[i_ins].base);
      i_ins++;
    }
    bv_new.set(ib + i_ins, (*p_bv)[ib]);
  }
  *p_bv = bv_new;
}



void substitutions_apply(const vec<Substitution> & subs, 
                         BaseVec * p_bv)
{
  const size_t n_subs = subs.size();
  const int nb = p_bv->size();
  size_t i_sub = 0;
  for (int ib = 0; ib != nb; ib++) {
    if (i_sub != n_subs && ib == subs[i_sub].ib) {
      p_bv->set(ib, subs[i_sub].base);
      i_sub++;
    }
  }
}





// =============================================
//  Tandem Repeats
// =============================================


class TandemRepeat
{
public:
  unsigned i0b;
  unsigned period;
  unsigned nb;
  
  TandemRepeat(const unsigned _i0b, 
               const unsigned _p, 
               const unsigned _nb) :
    i0b(_i0b), period(_p), nb(_nb) {}
  


  bool in_range(const unsigned i0b_range,
                const unsigned i1b_range)
  {
    ForceAssertGe(i1b_range, i0b_range);
    return (i0b_range + period < i0b + nb  &&  i1b_range - period > i0b);
  }  
  
};



void tandem_repeat_shrink(const unsigned n,
                          BaseVec * bv_p,
                          TandemRepeat * tr_p)
{
  const unsigned nb_left  = tr_p->i0b;  
  const unsigned nb_right = bv_p->size() - nb_left;

  const unsigned nb_del = n * tr_p->period;
  ForceAssertGe(tr_p->nb, nb_del);

  tr_p->nb -= nb_del;

  *bv_p = Cat(BaseVec(*bv_p,       0, nb_left),
              BaseVec(*bv_p, nb_left + nb_del, nb_right - nb_del));
} 
  
void tandem_repeat_expand(const unsigned n,
                          BaseVec * bv_p,
                          TandemRepeat * tr_p)
{
  const unsigned nb_left  = tr_p->i0b;   
  const unsigned nb_right = bv_p->size() - nb_left;

  const unsigned nb_ins = n * tr_p->period;
  BaseVec bv_ins(nb_ins);
  for (unsigned ib = 0; ib < nb_ins; ib++)
    bv_ins.set(ib, (*bv_p)[tr_p->i0b + ib % tr_p->period]);

  tr_p->nb += nb_ins;

  *bv_p = Cat(BaseVec(*bv_p, 0, nb_left),
              bv_ins,
              BaseVec(*bv_p, nb_left, nb_right));
} 








void tandem_repeats_print(const BaseVec & bv, 
                          const vec<TandemRepeat> & trs)
{
  const unsigned ntr = trs.size();
  unsigned itr = 0;
  bool is_repeat = false;
  for (unsigned ib = 0; ib < bv.size(); ib++) {
    if (itr < ntr && ib == trs[itr].i0b) is_repeat = true;
    if (itr < ntr && ib == trs[itr].i0b + trs[itr].nb) { is_repeat = false; itr++; }
    if (is_repeat) cout << as_base(bv[ib]);
    else           cout << ".";
    
  }     
}





void tandem_repeats_find(const BaseVec & bv, 
                         vec<TandemRepeat> * trs_p)
{
  trs_p->clear();
  const unsigned nb = bv.size();
  const unsigned period_min = 1;
  const unsigned period_max = 12;
  const unsigned n_copies_min = 5;
  const unsigned nb_rep_min = 12;
  for (unsigned ib = 0; ib < nb; ib++) {
    bool found = false;
    unsigned nb_rep = 0;
    unsigned period = period_min;
    do {
      unsigned nb_min = period * n_copies_min;
      if (nb_rep_min > nb_min) 
        nb_min = nb_rep_min;
      if (ib + nb_min < nb) {
        nb_rep = 1;
        while (ib + nb_rep < nb && bv[ib + nb_rep] == bv[ib + nb_rep % period]) 
          nb_rep++;
        if (nb_rep >= nb_min) {
          found = true;
          trs_p->push_back(TandemRepeat(ib, period, nb_rep));
        }
      }
    } while (!found && ++period <= period_max);
    if (found)
      ib += nb_rep - period;
  }
}
  
  


void tandem_repeats_optimize(const size_t i_gap,
                             const vec<BaseVec> & bvs,
                             const int dist0,
                             const int i0b,
                             int * i1b_p,
                             BaseVec * bv_p, 
                             vec<TandemRepeat> * trs_p,
                             const unsigned verbosity)
{
  const size_t ntr = trs_p->size();
  int dist_lowest = dist0;
  int n_ins = 0;


  for (unsigned itr = 0; itr < ntr; itr++) {
    (*trs_p)[itr].i0b += n_ins;      // add previous insertions (deletions are negative)
    TandemRepeat tr = (*trs_p)[itr]; // local copy
    BaseVec bv = *bv_p;         // local copy
    int dist = dist0;
    vec<AlignmentIndexes> alis;
    
    // ---- try deletions

    unsigned i_pass = 0;
    bool improved = false;
    do {
      if (tr.in_range(i0b, *i1b_p)) {

        tandem_repeat_shrink(1, &bv, &tr); 
        dist = dist_base_vecs_SW(bv, bvs, & alis, verbosity);

        if (dist < dist_lowest) {
          *bv_p = bv;
          (*trs_p)[itr] = tr;
          *i1b_p -= tr.period;
          n_ins  -= tr.period;
          dist_lowest = dist;
          improved = true;
        }
        if (verbosity >= 2) 
          #pragma omp critical
          cout << "i_gap= "   << i_gap 
               << " i_del= " << i_pass
               << " dist_lowest= " << dist_lowest
               << " dist= " << dist
               << " improved= " << improved
               << endl;
      }

    } while (i_pass++ < 2 && dist == dist_lowest);


    if (!improved) {   // ---- try insertions
      
      i_pass = 0;
      do {
        if (tr.in_range(i0b, *i1b_p)) {

          tandem_repeat_expand(1, &bv, &tr); 
          dist = dist_base_vecs_SW(bv, bvs, & alis, verbosity);

          if (dist < dist_lowest) {
            *bv_p = bv;
            (*trs_p)[itr] = tr;
            *i1b_p += tr.period;
            n_ins  += tr.period;
            dist_lowest = dist;
          }
          if (verbosity >= 2) 
            #pragma omp critical
            cout << "i_gap= "   << i_gap 
                 << " i_ins= " << i_pass
                 << " dist_lowest= " << dist_lowest
                 << " dist= " << dist
                 << " improved= " << improved
                 << endl;
        }
        
      } while (i_pass++ < 2 && dist < dist_lowest);
    }

  }
}
                             

























// ===================================
// patcher optimizer
// ===================================




void patcher_optimal(const BaseVecVec & unibases,
                     const vec<GapPatcher> & patchers,
                     const size_t ip_best, 
                     const int L,
                     const int nb_rad_SW,
                     const int sz_padding_min,
                     const unsigned i_gap,
                     GapPatcher0 * patcher0_opt_p,
                     const unsigned PATCH_VERBOSITY,
                     vec<double> * timers_p) 
{
  double & time_tot = (*timers_p)[0];
  double & time_short = (*timers_p)[1];
  double & time_large = (*timers_p)[2];
  time_tot = -WallClockTime();
  time_short = 0;
  time_large = 0;

  const size_t ibv1 = patchers[0].t1;
  const size_t ibv2 = patchers[0].t2;

  // ---- Test that all patchers refer to the same contigs
  const size_t n_patchers = patchers.size();
  for (size_t ip = 1; ip < n_patchers; ip++) {
    ForceAssertEq(patchers[0].t1, patchers[ip].t1);
    ForceAssertEq(patchers[0].t2, patchers[ip].t2);
  }

  const BaseVec & bv1 = unibases[ibv1];
  const BaseVec & bv2 = unibases[ibv2];

  const unsigned np = patchers.size();

  // ---- Find minimum tpos1 and maximum tpos2 and validate
  int tpos1_min = patchers[0].tpos1;
  int tpos2_max = patchers[0].tpos2;
  for (unsigned ip = 0; ip != np; ip++) {
    if (patchers[ip].tpos1 < tpos1_min) tpos1_min = patchers[ip].tpos1;
    if (patchers[ip].tpos2 > tpos2_max) tpos2_max = patchers[ip].tpos2;
  }

  ForceAssertGe(tpos1_min,               sz_padding_min);
  ForceAssertGe(bv2.isize() - tpos2_max, sz_padding_min);

  // ---- Expand by L and nb_rad_SW and validate
  tpos1_min -= nb_rad_SW + L;
  ForceAssertGt(tpos1_min, 0);

  tpos2_max += nb_rad_SW + L;
  ForceAssertLt(tpos2_max, bv2.isize());
  
  // ---- Output a whole bunch of stuff
  if (PATCH_VERBOSITY >= 1) {
    #pragma omp critical
    {
      cout << "i_gap= " << i_gap   // this allows grep for 'i_gap= <number>' 
           << "  ibv1= " << setw(10) << ibv1 
           << "  size= " << setw(10) << bv1.size() 
           << "  tpos1_min= " << setw(10) << tpos1_min << endl;
      cout << "i_gap= " << i_gap  // this allows grep for 'i_gap= <number>' 
           << "  ibv2= " << setw(10) << ibv2
           << "  size= " << setw(10) << bv2.size() 
           << "  tpos2_max= " << setw(10) << tpos2_max << endl;
      for (unsigned ip = 0; ip != np; ip++) {
        cout << "i_gap= " << i_gap  // this allows grep for 'i_gap= <number>' 
             << "  ip= " << setw(3) << ip
             << "  rid= "   << setw(10) << patchers[ip].rid 
             << "  rid/2= " << setw(10) << (patchers[ip].rid/2) 
             << "  tpos1= " << setw(6) << patchers[ip].tpos1 
             << "  tpos2= " << setw(6) << patchers[ip].tpos2 
             << "  r.size()= " << setw(8) << patchers[ip].r.size() 
             << "  gap= " << setw(6) 
             << (int(patchers[ip].r.size()) + patchers[ip].tpos1 - int(bv1.size()) - patchers[ip].tpos2)
             << (ip == ip_best ? "   BEST   " : "")
             << endl;
        cout << "i_gap= " << i_gap << " patcher= " << patchers[ip].r.ToString() << endl;
      }
      if (PATCH_VERBOSITY >= 3) {
        cout << "i_gap= " << i_gap << " bv1= " << bv1 << endl;
        cout << "i_gap= " << i_gap << " bv2= " << bv2 << endl;
      }
    }
  }
 

  // ---- Extend all BaseVecs to tpos1_min and tpos2_max
  vec<BaseVec> bvs(np);
  for (unsigned ip = 0; ip != np; ip++)
    bvs[ip] = Cat(BaseVec(bv1, tpos1_min, patchers[ip].tpos1 - tpos1_min), 
                  patchers[ip].r,
                  BaseVec(bv2, patchers[ip].tpos2, tpos2_max - patchers[ip].tpos2));


  // ---- Determine the consensus (best) read
  BaseVec bv_best = bvs[ip_best];
  const int nb_best = bv_best.size();

  // ---- ib0_best is the index, in bv_best, of the 1st base of the patching sequence
  int ib0_best = patchers[ip_best].tpos1 - tpos1_min;
  ForceAssertGe(ib0_best, nb_rad_SW);

  // ---- ib1_best is the index, in bv_best, of the last base of the patching sequence
  int ib1_best = nb_best - (tpos2_max - patchers[ip_best].tpos2);
  ForceAssertLe(ib1_best, nb_best - nb_rad_SW - 1);


  // ---- Compute all full alignments and initial distance
  vec<AlignmentIndexes> alis;
  int dist = dist_base_vecs_SW(bv_best, bvs, &alis);
  int dist_new = dist;
  
  if (PATCH_VERBOSITY >= 1) 
    #pragma omp critical
    cout << "i_gap= " << i_gap 
	 << " ib0_best= " << setw(5) << ib0_best 
	 << " ib1_best= " << setw(5) << ib1_best 
	 << " i_pass= -1 dist0= " << dist 
         << endl; 

  
  // ---- Declare some stuff
  vec<Deletion> dels;
  vec<Insertion> inss;
  vec<Substitution> subs;
  unsigned n_corr;
  
  const unsigned n_passes = 60;
  unsigned i_pass = 0;
  
  do {
    dist = dist_new;
    n_corr = 0;
  
    // -- Find deletions and correct them
    deletions_find(bv_best, ib0_best, ib1_best, bvs, alis, nb_rad_SW, & dels, 
                   & time_short, PATCH_VERBOSITY);
    deletions_apply(dels, & bv_best, & ib1_best, bvs);

    // -- Score corrections and get new alignments
    time_large -= WallClockTime();
    int dist_del = dist_new = dist_base_vecs_SW(bv_best, bvs, & alis);
    time_large += WallClockTime();



    // -- Find insertions and correct them
    insertions_find(bv_best, ib0_best, ib1_best, bvs, alis, nb_rad_SW, & inss, & time_short);
    insertions_apply(inss, & bv_best, & ib1_best);


    // -- Score corrections and get new alignments
    time_large -= WallClockTime();
    int dist_ins = dist_new = dist_base_vecs_SW(bv_best, bvs, & alis);
    time_large += WallClockTime();
    


    // -- Find substitutions and correct them
    substitutions_find(bv_best, ib0_best, ib1_best, bvs, alis, nb_rad_SW, & subs, & time_short);
    substitutions_apply(subs, & bv_best);

    
    // -- Score corrections and get new alignments
    time_large -= WallClockTime();
    int dist_sub = dist_new = dist_base_vecs_SW(bv_best, bvs, & alis);
    time_large += WallClockTime();


    n_corr = dels.size() + inss.size() + subs.size();


    if (PATCH_VERBOSITY >= 1)
      #pragma omp critical
      cout << "i_gap= "   << i_gap
           << " i_pass= " << setw(4) << i_pass 
           << " dist= "   << setw(5) << dist_new
	   << " ib0_best= " << setw(5) << ib0_best 
	   << " ib1_best= " << setw(5) << ib1_best 
           << " n_dels= " << setw(5) << dels.size() 
           << " "         << setw(5) << dist_del - dist
           << " n_inss= " << setw(5) << inss.size() 
           << " "         << setw(5) << dist_ins - dist_del
           << " n_subs= " << setw(5) << subs.size() 
           << " "         << setw(5) << dist_sub - dist_ins
           << " n_corr= " << setw(3) << n_corr
	   << endl;

    i_pass++;

  } while (n_corr   != 0 && 
           dist_new <  dist && 
           i_pass   <= n_passes);


  // ---- Tandem repeats testing

  if (true) {
    vec<TandemRepeat> repeats;
    tandem_repeats_find(bv_best, & repeats);
    if (PATCH_VERBOSITY >= 1) {
      #pragma omp critical
      {
        cout << "i_gap= "   << i_gap << " n_repeats= " << repeats.size() << " tandem= ";
        tandem_repeats_print(bv_best, repeats);
        cout << endl;
      }
    }
    
    tandem_repeats_optimize(i_gap, bvs, dist_new, 
                            ib0_best, & ib1_best, & bv_best, & repeats,
                            PATCH_VERBOSITY);
  }


  // ---- Trim bv_best down to its original start and end

  *patcher0_opt_p = GapPatcher0(BaseVec(bv_best, ib0_best, ib1_best - ib0_best), 
                                tpos1_min + ib0_best, 
                                tpos2_max - (bv_best.isize() - ib1_best));


  if (PATCH_VERBOSITY >= 1) { 
    #pragma omp critical
    cout << "i_gap= "   << i_gap
         << " patcher_opt= " << patcher0_opt_p->r << endl;
    if (PATCH_VERBOSITY >= 3) { 
      #pragma omp critical
      patcher_short_print(i_gap, *patcher0_opt_p, bv1, bv2);
    }
  }


  time_tot += WallClockTime();

}























