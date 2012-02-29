///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// author: Filipe Ribeiro      04/2011
// 
//
//

#ifndef _PATHS__OFFSET_DISTRIBUTION_H
#define _PATHS__OFFSET_DISTRIBUTION_H

#include "MainTools.h"


#include "math/IntDistribution.h"


//  ---- EXTERNAL (to class ContigReadIndex) read positions:
//
//       FORWARD read         |------------->
//       in [i0, i1]           i0         i1
//
//       BACKWARD read        <-------------|
//       in [i1, i0]           i1         i0
//
//
//  ---- INTERNAL (to class ContigReadIndex) read positions (represented by 'j0'): 
// 
//                            |---- contig 1 ----|  gap  |---- contig 2 ----|
//   ContigReadIndex._j{0,1} : 0 1 2 ...                  0 1 2 ...

//
// IMPORTANT: Do NOT call the ContigReadIndex() constructor directly.
//   Call, instead, the static function corresponding to the type of bridge,
//   using the EXTERNAL read positions (see above):
//
//     ContigReadIndex::contig1(sz_tig1, i0, i1);
//
//     ContigReadIndex::contig2(i0, i1);
//
//     ContigReadIndex::not_aligned(read_len);
//


class ContigReadIndex
{
  uint64_t _j0      : 31; // index of first read base on contig.
  uint64_t _j1      : 31; // index of last read base on contig.
  uint64_t _aligned :  2; // which contig it aligns to (0: not aligned). 

public:
  ContigReadIndex() : _j0(0), _j1(0), _aligned(3) {}
  ContigReadIndex(const int j0, 
                  const int j1,
                  const int aligned) 
    : _j0(j0), _j1(j1), _aligned(aligned) {}

  int  i_base0()        const { return _j0; }
  int  i_base1()        const { return _j1; }
  bool aligns()         const { return (_aligned == 1 || _aligned == 2); }
  bool aligns_tig1()    const { return (_aligned == 1); }
  bool aligns_tig2()    const { return (_aligned == 2); }
  bool is_forward()     const { return (_j1 > _j0); }
  bool is_backward()    const { return (_j1 < _j0); }
  int  read_length()    const { return (is_forward()) ? _j1 - _j0 + 1 : _j0 - _j1 + 1; }
  int  tig()            const { return _aligned; }

  static ContigReadIndex contig1(const int sz_tig1,
                                 const int i0, 
                                 const int i1) 
  {
    ForceAssertLt(i0, sz_tig1);
    ForceAssertLt(i1, sz_tig1);
    ForceAssertGe(i0, 0);
    ForceAssertGe(i1, 0);
    return ContigReadIndex(i0, i1, 1);
  }

  static ContigReadIndex contig2(const int sz_tig2,
                                 const int i0,
                                 const int i1) 
  {
    ForceAssertLt(i0, sz_tig2);
    ForceAssertLt(i1, sz_tig2);
    ForceAssertGe(i0, 0);
    ForceAssertGe(i1, 0);
    return ContigReadIndex(i0, i1, 2); 
  }

  static ContigReadIndex not_aligned(const int len) 
  {
    return ContigReadIndex(0, len, 0);
  }
};











// IMPORTANT: Examples of GapBridge declarations:
//
//
//  1. FORWARD  read    contig 1     [i0_f, i1_f] 
//     BACKWARD read    contig 2     [i1_b, i0_b]
//
//       const GapBridge gl(i_dist, sz_tig1, sz_tig2,
//                          ContigReadIndex::contig1(sz_tig1, i0_f, i1_f), 
//                          ContigReadIndex::contig2(sz_tig2, i0_b, i1_b));
//
//
//  2. FORWARD  read    contig 2     [i0_f, i1_f]
//     BACKWARD read    not aligned  (must still provide read length)
//
//       const GapBridge gl(i_dist, sz_tig1, sz_tig2,
//                          ContigReadIndex::contig2(sz_tig2, i0_f, i1_f), 
//                          ContigReadIndex::not_aligned(len_b));
//
//
//  3. FORWARD  read    contig 1     [i0_f, i1_f] 
//     BACKWARD read    contig 1     [i1_b, i0_b]
//
//       const GapBridge gl(i_dist, sz_tig1, sz_tig2,
//                          ContigReadIndex::contig1(sz_tig1, i0_f, i1_f), 
//                          ContigReadIndex::contig1(sz_tig1, i0_b, i1_b));
//
//
//  4. FORWARD  read    not_aligned  (must still provide read length)
//     BACKWARD read    contig 1     [i1_b, i0_b]
//
//       const GapBridge gl(i_dist, sz_tig1, sz_tig2,
//                          ContigReadIndex::not_aligned(len_f),
//                          ContigReadIndex::contig1(sz_tig1, i0_b, i1_b));
//
//
//
//  NOTES: 
//
//     - the forward read must always come first.
//
//     - 'i_dist' is the index of the invariant size distribution.  



class GapBridge
{
private:
  size_t                  _i_dist;   // distribution index 
  int                     _sz_tig1;
  int                     _sz_tig2;
  ContigReadIndex         _fw;       // FORWARD read alignment info
  ContigReadIndex         _bw;       // BACKWARD read alignment info
   
public:
  GapBridge(const size_t i_dist,
            const int sz_tig1,
            const int sz_tig2,
            const ContigReadIndex & fw, 
            const ContigReadIndex & bw) 
    : _i_dist(i_dist), 
      _sz_tig1(sz_tig1), 
      _sz_tig2(sz_tig2), 
      _fw(fw), 
      _bw(bw) 
  {
    if (_fw.aligns()) ForceAssert(_fw.is_forward());
    if (_bw.aligns()) ForceAssert(_bw.is_backward());
  }

  size_t i_lib() const { return _i_dist; }

  int size_contig1() const { return _sz_tig1; }
  int size_contig2() const { return _sz_tig2; }

  int i0_fw() const { return _fw.i_base0(); }
  int i0_bw() const { return _bw.i_base0(); }
  int i1_fw() const { return _fw.i_base1(); }
  int i1_bw() const { return _bw.i_base1(); }

  int size_fw_read() const { return _fw.read_length(); }
  int size_bw_read() const { return _bw.read_length(); }

  int tig_fw() const { return _fw.tig(); }
  int tig_bw() const { return _bw.tig(); }

  bool aligns_fw() const { return _fw.aligns(); }
  bool aligns_bw() const { return _bw.aligns(); }
  
  bool aligns_fw0() const { return !_fw.aligns(); }
  bool aligns_bw0() const { return !_bw.aligns(); }

  bool aligns_fw1() const { return _fw.aligns_tig1(); }
  bool aligns_bw1() const { return _bw.aligns_tig1(); }

  bool aligns_fw2() const { return _fw.aligns_tig2(); }
  bool aligns_bw2() const { return _bw.aligns_tig2(); }
  
  // 0: no alignment
  // 1: alignment to contig 1
  // 2: alignment to contig 2
  
  bool bridges_fw0_bw0()  const { return !_fw.aligns()       &&  !_bw.aligns(); }
  bool bridges_fw0_bw1()  const { return !_fw.aligns()       &&   _bw.aligns_tig1(); }
  bool bridges_fw0_bw2()  const { return !_fw.aligns()       &&   _bw.aligns_tig2(); }

  bool bridges_fw1_bw0()  const { return  _fw.aligns_tig1()  &&  !_bw.aligns(); }
  bool bridges_fw1_bw1()  const { return  _fw.aligns_tig1()  &&   _bw.aligns_tig1(); }
  bool bridges_fw1_bw2()  const { return  _fw.aligns_tig1()  &&   _bw.aligns_tig2(); }

  bool bridges_fw2_bw0()  const { return  _fw.aligns_tig2()  &&  !_bw.aligns(); }
  bool bridges_fw2_bw1()  const { return  _fw.aligns_tig2()  &&   _bw.aligns_tig1(); }
  bool bridges_fw2_bw2()  const { return  _fw.aligns_tig2()  &&   _bw.aligns_tig2(); }

  bool bridges_fw01_bw1() const { return !_fw.aligns_tig2()  &&   _bw.aligns_tig1(); }
  bool bridges_fw02_bw2() const { return !_fw.aligns_tig1()  &&   _bw.aligns_tig2(); }

  bool bridges_fw1_bw01() const { return  _fw.aligns_tig1()  &&   !_bw.aligns_tig2(); }
  bool bridges_fw2_bw02() const { return  _fw.aligns_tig2()  &&   !_bw.aligns_tig1(); }


  friend bool less_than(const GapBridge & a, const GapBridge & b)
  {
    if (a._i_dist < b._i_dist) return true;
    if (a._i_dist > b._i_dist) return false;

    if (   a.bridges_fw1_bw2()  &&  ! b.bridges_fw1_bw2()) return false;
    if ( ! a.bridges_fw1_bw2()  &&    b.bridges_fw1_bw2()) return false;
    if (   a.bridges_fw1_bw2()  &&    b.bridges_fw1_bw2()) {
      return ((a.i0_bw() - a.i0_fw()) > (b.i0_bw() - b.i0_fw()));
    }

    if (   a.bridges_fw2_bw1()  &&  ! b.bridges_fw2_bw1()) return false;
    if ( ! a.bridges_fw2_bw1()  &&    b.bridges_fw2_bw1()) return false;
    if (   a.bridges_fw2_bw1()  &&    b.bridges_fw2_bw1()) {
      return ((a.i0_bw() - a.i0_fw()) > (b.i0_bw() - b.i0_fw()));
    }

    if (   a.bridges_fw1_bw01()  &&  ! b.bridges_fw1_bw01()) return true;
    if ( ! a.bridges_fw1_bw01()  &&    b.bridges_fw1_bw01()) return false;
    if (   a.bridges_fw1_bw01()  &&    b.bridges_fw1_bw01()) {
      return (a.i0_fw() < b.i0_fw());
    }

    if (   a.bridges_fw2_bw02()  &&  ! b.bridges_fw2_bw02()) return true;
    if ( ! a.bridges_fw2_bw02()  &&    b.bridges_fw2_bw02()) return false;
    if (   a.bridges_fw2_bw02()  &&    b.bridges_fw2_bw02()) {
      return (a.i0_fw() < b.i0_fw());
    }

    if (   a.bridges_fw01_bw1()  &&  ! b.bridges_fw01_bw1()) return true;
    if ( ! a.bridges_fw01_bw1()  &&    b.bridges_fw01_bw1()) return false;
    if (   a.bridges_fw01_bw1()  &&    b.bridges_fw01_bw1()) {
      return (a.i0_bw() < b.i0_bw());
    }

    if (   a.bridges_fw02_bw2()  &&  ! b.bridges_fw02_bw2()) return true;
    if ( ! a.bridges_fw02_bw2()  &&    b.bridges_fw02_bw2()) return false;
    if (   a.bridges_fw02_bw2()  &&    b.bridges_fw02_bw2()) {
      return (a.i0_bw() < b.i0_bw());
    }

    return false; // got this far, then the bridges are equal.
  }

};





IntDistribution offset_distribution_compute(const vec<IntDistribution> & dists,
                                            const vec<GapBridge> & bridges,
                                            ostream * p_log = 0,
                                            const int flag = 0);









IntDistribution distribution_bridge_fw_given_offset(const IntDistribution dist_inv,
                                                    const int sz_tig1,
                                                    const int sz_tig2,
                                                    const int sz_fw,
                                                    const int sz_bw,
                                                    const int os_min,
                                                    const int os_max,
                                                    const int flag);





#endif
