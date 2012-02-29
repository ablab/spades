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



#include <math.h>

#include "paths/OffsetDistribution.h"



void bridges_output(const vec<IntDistribution> & dists, 
                    const vec<GapBridge> & bridges)
{
  const size_t nd = dists.size();
  const size_t nb = bridges.size();

  for (size_t ib = 0; ib != nb; ib++) {
    const GapBridge & b = bridges[ib];
    
    cout << "lib" << b.i_lib()
	 << " fw" << b.tig_fw() 
	 << " bw" << b.tig_bw()   
	 << " " << setw(10) << b.size_contig1() 
	 << " " << setw(10) << b.size_contig2() 
	 << " " << setw(10) << b.i0_fw() 
	 << " " << setw(10) << b.i0_bw() 
	 << endl;
  }
}






void offset_range_update(const vec<IntDistribution> & dists,
                         const int sz_tig1,
                         const int sz_tig2,
                         const vec<int> & sz_fw,
                         const vec<int> & sz_bw,
                         int * p_os_min,
                         int * p_os_max) 
{
  const size_t n_libs = dists.size();
  for (size_t i_l = 0; i_l != n_libs; i_l++) {

    const int sz_tigs = sz_tig1 + sz_tig2;
    const int sz_reads = sz_fw[i_l] + sz_bw[i_l];
    
    const int l_min = dists[i_l].x_min();
    const int l_max = dists[i_l].x_max();
    
    //cout << "l_range = " << l_max - l_min << endl; 
    
    // the offset range for l_max
    const int os12_max =   l_max + sz_tig1 - sz_reads;
    const int os21_max = -(l_max + sz_tig2 - sz_reads);
    
    //cout << "l_max = " << setw(6) << l_max 
    //     << "  os12 = " << setw(6) << os12_max
    //     << "  os21 = " << setw(6) << os21_max
    //     << endl;
    
    // the gap range for l_min
    const int os12_min =   l_min - sz_tig2;
    const int os21_min = -(l_min - sz_tig1);
    
    //cout << "l_min = " << setw(6) << l_min 
    //     << "  os12 = " << setw(6) << os12_min
    //     << "  os21 = " << setw(6) << os21_min
    //     << endl;
    
    
    const int os_max = (os12_max > os21_min) ? os12_max : os21_min;
    const int os_min = (os21_max < os12_min) ? os21_max : os12_min;
    
    if (*p_os_min > os_min) *p_os_min = os_min;
    if (*p_os_max < os_max) *p_os_max = os_max;
  }

  //cout << "offset range = [" << *p_os_min
  //     << ", " << *p_os_max << "]"
  //     << "  " << *p_os_max - *p_os_min
  //     << endl;
}


void log_dists_bw_not_tig2_compute(const vec<IntDistribution> & dists,
                                   const vec<int> & szs_fw,
                                   const vec<int> & szs_bw,
                                   const int sz_tig1,
                                   const int sz_tig2,
                                   const int os_min,
                                   const int os_max,
                                   vec<IntLogDistribution> * p_log_dists)
{
  const size_t n_libs = dists.size();

  for (size_t i_l = 0; i_l != n_libs; i_l++) {
    const IntDistribution & dist = dists[i_l];
    //dist.to_text_file("dist_not2.l" + ToString(i_l));
    const int sz_fw = szs_fw[i_l];
    const int sz_bw = szs_bw[i_l];

    const int i_fw_max = sz_tig1 - sz_fw;
    const int i_fw_min = 0;
    
    const int x_min = os_min - i_fw_max;
    const int x_max = os_max - i_fw_min;
       
    IntFunction<double> d_not2(x_min, x_max);
    for (int x = x_min; x <= x_max; x++)
      d_not2[x] = 1.0 - dist.prob_in(x + sz_bw, x + sz_tig2);

    p_log_dists->push_back(d_not2);

    //d_not2.to_text_file("d_not2.l" + ToString(i_l));
    //p_log_dists->back().to_text_file("d_not2.log.l" + ToString(i_l));
  }
}


void log_dists_bw_not_tig1_compute(const vec<IntDistribution> & dists,
                                   const vec<int> & szs_fw,
                                   const vec<int> & szs_bw,
                                   const int sz_tig1,
                                   const int sz_tig2,
                                   const int os_min,
                                   const int os_max,
                                   vec<IntLogDistribution> * p_log_dists)
{
  const size_t n_libs = dists.size();

  for (size_t i_l = 0; i_l != n_libs; i_l++) {
    const IntDistribution & dist = dists[i_l];
    const int sz_fw = szs_fw[i_l];
    const int sz_bw = szs_bw[i_l];
 
    const int i_fw_max = sz_tig2 - sz_fw;
    const int i_fw_min = 0;
    
    const int x_min = - os_max - i_fw_max;
    const int x_max = - os_min - i_fw_min;
    
    IntFunction<double> d_not1(x_min, x_max);
    for (int x = x_min; x <= x_max; x++)
      d_not1[x] = 1.0 - dist.prob_in(x + sz_bw, x + sz_tig1);

    p_log_dists->push_back(d_not1);
      
    //d_not1.to_text_file("1.txt");
    //p_log_dists->back().to_text_file("log1.txt");
  }
}





IntDistribution offset_distribution_compute(const vec<IntDistribution> & dists, 
                                            const vec<GapBridge> & bridges,
                                            ostream * p_log,
                                            const int flag)
{
  typedef vec<GapBridge>::const_iterator BridgeIter;


  ofstream devnull("/dev/null");
  ostream &ostrm = p_log ? *p_log : devnull;

  //if (p_log) output_bridges(dists, bridges);

  const size_t n_libs = dists.size();

  const int sz_tig1 = bridges.begin()->size_contig1();
  const int sz_tig2 = bridges.begin()->size_contig2();

  // ---- obtain the mean read lengths for fw and bw reads for each library
  vec<int> sz_mean_fw_lib(n_libs, 0);
  vec<int> sz_mean_bw_lib(n_libs, 0);
  vec<int> n_bridges_lib(n_libs, 0);
  for (BridgeIter it_b = bridges.begin(); it_b != bridges.end(); it_b++) {
    ForceAssertEq(it_b->size_contig1(), sz_tig1);  
    ForceAssertEq(it_b->size_contig2(), sz_tig2);
    const size_t i_l = it_b->i_lib();
    n_bridges_lib[i_l]++;
    sz_mean_fw_lib[i_l] += it_b->size_fw_read();
    sz_mean_bw_lib[i_l] += it_b->size_bw_read();
  }
  for (size_t i_l = 0; i_l != n_libs; i_l++) {
    sz_mean_fw_lib[i_l] /= n_bridges_lib[i_l];
    sz_mean_bw_lib[i_l] /= n_bridges_lib[i_l];
  }


  // ---- Find the maximum possible offset from all the distributions.
  int os_min = 0;
  int os_max = 0;
  offset_range_update(dists, sz_tig1, sz_tig2,
                      sz_mean_fw_lib, sz_mean_bw_lib,
                      & os_min, &os_max);
  





  // ---- Logarithms of invariant size probability
  vec<IntLogDistribution> log_dists_inv;
  for (size_t i_l = 0; i_l != n_libs; i_l++)
    log_dists_inv.push_back(dists[i_l]);
  

  //dists[0].to_text_file("0.txt");
  //log_dists_inv.back().to_text_file("log0.txt");
  
  // ---- Logarithms of probability of starting on tig 1 and NOT ending on tig 2
  vec<IntLogDistribution> log_dists_bw_not_tig2;
  log_dists_bw_not_tig2_compute(dists, sz_mean_fw_lib, sz_mean_bw_lib,
                                sz_tig1, sz_tig2,
                                os_min, os_max,
                                &log_dists_bw_not_tig2);


  // ---- Logarithms of probability of starting on tig 2 and NOT ending on tig 1
  vec<IntLogDistribution> log_dists_bw_not_tig1;
  log_dists_bw_not_tig1_compute(dists, sz_mean_fw_lib, sz_mean_bw_lib,
                                sz_tig1, sz_tig2,
                                os_min, os_max,
                                &log_dists_bw_not_tig1);

  /*
  for (size_t i_l = 0; i_l != n_libs; i_l++) {
    log_dists_inv[i_l].to_text_file("logd_inv.l" + ToString(i_l));
    log_dists_bw_not_tig2[i_l].to_text_file("logd_bw_n2.l" + ToString(i_l));
    log_dists_bw_not_tig1[i_l].to_text_file("logd_bw_n1.l" + ToString(i_l));
  }
  */




  // ---- The prior distribution.  Assume uniform.
  ostrm << "  os: [" << os_min << ", " << os_max << "] " << flush;
  vec< vec<IntLogDistribution> > ldos(n_libs, 
				      vec<IntLogDistribution>
				      (4, IntLogDistribution(os_min, os_max, 0.0)));
  

  // ---- Go through bridges and update offset distribution.
  size_t i_b = 0;
  const size_t n_b = bridges.end() - bridges.begin();
  vec<size_t> ntype(4,0);
  ostrm << " " << flush;

  for (BridgeIter it_b = bridges.begin(); it_b != bridges.end(); it_b++) {
    if ((10 * i_b) / n_b != (10 * (i_b+1) / n_b))
      ostrm << "." << flush;
    

    const size_t i_l = it_b->i_lib();
    
    bool done = false;

    if (it_b->bridges_fw1_bw2() && (!flag || flag == 1)) {
      ntype[0]++;
      const IntLogDistribution & log_dist = log_dists_inv[i_l];
      const int delta = it_b->i0_bw() + 1 - it_b->i0_fw();    // verified!
      for (int os = os_min; os <= os_max; os++) 
        ldos[i_l][0][os] += log_dist[delta + os];
      done = true;
    }

    //if ((it_b->bridges_fw1_bw0() || it_b->bridges_fw1_bw1()) && (!flag || flag == 2)) {
    if (it_b->bridges_fw1_bw0() && (!flag || flag == 2)) {
      const IntLogDistribution & log_dist = log_dists_bw_not_tig2[i_l];
      const int delta = - it_b->i0_fw();  // verified!
      if (1) {
	ntype[1]++;
        for (int os = os_min; os <= os_max; os++) 
          ldos[i_l][1][os] += log_dist[delta + os];
      }
      done = true;
    }


    if (it_b->bridges_fw2_bw1() && (!flag || flag == 3)) {
      ntype[2]++;
      const IntLogDistribution & log_dist = log_dists_inv[i_l];
      const int delta = it_b->i0_bw() + 1 - it_b->i0_fw();    // verified!
      for (int os = os_min; os <= os_max; os++) 
        ldos[i_l][2][os] += log_dist[delta - os];
      done = true;
    }
  

    //if ((it_b->bridges_fw2_bw0() || it_b->bridges_fw2_bw2()) && (!flag || flag == 4)) {
    if (it_b->bridges_fw2_bw0() && (!flag || flag == 4)) {
      const IntLogDistribution & log_dist = log_dists_bw_not_tig1[i_l];
      const int delta = - it_b->i0_fw();  // verified!
      if (1) {
	ntype[3]++;
        for (int os = os_min; os <= os_max; os++) 
          ldos[i_l][3][os] += log_dist[delta - os];
      }
      done = true;
    }
    i_b++;
    //ForceAssert(i_b <= 10); 
  }
  ostrm << " [" << ntype[0] 
        << ", " << ntype[1] 
        << ", " << ntype[2] 
        << ", " << ntype[3] 
        << "]" << flush;
    
  /* 
    for (size_t i_l = 0; i_l != n_libs; i_l++)
      for (size_t j = 0; j != 4; j++)
	ldos[i_l][j].to_text_file("logd.l" + ToString(i_l) + "." + ToString(j));
  */
				    
  IntLogDistribution log_dist_offset(os_min, os_max, 0.0);
  for (size_t i_l = 0; i_l != n_libs; i_l++)
    for (size_t j = 0; j != 4; j++)
      for (int os = os_min; os <= os_max; os++) 
	log_dist_offset[os] += ldos[i_l][j][os];

  return log_dist_offset.int_distribution();
}
















IntDistribution distribution_bridge_fw_given_offset(const IntDistribution dist_inv,
                                                    const int sz_tig1,
                                                    const int sz_tig2,
                                                    const int sz_fw,
                                                    const int sz_bw,
                                                    const int os_min,
                                                    const int os_max,
                                                    const int flag)
{
  IntFunction<double> func(os_min, os_max);

  if (flag == 1) {   // fw1 -> bw2
    for (int os = os_min; os <= os_max; os++) {
      const int j0_bw_min = os + sz_bw   - 1;
      const int j0_bw_max = os + sz_tig2 - 1;
      const int j0_fw_min = 0;
      const int j0_fw_max = sz_tig1 - sz_fw;
      double p = 0.0;
      for (int j0_fw = j0_fw_min; j0_fw <= j0_fw_max; j0_fw++) 
        p += dist_inv.prob_in(j0_bw_min - j0_fw + 1, j0_bw_max - j0_fw + 1);
      
      func[os] = p;  // no need to divide by (sz_tig1 - sz_fw + 1), just a constant
    }
  }
  else if (flag == 2) { // fw1 -> ?
    for (int os = os_min; os <= os_max; os++) {
      const int j0_bw_min = os + sz_bw   - 1;
      const int j0_bw_max = os + sz_tig2 - 1;
      const int j0_fw_min = 0;
      const int j0_fw_max = sz_tig1 - sz_fw;
      double p = 0.0;
      for (int j0_fw = j0_fw_min; j0_fw <= j0_fw_max; j0_fw++) 
        p += dist_inv.prob_not_in(j0_bw_min - j0_fw + 1, j0_bw_max - j0_fw + 1);
      
      func[os] = p;  // no need to divide by (sz_tig1 - sz_fw + 1), just a constant
    }
  }
  else if (flag == 3) {  // fw2 -> bw1
    for (int os = os_min; os <= os_max; os++) {
      const int j0_bw_min = sz_bw   - 1;
      const int j0_bw_max = sz_tig1 - 1;
      const int j0_fw_min = os;
      const int j0_fw_max = os + sz_tig2 - sz_fw;
      double p = 0.0;
      for (int j0_fw = j0_fw_min; j0_fw <= j0_fw_max; j0_fw++) 
        p += dist_inv.prob_in(j0_bw_min - j0_fw + 1, j0_bw_max - j0_fw + 1);
      
      func[os] = p;  // no need to divide by (sz_tig2 - sz_fw + 1), just a constant
    }
  } 
  else if (flag == 4) {  // fw2 -> ?
    for (int os = os_min; os <= os_max; os++) {
      const int j0_bw_min = sz_bw   - 1;
      const int j0_bw_max = sz_tig1 - 1;
      const int j0_fw_min = os;
      const int j0_fw_max = os + sz_tig2 - sz_fw;
      double p = 0.0;
      for (int j0_fw = j0_fw_min; j0_fw <= j0_fw_max; j0_fw++) 
        p += dist_inv.prob_not_in(j0_bw_min - j0_fw + 1, j0_bw_max - j0_fw + 1);
      
      func[os] = p;  // no need to divide by (sz_tig2 - sz_fw + 1), just a constant
    }
  }
  else {
    ForceAssert(0 == 1);
  }

  // here I rely on the IntDistribution constructor that accepts an IntFunction
  return func;
}
  








