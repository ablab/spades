///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Displays basic statistics for a fastb file.";

#include "Basevector.h"
#include "Vec.h"
#include "feudal/VirtualMasterVec.h"
//#include "math/Functions.h"
#include "system/ParsedArgs.h"
#include "system/RunTime.h"
#include "system/System.h"

#include "random/Random.h"

typedef VirtualMasterVec<BaseVec> VBaseVecVec;


void quick_stats(const String & fastb_fn) 
{
  size_t n_bv = MastervecFileObjectCount(fastb_fn);
  size_t objects =  MastervecFileRawCount(fastb_fn);
  size_t uncertainty = n_bv * 15;
  size_t max_size = objects * 16;
  cout << " total objects count: " << n_bv << "\n";
  cout << " total length estimate: " << (2 * max_size - uncertainty) / 2
       << " +/- " << uncertainty / 2
       << " (" << static_cast<double>(100 * uncertainty / 2) / (2 * max_size - uncertainty)
       << " %)" << "\n";
  cout << " total length range: [" << max_size - uncertainty << " to "
       << max_size << "]" << "\n\n";

  if (n_bv > 0) {
    BaseVecVec first_bv;
    first_bv.ReadOne(fastb_fn, 0);
    cout << " First object length: " << first_bv[0].size() << endl;
    cout << " Projected total length: " << first_bv[0].size() * n_bv << endl;
  }
}




class FastbStats
{
  bool _sampled;
  const VBaseVecVec & _vbvv;
  const size_t _n_bv;
  size_t _len0;
  const size_t _min_read_len;
  const size_t _max_read_len;

  size_t _n_bases;
  size_t _n_bv_empty;
  size_t _n_bv_short;
  size_t _n_bv_long;
  size_t _n_bv_good;

  size_t _min_len;
  size_t _max_len;
  vec<size_t> _lens;
  
  size_t _n_gc_bases;
  vec<size_t> _gc_spec;

public:
 
  FastbStats(const VBaseVecVec & vbvv, 
	     const size_t min_read_len,
	     const size_t max_read_len)
    : _sampled(true), 
      _vbvv(vbvv),
      _n_bv(_vbvv.size()),
      _len0(0),
      _min_read_len(min_read_len),
      _max_read_len(max_read_len),
      _n_bases(0),
      _n_bv_empty(0),
      _n_bv_short(0),
      _n_bv_long(0),
      _n_bv_good(0),
      _min_len(INT_MAX),
      _max_len(0),
      _n_gc_bases(0)
  {
    if (_n_bv == 0)
      cout << " total objects count: 0" << endl << endl;
  }

  

  void scale_stats()
  {
    const size_t n_bv_tot = _n_bv_empty + _n_bv_short + _n_bv_long + _n_bv_good;
    
    _n_bases    = 0.5 + (double(_n_bases)    / double(n_bv_tot)) * _n_bv;
    _n_gc_bases = 0.5 + (double(_n_gc_bases) / double(n_bv_tot)) * _n_bv;
    _n_bv_empty = 0.5 + (double(_n_bv_empty) / double(n_bv_tot)) * _n_bv;
    _n_bv_short = 0.5 + (double(_n_bv_short) / double(n_bv_tot)) * _n_bv;
    _n_bv_long  = 0.5 + (double(_n_bv_long)  / double(n_bv_tot)) * _n_bv;
    _n_bv_good  = 0.5 + (double(_n_bv_good)  / double(n_bv_tot)) * _n_bv;

  }



  void print_lengths_stats()
  {
    cout << " total objects count: " << _n_bv << endl;
    if (_sampled) 
      cout << " ---- estimated from sample ----" << endl;

    cout << " total length: " << _n_bases << endl;

    if (_lens.size() > 0) {
      sort(_lens.begin(), _lens.end());
      cout << " N50 length: " << N50(_lens) << endl;
    }

    cout << " avg length: " << fixed << setprecision(1) 
	 << double(_n_bases)/double(_n_bv_good) << endl;
    cout << " range: [" << _min_len << ", " << _max_len << "]" << endl;

    if (_min_read_len > 0) 
      cout << " reads shorter than " << _min_read_len << ": " << _n_bv_short 
	   << " (" << 100.0 * double(_n_bv_short)/double(_n_bv) << " %)" << endl;

    if (_max_read_len < INT_MAX) 
      cout << " reads longer than " << _max_read_len << ": " << _n_bv_long 
	   << " (" << 100.0 * double(_n_bv_long)/double(_n_bv) << " %)" << endl;

    if (_n_bv_empty > 0) 
      cout << " Warning! " << _n_bv_empty << " objects had zero length." << endl << endl;
  }  
  


  void print_gc_stats(const bool do_gc,
		      const bool do_gc_spec,
		      const String & gc_spec_head)
  {
    if (do_gc)
      cout << " GC content: " << fixed << setprecision(1) 
           << 100.0 * double(_n_gc_bases)/double(_n_bases) << " %" << endl; 
 
    if (gc_spec_head != "") {
      if (!do_gc_spec) {
        cout << " Warning! Can't compute GC spectrum because reads have different lengths." 
             << endl << endl;
      }
      else {
        const String gc_fn = gc_spec_head + ".gc_spec";
        ofstream os;
        os.open(gc_fn.c_str());

        os << "# 1:gc 2:freq 3:freq_cum 4:frac 5:frac_cum" << endl;
        os << fixed;
        
        size_t freq_cum = 0;
        double frac;
        double frac_cum = 0.0;
        
        for (unsigned i = 0; i != _gc_spec.size(); i++) {
          size_t freq = _gc_spec[i];
          freq_cum += freq;
          frac = double(freq) / double(_n_bv_good);
          frac_cum = double(freq_cum) / double(_n_bv_good);
          os << " " << setw(5) << i  
             << " " << setw(5) << freq
             << " " << setw(5) << freq_cum
             << " " << setw(5) << frac
             << " " << setw(5) << frac_cum
             << endl;
        }
        os.close();
        cout << "Wrote gc spectrum to " << gc_fn << endl;
      }
    }    
  }



  void add_base_vec(const BaseVec & bv, 
		    const bool size_only,
		    const bool do_gc,
		    bool & do_gc_spec) 
  {
    const size_t len = bv.size();
    
    if      (len == 0)            _n_bv_empty++;
    else if (len < _min_read_len) _n_bv_short++;
    else if (len > _max_read_len) _n_bv_long++;
    else {
      _n_bv_good++;
      _n_bases += len;

      if (len < _min_len) _min_len = len;
      if (len > _max_len) _max_len = len;
     
      if (!size_only) { 
	_lens.push_back(len);
	
	if (do_gc) {        
	  size_t n_gc = 0;
	  for (size_t i = 0; i != len; i++) {
	    const unsigned base = bv[i];
	    if (base == 1 || base == 2) n_gc++;
	  }
	  _n_gc_bases += n_gc;
	  if (do_gc_spec) {
	    if (_len0 == 0) {
	      _len0 = len;
	      _gc_spec.resize(_len0 + 1, 0);
	    }
	    if (len == _len0) _gc_spec[n_gc]++;
	    else              do_gc_spec = false; // different read lens => no GC spec
	  }
	}
      }
    }
  }
  

  void complete_stats(const bool size_only,
		      const bool do_gc,
		      const String & gc_spec_head)
  {
    _sampled = false;

    if (_n_bv > 0) {
      
      // ---- parse
      cout << "Streaming data:" << endl;

      _lens.reserve(_n_bv);
      bool do_gc_spec = (gc_spec_head != "");
      size_t i_bv = 0;
      for (VBaseVecVec::const_iterator it = _vbvv.begin(); 
	   it != _vbvv.end(); i_bv++, it++) {

	add_base_vec(*it, size_only, do_gc, do_gc_spec);
	dots_pct(i_bv, _n_bv);
      }
      cout << endl << endl;
 

      // ---- print
     
      print_lengths_stats();
      print_gc_stats(do_gc, do_gc_spec, gc_spec_head);
      cout << endl;
    }
  }




  void sample_indices(const size_t n_blocks,
		      const size_t n_bv_per_block,
		      vec<size_t> & indices)
  {
    const size_t n_indices = n_blocks * n_bv_per_block;
    indices.reserve(n_indices);

    const size_t block_size = _n_bv / n_blocks;
    const size_t d_ind = (block_size - n_bv_per_block)/2;

    for (size_t i = 0; i != n_blocks; i++)
      for (size_t j = 0; j != n_bv_per_block; j++) 
	indices.push_back(i * block_size + d_ind + j);
  }




  void sample_stats(const bool size_only,
		    const bool do_gc,
		    const String & gc_spec_head)
  {
    _sampled = true;

    if (_n_bv > 0) {
      const size_t n_blocks = 100;
      const size_t n_reads_per_block = 10000;
      const size_t n_indices = n_blocks * n_reads_per_block;

      if (n_indices > 0.8 * _n_bv) {
	complete_stats(size_only, do_gc, gc_spec_head);
      }
      else {

	vec<size_t> indices;
	sample_indices(n_blocks, n_reads_per_block, indices);


	// ---- parse
	cout << "Streaming data:" << endl;

	_lens.reserve(n_indices);
	bool do_gc_spec = (gc_spec_head != "");

	for (size_t i_ind = 0; i_ind != n_indices; i_ind++) {
	  add_base_vec(_vbvv[indices[i_ind]], size_only, do_gc, do_gc_spec);
	  dots_pct(i_ind, n_indices);
	}
	cout << endl << endl;

	scale_stats();

	// ---- print

	print_lengths_stats();
	print_gc_stats(do_gc, do_gc_spec, gc_spec_head);
	cout << endl;
      }
    }
  }

};




int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc
    (FASTB, "File for analysis");
  CommandArgument_Bool_OrDefault_Doc
    (QUICK, False, "Quick estimate of size without loading file");
  CommandArgument_Bool_OrDefault_Doc
    (SAMPLE, False, "Estimate stats by sample a subset of the data");
  CommandArgument_Bool_OrDefault_Doc
    (SIZE_ONLY, False, "Don't compute MIN, MAX or N50");
  CommandArgument_Int_OrDefault_Doc
    (MIN_READ_LEN, 0, "Only consider reads of length >= than MIN_READ_LEN");
  CommandArgument_Int_OrDefault_Doc
    (MAX_READ_LEN, INT_MAX, "Only consider reads of length <= than MAX_READ_LEN");
  CommandArgument_Bool_OrDefault_Doc
    (GC, False, "if True, GC content is computed.");
  CommandArgument_String_OrDefault_Doc
    (GC_SPEC, "", "if specified, GC spectrum goes to <GC_SPEC>.gc_spec.");
  EndCommandArguments;

  // Figure out the correct filename if the fastb extension isn't given

  if ( !IsRegularFile(FASTB) ) {
    if ( FASTB.EndsWith(".") )
      FASTB += "fastb";
    else if ( !FASTB.EndsWith(".fastb") )
      FASTB += ".fastb";
    if (!IsRegularFile(FASTB) )
      FatalErr("Could not find file: " + FASTB);
  }


  if (QUICK) { // 'n' dirty
    quick_stats(FASTB);
  }
  else {

    // Use VirtualMasterVec to avoid loading in the full BaseVecVec
    VBaseVecVec vbvv(FASTB.c_str());

    FastbStats fs(vbvv, MIN_READ_LEN, MAX_READ_LEN);
    if (SAMPLE) {
      fs.sample_stats(SIZE_ONLY, GC = true, GC_SPEC);
    }
    else {
      fs.complete_stats(SIZE_ONLY, GC, GC_SPEC);
    }
  }
}
