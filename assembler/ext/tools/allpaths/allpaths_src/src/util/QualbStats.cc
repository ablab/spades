///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
 * QualbStats
 *
 *
 * Filipe Ribeiro 2011-01 
 * 
 ******************************************************************************/


#include "MainTools.h"

#include "Qualvector.h"
#include "Basevector.h"


#define NQ 256

typedef VirtualMasterVec<QualVec> VQualVecVec;



float error_probability_from_quality(const uint8_t q)
{
  static vec<float> prob(NQ, 0.0);
  static const float factor =  -0.1 * log(10.0);

  if (prob[q] == 0.0) 
    prob[q] = exp(factor * q);
  return prob[q];
}


float quality_from_error_probability(const float pe)
{
  return - 10.0 * log(pe) / log(10.0); 
}  






size_t quantile_from_freqs(const vec<size_t> & freq,
			   const float quantile)
{
  const size_t n = freq.size();
  size_t tot = 0;
  for (size_t i = 0; i != n; i++)
    tot += freq[i];

  const size_t Nq = 0.5 + tot * quantile;
  size_t F = 0;        // cumulative of freq[i]
  for (size_t i = 0; i != n; i++) {
    F += freq[i];
    if (F >= Nq) return i;
  }
  return n;
}


float mean_from_freqs(const vec<size_t> & freq)
{
  float tot = 0.0;
  float sum = 0.0;
  const size_t n = freq.size();
  for (size_t i = 0; i != n; i++) {
    sum += freq[i] * i;
    tot += freq[i];
  }

  return sum / tot;
}




size_t min_from_freqs(const vec<size_t> & freq)
{
  size_t min = 0;
  while (freq[min] == 0) min++;
  return min;
}

size_t max_from_freqs(const vec<size_t> & freq)
{
  size_t max = freq.size() - 1;
  while (freq[max] == 0) max--;
  return max;
}




float effective_qual_from_freqs(const vec<size_t> & q_freq)
{

  float ne = 0.0;    // expected number of errors in the data
  float nq = 0.0;
  for (size_t q = 0; q != q_freq.size(); q++) {
    ne += q_freq[q] * error_probability_from_quality(q);
    nq += q_freq[q];
  }

  return quality_from_error_probability(ne / nq);
}







class QualbStats
{
private: 

  size_t _nc;
  
  // ---- quality frequencies for each cycle (or position in a read)
  vec< vec<size_t> > _q_cycle_freqs;

 // ---- read and base quality frequencies 
  vec<size_t> _q_read_freqs;
  vec<size_t> _q_base_freqs;

public:
  QualbStats() 
    : _nc(0),
      _q_cycle_freqs(),
      _q_read_freqs(NQ, 0),
      _q_base_freqs(NQ, 0)
  {
  }

  QualbStats operator += (const QualbStats & qs) 
  {
    for (unsigned q = 0; q != NQ; q++) {
      _q_read_freqs[q] += qs._q_read_freqs[q];
      _q_base_freqs[q] += qs._q_base_freqs[q];
    }
    
    if (qs._nc > _nc) {
      _nc = qs._nc;
      _q_cycle_freqs.resize(_nc, vec<size_t>(NQ, 0));
    }
    for (size_t c = 0; c != qs._nc; c++) {
      for (unsigned q = 0; q != NQ; q++) {
        _q_cycle_freqs[c][q] +=  qs._q_cycle_freqs[c][q];
      }
    }

    return *this;
  }

  size_t number_of_cycles() { return _q_cycle_freqs.size(); }


  void add_qual_vec(const QualVec & qv)
  {
    const size_t nc = qv.size();
    if (nc > _nc) {
      _nc = nc;
      _q_cycle_freqs.resize(_nc, vec<size_t>(NQ, 0));
    }

    float ne = 0.0;  // expected number of errors in read
    for (size_t c = 0; c != nc; c++) {
      const unsigned q = qv[c];
      
      _q_cycle_freqs[c][q]++;
      //_q_base_freqs[q]++;
      
      ne += error_probability_from_quality(q);
    }
    const unsigned q_read = 0.5 + quality_from_error_probability(ne / float(nc));
    _q_read_freqs[q_read]++;
  }





  void summarize()
  {
    for (unsigned c = 0; c != _nc; c++) 
      for (unsigned q = 0; q != NQ; q++)
        _q_base_freqs[q] += _q_cycle_freqs[c][q];
  }



  float    q_eff() { return effective_qual_from_freqs(_q_base_freqs); }
  float    q_mean() { return mean_from_freqs(_q_base_freqs); }
  float    q_min() { return min_from_freqs(_q_base_freqs); }
  float    q_max() { return max_from_freqs(_q_base_freqs); }
  unsigned q_Q1()  { return quantile_from_freqs(_q_base_freqs, 0.25); }
  unsigned q_Q2()  { return quantile_from_freqs(_q_base_freqs, 0.50); }
  unsigned q_Q3()  { return quantile_from_freqs(_q_base_freqs, 0.75); }

  void global_stats_to_ostream(ostream & os)
  {
    os << "# 1:limits 2:q_eff 3:q_mean 4:q_Q1 5:q_Q2 6:q_Q3" << endl;
    os << fixed;
    
    for (unsigned i = 0; i != 2; i++) {
      os << " " << setw(5) << (i * (_nc - 1))  
         << " " << setw(8) << setprecision(2) << q_eff()
         << " " << setw(8) << setprecision(2) << q_mean()
         << " " << setw(5) << q_Q1()
         << " " << setw(5) << q_Q2()
         << " " << setw(5) << q_Q3()
         << endl;
    }
  }


  void cycle_qual_freq_to_ostream(ostream & os, const bool RC)
  {
    for (unsigned c = 0; c != _nc; c++) {
      os << "# 1:cycle 2:q 3:freq" << endl;
      const unsigned cc = ((RC) ? _nc - c - 1 : c);
      for (unsigned q = 1; q != NQ; q++) {
	const size_t f = _q_cycle_freqs[cc][q];
	if (f != 0) {
	  os << " " << setw(5) << c
	     << " " << setw(5) << q
	     << " " << setw(8) << f
	     << endl;
	}
      }
      os << endl;
    }
  }



  void cycle_quals_to_ostream(ostream & os, const bool RC)
  {
    os << "# 1:cycle 2:q_eff 3:q_mean 4:q_Q1 5:q_Q2 6:q_Q3" << endl;
    os << fixed;

    for (unsigned c = 0; c != _nc; c++) {

      const unsigned cc = ((RC) ? _nc - c - 1 : c);
     
      const float q_cycle_eff = effective_qual_from_freqs(_q_cycle_freqs[cc]);
      const float q_cycle_mean = mean_from_freqs(_q_cycle_freqs[cc]);

      // ---- get quantiles
      const unsigned qp25 = quantile_from_freqs(_q_cycle_freqs[cc], 0.25);
      const unsigned qp50 = quantile_from_freqs(_q_cycle_freqs[cc], 0.50);
      const unsigned qp75 = quantile_from_freqs(_q_cycle_freqs[cc], 0.75);
	
      os << " " << setw(5) << c
         << " " << setw(8) << setprecision(2) << q_cycle_eff
         << " " << setw(8) << setprecision(2) << q_cycle_mean
         << " " << setw(5) << qp25
         << " " << setw(5) << qp50
         << " " << setw(5) << qp75
         << endl;
    }
  }


  void qual_freqs_to_ostream(ostream & os) 
  {
    os << "# 1:q 2:read_freq 3:read_freq_cum 4:bases_freq 5:bases_freq_cum" << endl;
    os << scientific;
      
    size_t qfr_tot = 0;
    size_t qfb_tot = 0;
    for (unsigned q = 0; q != NQ; q++) {
      qfr_tot += _q_read_freqs[q];
      qfb_tot += _q_base_freqs[q];
    }
      
    double qfr_cum = 0;
    double qfb_cum = 0;
    double dqfr = 1.0 / double(qfr_tot);
    double dqfb = 1.0 / double(qfb_tot);

    for (unsigned q = 0; q != NQ; q++) {
      const size_t qfr = _q_read_freqs[q];
      const size_t qfb = _q_base_freqs[q];
      if (qfr != 0 || qfb != 0) {

	qfr_cum += qfr * dqfr;
	qfb_cum += qfb * dqfb;
	
        os << " " << setw(5) << q
           << " " << setw(15) << setprecision(6) << qfr
           << " " << setw(15) << setprecision(6) << qfr_cum
           << " " << setw(15) << setprecision(6) << qfb
           << " " << setw(15) << setprecision(6) << qfb_cum
           << endl;
      }
    }
  }

      
};    











class QualbStatsParser
{
  bool _sampled;
  const VQualVecVec & _vqvv;
  const size_t _n_qv;
  const size_t _min_read_len;
  const size_t _max_read_len;

  size_t _n_bases;
  size_t _n_qv_empty;
  size_t _n_qv_short;
  size_t _n_qv_long;
  size_t _n_qv_good;

  vec<size_t> _lens;
  
  QualbStats _q_stats[2];

public:
 
  QualbStatsParser(const VQualVecVec & vqvv, 
		   const size_t min_read_len,
		   const size_t max_read_len)
    : _sampled(true), 
      _vqvv(vqvv),
      _n_qv(_vqvv.size()),
      _min_read_len(min_read_len),
      _max_read_len(max_read_len),
      _n_bases(0),
      _n_qv_empty(0),
      _n_qv_short(0),
      _n_qv_long(0),
      _n_qv_good(0)
  {
    if (_n_qv == 0)
      cout << " total objects count: 0" << endl << endl;
  }

  void scale_stats()
  {
    const size_t n_qv_tot = _n_qv_empty + _n_qv_short + _n_qv_long + _n_qv_good;
    
    _n_bases    = 0.5 + (double(_n_bases)    / double(n_qv_tot)) * _n_qv;
    _n_qv_empty = 0.5 + (double(_n_qv_empty) / double(n_qv_tot)) * _n_qv;
    _n_qv_short = 0.5 + (double(_n_qv_short) / double(n_qv_tot)) * _n_qv;
    _n_qv_long  = 0.5 + (double(_n_qv_long)  / double(n_qv_tot)) * _n_qv;
    _n_qv_good  = 0.5 + (double(_n_qv_good)  / double(n_qv_tot)) * _n_qv;
  }


  void print_lengths_stats()
  {
    sort(_lens.begin(), _lens.end());

    cout << " overall statistics:" << endl << endl;
    cout << fixed;
    cout << setw(13) << _n_qv << "  reads" << endl << endl;

    if (_n_qv_good > 0) {

      cout << ((_sampled) ? 
	       "   statistics estimated from data set sample:" :
	       "   statistics from full data set:")
	 << endl << endl;
      
      cout << setw(13) << _n_bases << "  total length" << endl;
      cout << setw(13) << _lens[0] << "  minimum length" << endl;
      cout << setw(13) << _lens.back() << "  maximum length" << endl;
      cout << setw(13) << setprecision(1) << double(_n_bases) / double(_n_qv_good)
	   << "  mean length" << endl;
      cout << setw(13) << N50(_lens) << "  N50 length" << endl;
    }
    else {
      cout << " no reads found for criteria." << endl << endl;
    }
      
    if (_min_read_len > 0) 
      cout << setw(13) << _n_qv_short 
	   << "  (" << setprecision(1) << 100.0 * double(_n_qv_short)/double(_n_qv) << " %)" 
	   << " reads shorter than " << _min_read_len
	   << endl;
    
    if (_max_read_len < INT_MAX)  
      cout << setw(13) << _n_qv_long 
	   << "  (" << setprecision(1) << 100.0 * double(_n_qv_long)/double(_n_qv) << " %)" 
	   << " reads longer than " << _max_read_len
	   << endl;
    
    if (_n_qv_empty > 0) 
      cout << " Warning! " << _n_qv_empty << " objects had zero length." << endl << endl;


  }


  void print_qual_stats()
  {
    QualbStats q_stats_tot = _q_stats[0];
    q_stats_tot += _q_stats[1];

    cout << endl;
    if (q_stats_tot.number_of_cycles() == 0) {
      cout << endl << " no quality data found for criteria." << endl << endl;
    }
    else {
      cout << " overall quality statistics:" << endl << endl;
      
      cout << setw(13) << setprecision(1) <<  q_stats_tot.q_mean() << "  mean Q" << endl;
      cout << setw(13) << setprecision(1) <<  q_stats_tot.q_eff() << "  effective Q" << endl;
      cout << setw(13) << q_stats_tot.q_Q1() << "  Q1 Q (1st quartile)" << endl;
      cout << setw(13) << q_stats_tot.q_Q2() << "  Q2 Q (median)" << endl;
      cout << setw(13) << q_stats_tot.q_Q3() << "  Q3 Q (3rd quartile)" << endl;
      cout << setw(13) << q_stats_tot.q_min() << "  min Q" << endl;
      cout << setw(13) << q_stats_tot.q_max() << "  max Q" << endl;
      
      if (_q_stats[0].number_of_cycles() == 0) {
	cout << endl << " no first pair found." << endl << endl;
      }
      else if (_q_stats[1].number_of_cycles() == 0) {
	cout << endl << " no second pair found." << endl << endl;
      }
      else {
	cout << endl;
	cout << " pair quality statistics:" << endl << endl;
	cout << setw(6) << setprecision(1) << _q_stats[0].q_mean() << " "
	     << setw(6) << setprecision(1) << _q_stats[1].q_mean() << "  mean Q" << endl;
	
	cout << setw(6) << setprecision(1) << _q_stats[0].q_eff() << " "
	     << setw(6) << setprecision(1) << _q_stats[1].q_eff() << "  effective Q" << endl;
	
	cout << setw(6) << _q_stats[0].q_Q1() << " " 
	     << setw(6) << _q_stats[1].q_Q1() << "  Q1 Q (1st quartile)" << endl;
	cout << setw(6) << _q_stats[0].q_Q2() << " "
	     << setw(6) << _q_stats[1].q_Q2() << "  Q2 Q (median)" << endl;
	cout << setw(6) << _q_stats[0].q_Q3() << " " 
	     << setw(6) << _q_stats[1].q_Q3() << "  Q3 Q (3rd quartile)" << endl;
	cout << setw(6) << _q_stats[0].q_min() << " "
	     << setw(6) << _q_stats[1].q_min() << "  min Q" << endl;
	cout << setw(6) << _q_stats[0].q_max() << " " 
	     << setw(6) << _q_stats[1].q_max() << "  max Q" << endl;
      }
    }
    cout << endl;
  }


  void output_qual_spec(const bool rc,
			const String & out_head)
  {
  
    // ---- output to file
    if (out_head != "") {
      {
	const String stats_fn = out_head + ".stats";
	cout << Date() << ": writing stats to " << stats_fn << endl;
    
	ofstream os;
	os.open(stats_fn.c_str());

	//os << "# stats for '" << HEAD << "'" << endl << endl;
	os << "# n_reads = " << _n_qv << endl;
	os << "# n_bases = " << _n_bases << endl;
	os << "# min_read_len = " << _lens[0] << endl;
	os << "# max_read_len = " << _lens.back() << endl;
	os << endl;


	os << "# 1st reads: global stats" << endl;
	_q_stats[0].global_stats_to_ostream(os);
	os << endl << endl;
	os << "# 2nd reads: global stats" << endl;
	_q_stats[1].global_stats_to_ostream(os);
	os << endl << endl;

	os << "# 1st reads: per cycle stats" << endl;
	_q_stats[0].cycle_quals_to_ostream(os, rc);
	os << endl << endl;
	os << "# 2nd reads: per cycle stats" << endl;
	_q_stats[1].cycle_quals_to_ostream(os, rc);
	os << endl << endl;

	os << "# 1st reads: quality frequencies" << endl;
	_q_stats[0].qual_freqs_to_ostream(os);
	os << endl << endl;
	os << "# 2nd reads: quality frequencies" << endl;
	_q_stats[1].qual_freqs_to_ostream(os);
	os << endl << endl;

	os.close();
      }
      {
	const String cq_freqs_fn = out_head + ".cq_freqs";
	cout << Date() << ": writing cycle quals frequencies to " << cq_freqs_fn << endl;
	ofstream os;
	os.open(cq_freqs_fn.c_str());
	os << "# 1st reads:" << endl;
	_q_stats[0].cycle_qual_freq_to_ostream(os, rc);
	os << endl;
	os << "# 2nd reads:" << endl;
	_q_stats[1].cycle_qual_freq_to_ostream(os, rc);
	os << endl;
	os.close();
      }

    }
    
  }



  
  void add_qual_vec(const QualVec & qv, 
		    QualbStats & qs) 
  {
    const size_t len = qv.size();
    
    if      (len == 0)            _n_qv_empty++;
    else if (len < _min_read_len) _n_qv_short++;
    else if (len > _max_read_len) _n_qv_long++;
    else {
      _n_qv_good++;
      _n_bases += len;

      _lens.push_back(len);
      
      qs.add_qual_vec(qv);	
    }
  }
  

  void sample_indices(const size_t n_blocks,
		      const size_t n_qv_per_block,
		      vec<size_t> & indices)
  {
    const size_t n_indices = n_blocks * n_qv_per_block;
    indices.reserve(n_indices);

    const size_t block_size = _n_qv / n_blocks;
    const size_t d_ind = (block_size - n_qv_per_block)/2;

    for (size_t i = 0; i != n_blocks; i++)
      for (size_t j = 0; j != n_qv_per_block; j++) 
	indices.push_back(i * block_size + d_ind + j);
  }



  void sample_stats(const bool rc,
		    const String & out_head)
  {
    _sampled = true;

    if (_n_qv > 0) {
      const size_t n_blocks = 100;
      const size_t n_reads_per_block = 10000;
      const size_t n_indices = n_blocks * n_reads_per_block;

      if (n_indices > 0.8 * _n_qv) {
	complete_stats(rc, out_head);
      }
      else {

	vec<size_t> indices;
	sample_indices(n_blocks, n_reads_per_block, indices);


	// ---- parse
	cout << " streaming data:" << flush;

	_lens.reserve(n_indices);
	bool do_spec = (out_head != "");

	for (size_t i_ind = 0; i_ind != n_indices; i_ind++) {

	  const size_t i_qv = indices[i_ind];
	  add_qual_vec(_vqvv[i_qv], _q_stats[i_qv & 1]);

	  dots_pct(i_ind, n_indices);
	}
	cout << endl << endl;

	scale_stats();

	_q_stats[0].summarize();
	_q_stats[1].summarize();

	// ---- print

	print_lengths_stats();
	print_qual_stats();
	output_qual_spec(rc, out_head);
	cout << endl;
      }
    }
  }


  void complete_stats(const bool rc,
		      const String & out_head)
  {
    _sampled = false;

    if (_n_qv > 0) {
      // ---- parse
      cout << "Streaming data:" << endl;

      _lens.reserve(_n_qv);
      size_t i_qv = 0;
      for (VQualVecVec::const_iterator it = _vqvv.begin(); 
	   it != _vqvv.end(); i_qv++, it++) {

	add_qual_vec(*it, _q_stats[i_qv & 1]);

	dots_pct(i_qv, _n_qv);
      }
      cout << endl << endl;
 

      _q_stats[0].summarize();
      _q_stats[1].summarize();

      // ---- print
     
      print_lengths_stats();
      print_qual_stats();
      output_qual_spec(rc, out_head);
      cout << endl;
    }
  }



};






















int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String_Doc
    (HEAD, "Input is in <HEAD>.qualb");
  CommandArgument_Int_OrDefault_Doc
    (MIN_READ_LEN, 0, "Only consider reads of length >= than MIN_READ_LEN");
  CommandArgument_Bool_OrDefault_Doc
    (SAMPLE, False, "Estimate stats by sample a subset of the data");
  CommandArgument_Int_OrDefault_Doc
    (MAX_READ_LEN, INT_MAX, "Only consider reads of length <= than MAX_READ_LEN");
  CommandArgument_String_OrDefault_Doc
    (OUT_HEAD, "", "if specified, stats output goes to <OUT_HEAD>.{stats,cq_freq}.");
  CommandArgument_Bool_OrDefault_Doc
    (RC, False, "True if the reads have been RC'ed."); 
  
  EndCommandArguments;



  // ---- reading quals
  const String qualb_fn = HEAD + ".qualb";
  if (!IsRegularFile(qualb_fn)) FatalErr("Could not find file: " + qualb_fn);
    
  // Use VirtualMasterVec to avoid loading in the full QualVecVec
  VQualVecVec vqvv(qualb_fn.c_str());


  QualbStatsParser qsp(vqvv, MIN_READ_LEN, MAX_READ_LEN);

  if (SAMPLE) qsp.sample_stats(RC, OUT_HEAD);
  else        qsp.complete_stats(RC, OUT_HEAD);

}

