//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "Read.h"
#include "bithash.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <set>
#include <queue>

#define TESTING false

int bithash::k;

////////////////////////////////////////////////////////////
// corrections_compare
//
// Simple class to compare to corrected_read's in the
// priority queue
////////////////////////////////////////////////////////////
class corrections_compare {
public:
  //  corrections_compare() {};
  bool operator() (const corrected_read* lhs, const corrected_read* rhs) const {
    //return lhs->likelihood < rhs->likelihood;
    if(lhs->likelihood < rhs->likelihood)
      return true;
    else if(lhs->likelihood > rhs->likelihood)
      return false;
    else
      return lhs->region_edits > rhs->region_edits;
  }
};

const float Read::trust_spread_t = .1;
const float Read::correct_min_t = .000001;
const float Read::learning_min_t = .00001;

////////////////////////////////////////////////////////////
// Read (constructor)
//
// Make shallow copies of sequence and untrusted, and
// convert quality value string to array of probabilities
////////////////////////////////////////////////////////////
Read::Read(const string & h, const unsigned int* s, const string & q, vector<int> & u, const int rl)
  :untrusted(u) {

  header = h;
  read_length = rl;
  trim_length = rl;
  seq = new unsigned int[read_length];
  quals = new unsigned int[read_length];
  prob = new float[read_length];
  for(int i = 0; i < read_length; i++) {
    seq[i] = s[i];
    // quality values of 0,1 lead to p < .25
    quals[i] = q[i] - quality_scale;
    if(quals[i] >= max_qual) {
     cerr << "Quality value " << quals[i] << "larger than maximum allowed quality value " << max_qual << ". Increase the variable 'max_qual' in Read.h." << endl;
     exit(EXIT_FAILURE);
    }     
    prob[i] = max(.25, 1.0-pow(10.0,-(quals[i]/10.0)));
  }
  trusted_read = 0;
  global_like = 1.0;
}

Read::~Read() {
  delete[] seq;
  delete[] quals;
  delete[] prob;
  if(trusted_read != 0)
    delete trusted_read;
}

////////////////////////////////////////////////////////////
// trim
//
// Trim the end of the read the way BWA does it.
// Removes affected untrusted k-mers.
// Returns the trimmed read as a string.
////////////////////////////////////////////////////////////
string Read::trim(int t) {
  // find trim index
  int phredq;
  int current_trimfunc = 0;
  int max_trimfunc = 0;
  trim_length = read_length; // already set in constructor but ok
  for(int i = read_length-1; i >= 0; i--) {
    phredq = floor(.5-10*log(1.0 - prob[i])/log(10));
    current_trimfunc += (t - phredq);
    if(current_trimfunc > max_trimfunc) {
      max_trimfunc = current_trimfunc;
      trim_length = i;
    }
  }

  // update untrusted
  for(int i = untrusted.size()-1; i >= 0; i--) {
    if(untrusted[i] > trim_length - bithash::k)
      untrusted.pop_back();
  }

  vector<correction> no_cors;
  return print_corrected(no_cors);
}


////////////////////////////////////////////////////////////
// single_correct
//
// Find the set of corrections with maximum likelihood
// that result in all trusted kmers.
//
// Assumes a short read so obsolete.
////////////////////////////////////////////////////////////
/*
bool Read::single_correct(bithash *trusted, ofstream & out, double (&ntnt_prob)[4][4], bool learning) {
  if(correct_subset(untrusted, trusted, out, ntnt_prob, learning)) {
    out << header << "\t" << print_seq() << "\t" << print_corrected(trusted_read->corrections) << endl;
    return true;
  } else
    return false;
}
*/

////////////////////////////////////////////////////////////
// correct_cc
//
// Find the set of corrections with maximum likelihood that
// result in all trusted kmers in the region defined by
// the given untrusted kmer indices.
//
// Will print output if the read will not be corrected,
// but otherwise abstains.
//
// Corrections can be accessed through 'trusted_read'
//
// Return codes are:
// 0: corrected
// 1: ambiguous
// 2: full queue or quit early
// 3: empty queue or empty region
////////////////////////////////////////////////////////////
//bool Read::correct_cc(vector<short> region, vector<int> untrusted_subset, bithash *trusted, double (&ntnt_prob)[4][4], double prior_prob[4], bool learning) {  
int Read::correct_cc(vector<short> region, vector<int> untrusted_subset, bithash *trusted, double ntnt_prob[][4][4], double prior_prob[4], bool learning) {  

  unsigned int max_queue_size = 400000;

  /*
  if(header == "@read3") {
    cout << "Untrusted: " << untrusted_subset.size() << endl;
    for(int i = 0; i < untrusted_subset.size(); i++)
      cout << untrusted_subset[i] << " ";
    cout << endl;
    cout << "Region: " << region.size() << endl;
    for(int i = 0; i < region.size(); i++)
      cout << region[i] << " ";
    cout << endl << endl;
  } 
  */

  ////////////////////////////////////////
  // region
  ////////////////////////////////////////
  // sort by quality
  if(region.size() > 0)
    quality_quicksort(region, 0, region.size()-1);
  else
    // die quietly and try again with bigger region
    return 3;

  ////////////////////////////////////////
  // stats
  ////////////////////////////////////////
  unsigned int cpq_adds = 0;
  unsigned int check_count = 0;
  float exp_errors = 0;
  int nt90 = 0;
  int nt99 = 0;
  int non_acgt = 0;
  for(int i = 0; i < region.size(); i++) {
    exp_errors += (1-prob[region[i]]);
    if(prob[region[i]] < .9)
      nt90++;
    if(prob[region[i]] < .99)
      nt99++;
    if(seq[region[i]] >= 4)
      non_acgt++;
  }  

  ////////////////////////////////////////
  // filter
  ////////////////////////////////////////
  double mylike_t = correct_min_t;
  double myglobal_t = correct_min_t;
  double myspread_t = trust_spread_t;
  if(learning) {
    if(nt99 >= 8 || non_acgt > 1) {
      //out << header << "\t" << print_seq() << "\t." << endl;
      return 2;
    }
    mylike_t = learning_min_t;
    myglobal_t = learning_min_t;
    myspread_t = trust_spread_t / 2.0;
  } else if(nt99 >= 13 || non_acgt > 2) {
    // just quit
    if(TESTING)
    cerr << header << "\t" << print_seq() << "\t." << endl;
    //cerr << header << "\t" << region.size() << "\t" << untrusted_subset.size() << "\t" << nt90 << "\t" << nt99 << "\t" << exp_errors << "\t0\t0\t0\t0" << endl;
    return 2;

  } else if(nt99 >= 11) {  
    // proceed very cautiously
    if(aggressive)
      mylike_t = .05;
    else
      mylike_t = .1;

  } else if(nt99 >= 9) {
    //proceed cautiously
    if(aggressive)
      mylike_t = .001;
    else
      mylike_t = .03;
  }

  ////////////////////////////////////////
  // priority queue
  ////////////////////////////////////////
  // data structure for corrected_reads sorted by likelihood
  //priority_queue< corrected_read*, vector<corrected_read*>, corrections_compare > cpq;
  vector<corrected_read*> cpq;
  corrections_compare cpq_comp;

  ////////////////////////////////////////
  // initialize
  ////////////////////////////////////////
  corrected_read *cr, *next_cr;
  short edit_i;
  float like;
  bitset<bitsize> bituntrusted;
  for(int i = 0; i < untrusted_subset.size(); i++) {
       if(untrusted_subset[i] >= bitsize) {
        cerr << "These reads must be longer than assumed. Increase the variable 'bitsize' in 'Read.h' to the read length." << endl;
        exit(1);
       } else
        bituntrusted.set(untrusted_subset[i]);
  }

  bool cr_added = true;  // once an iteration passes w/ no corrected reads added, we can stop
  for(short region_edit = 0; region_edit < region.size() && cr_added; region_edit++) {
    edit_i = region[region_edit];
    cr_added = false;

    for(short nt = 0; nt < 4; nt++) {
      if(seq[edit_i] != nt) {
    // P(obs=o|actual=a)*P(actual=a) for Bayes
    if(seq[edit_i] < 4)
      like = (1.0-prob[edit_i]) * ntnt_prob[quals[edit_i]][nt][seq[edit_i]] * prior_prob[nt] / (prob[edit_i] * prior_prob[seq[edit_i]]);
    else
      // non-ACGT
      like = prior_prob[nt] / (1.0/3.0);

    // P(actual=a|obs=o)
    //like = (1.0-prob[edit_i]) * ntnt_prob[seq[edit_i]][nt] * / prob[edit_i];
    
    next_cr = new corrected_read(bituntrusted, like, region_edit+1);
    next_cr->corrections.push_back(correction(edit_i, nt));
      
    // add to priority queue
    //cpq.push(next_cr);
    cpq.push_back(next_cr);
    push_heap(cpq.begin(), cpq.end(), cpq_comp);
    cpq_adds++;
    cr_added = true;
      }
    }
  }

  ////////////////////////////////////////
  // process corrected reads
  ////////////////////////////////////////
  // initialize likelihood parameters
  trusted_read = 0;
  float trusted_likelihood;
  signed int untrusted_count;  // trust me
  bool ambiguous_flag = false;

  while(cpq.size() > 0) {    

    /////////////////////////
    // quit if pq is too big
    /////////////////////////
    if(cpq.size() > max_queue_size) {
      //cout << "queue is too large for " << header << endl;
      if(TESTING)
          cerr << header << "\t" << print_seq() << "\t." << endl;

      if(trusted_read != 0) {
    delete trusted_read;
    trusted_read = 0;
      }
      break;
    }

    /////////////////////////
    // pop next
    /////////////////////////
    cr = cpq[0];
    pop_heap(cpq.begin(), cpq.end(), cpq_comp);
    cpq.pop_back();
    
    /////////////////////////
    // check likelihood
    /////////////////////////
    if(trusted_read != 0) {
      // if a corrected read exists, compare likelihoods and if likelihood is too low, break loop return true
      if(cr->likelihood < trusted_likelihood*myspread_t) {
    delete cr;
    break;
      }
    } else {
      // if no corrected read exists and likelihood is too low, break loop return false
      if(cr->likelihood < mylike_t || global_like*cr->likelihood < myglobal_t) {
    delete cr;
    break;
      }
    }
    
    /////////////////////////
    // check trust
    /////////////////////////
    // save for later comparison
    untrusted_count = (signed int)cr->untrusted.count();    
    if(check_trust(cr, trusted, check_count)) {
      if(trusted_read == 0) {
    // if yes, and first trusted read, save
    trusted_read = cr;
    trusted_likelihood = cr->likelihood;
      } else {
    // if yes, and if trusted read exists
    ambiguous_flag = true;

    // output ambiguous corrections for testing
    if(TESTING)
      cerr << header << "\t" << print_seq() << "\t" << print_corrected(trusted_read->corrections) << "\t" << print_corrected(cr->corrections) << endl;
    
    // delete trusted_read, break loop
    delete trusted_read;
    delete cr;
    trusted_read = 0;
    break;
      }
    }

    /*
    if(header == "@read3") {
      cout << cr->likelihood << "\t";
      for(int c = 0; c < cr->corrections.size(); c++) {
    cout << " (" << cr->corrections[c].index << "," << cr->corrections[c].to << ")";
      }
      cout << "\t";
      for(int c = 0; c < trim_length-bithash::k+1; c++) {
    if(cr->untrusted[c])
      cout << 1;
    else
      cout << 0;
      }
      cout << endl;
    }
    */
    
    // if untrusted sharply increases, just bail
    if(((signed int)cr->untrusted.count() - untrusted_count)*3 < bithash::k) {

      /////////////////////////
      // add next correction
      /////////////////////////
      bool cr_added = true;  // once an iteration passes w/ no corrected reads added, we can stop
      for(short region_edit = cr->region_edits; region_edit < region.size() && cr_added; region_edit++) {
    edit_i = region[region_edit];
    cr_added = false;
    
    // add relatives
    for(short nt = 0; nt < 4; nt++) {
      // if actual edit, 
      if(seq[edit_i] != nt) {
        // calculate new likelihood

        // P(obs=o|actual=a)*P(actual=a) for Bayes
        if(seq[edit_i] < 4)
          like = cr->likelihood * (1.0-prob[edit_i]) * ntnt_prob[quals[edit_i]][nt][seq[edit_i]] * prior_prob[nt] / (prob[edit_i] * prior_prob[seq[edit_i]]);
        else
          // non-ACGT
          like = cr->likelihood * prior_prob[nt] / (1.0/3.0);

        // P(actual=a|obs=o)        
        //like = cr->likelihood * (1.0-prob[edit_i]) * ntnt_prob[seq[edit_i]][nt] / prob[edit_i];    
        
        // if thresholds ok, add new correction
        if(trusted_read != 0) {
          if(like < trusted_likelihood*myspread_t)
        continue;
        } else {
          // must consider spread or risk missing a case of ambiguity
          if(like < mylike_t*myspread_t || global_like*like < myglobal_t*myspread_t)
        continue;
        }
        
        next_cr = new corrected_read(cr->corrections, cr->untrusted, like, region_edit+1);
        next_cr->corrections.push_back(correction(edit_i, nt));
      
        // add to priority queue
        cpq.push_back(next_cr);
        push_heap(cpq.begin(), cpq.end(), cpq_comp);
        cpq_adds++;
        cr_added = true;
      }
    }
      }
    }

    // if not the saved max trusted, delete
    if(trusted_read != cr) {
      delete cr;
    }
  }

  // clean up priority queue
  for(int i = 0; i < cpq.size(); i++)
    delete cpq[i];
  
  if(trusted_read != 0) {
    //cerr << header << "\t" << region.size() << "\t" << untrusted_subset.size() << "\t" << nt90 << "\t" << nt99 << "\t" << exp_errors << "\t" << cpq_adds << "\t" << check_count << "\t1\t" << trusted_read->likelihood << endl;
    return 0;
  } else {
    if(TESTING && mylike_t > correct_min_t)
      cerr << header << "\t" << print_seq() << "\t." << endl;
    //cerr << header << "\t" << region.size() << "\t" << untrusted_subset.size() << "\t" << nt90 << "\t" << nt99 << "\t" << exp_errors << "\t" << cpq_adds << "\t" << check_count << "\t0\t0" << endl;

    if(ambiguous_flag)
      return 1;
    else if(cpq.size() > max_queue_size)
      return 2;
    else
      return 3;
  }
}

////////////////////////////////////////////////////////////
// print_seq
////////////////////////////////////////////////////////////
string Read::print_seq() {
  char nts[5] = {'A','C','G','T','N'};
  string sseq;
  for(int i = 0; i < read_length; i++)
    sseq.push_back(nts[seq[i]]);
  return sseq;
}

////////////////////////////////////////////////////////////
// print_corrected
//
// Print read with corrections and trimming.
////////////////////////////////////////////////////////////
string Read::print_corrected(vector<correction> & cor) {
  return print_corrected(cor, trim_length);
}
string Read::print_corrected(vector<correction> & cor, int print_nt) {
  char nts[5] = {'A','C','G','T','N'};
  string sseq;
  int correct_i;
  for(int i = 0; i < print_nt; i++) {
    correct_i = -1;
    for(int c = 0; c < cor.size(); c++) {
      if(cor[c].index == i)
    correct_i = c;
    }
    if(correct_i != -1)
      sseq.push_back(nts[cor[correct_i].to]);
    else
      sseq.push_back(nts[seq[i]]);
  }
  return sseq;
}
  

////////////////////////////////////////////////////////////
// correct
//
// Perform correction by breaking up untrusted kmers
// into connected components and correcting them
// independently.
////////////////////////////////////////////////////////////
//string Read::correct(bithash *trusted, double (&ntnt_prob)[4][4], double prior_prob[4], bool learning) {
string Read::correct(bithash *trusted, double ntnt_prob[][4][4], double prior_prob[4], bool learning) {
  ////////////////////////////////////////
  // find connected components
  ////////////////////////////////////////
  vector< vector<int> > cc_untrusted;

  // add first
  cc_untrusted.push_back(vector<int>());
  int cc = 0;
  cc_untrusted[cc].push_back(untrusted[0]);

  for(int i = 1; i < untrusted.size(); i++) {
    // if kmer from last untrusted doesn't reach next
    if(untrusted[i-1]+bithash::k-1 < untrusted[i]) {
      cc++;
      cc_untrusted.push_back(vector<int>());
    }
    cc_untrusted[cc].push_back(untrusted[i]);
  }

  ////////////////////////////////////////
  // process connected components
  ////////////////////////////////////////
  vector<correction> multi_cors;
  vector<short> chop_region;
  vector<short> big_region;
  int chop_correct_code, big_correct_code;
  for(cc = 0; cc < cc_untrusted.size(); cc++) {
    // try chopped error region
    chop_region = error_region_chop(cc_untrusted[cc]);
    chop_correct_code = correct_cc(chop_region, cc_untrusted[cc], trusted, ntnt_prob, prior_prob, learning);
    if(chop_correct_code > 0) {
      // try bigger error region
      big_region = error_region(cc_untrusted[cc]);
      if(chop_region.size() == big_region.size()) {
    // cannot correct, and nothing found so trim to untrusted
    if(chop_correct_code == 1)
      return print_corrected(multi_cors, chop_region.front());
    else
      return print_corrected(multi_cors, cc_untrusted[cc].front());

      } else {
    big_correct_code = correct_cc(big_region, cc_untrusted[cc], trusted, ntnt_prob, prior_prob, learning);

    if(big_correct_code == 1) {
      // ambiguous
      // cannot correct, but trim to region
      if(chop_correct_code == 1)
        return print_corrected(multi_cors, chop_region.front());
      else
        return print_corrected(multi_cors, big_region.front());

    } else if(big_correct_code == 2 || big_correct_code == 3) {
      // cannot correct, and chaotic or nothing found so trim to untrusted
      return print_corrected(multi_cors, cc_untrusted[cc].front());
    }
      }
    }
    // else, corrected!

    // corrected
    global_like *= trusted_read->likelihood;

    // store
    for(int c = 0; c < trusted_read->corrections.size(); c++)
      multi_cors.push_back(trusted_read->corrections[c]);
  }  

  // create new trusted read (mostly for learn_errors)
  corrected_read * tmp = trusted_read;
  trusted_read = new corrected_read(multi_cors, tmp->untrusted, global_like, 0);
  delete tmp;

  // print read with all corrections
  return print_corrected(multi_cors);
}


////////////////////////////////////////////////////////////
// error_region
//
// Find region of the read to consider for errors based
// on the pattern of untrusted kmers
////////////////////////////////////////////////////////////
vector<short> Read::error_region(vector<int> untrusted_subset) {
  // find intersection, or union
  vector<short> region;
  if(!untrusted_intersect(untrusted_subset, region))
    untrusted_union(untrusted_subset, region);

  // if front kmer can reach region, there may be more
  // errors in the front
  short f = region.front();
  short b = region.back();
  
  if(bithash::k-1 >= f) {
    // extend to front
    for(short i = f-1; i >= 0; i--)
      region.push_back(i);
  }
  if(trim_length-bithash::k <= b) {
    // extend to back
    for(short i = b+1; i < trim_length; i++)
      region.push_back(i);
  }

  return region;
}

////////////////////////////////////////////////////////////
// error_region_chop
//
// Find region of the read to consider for errors based
// on the pattern of untrusted kmers, using trusted kmers
// to further trim the area.
////////////////////////////////////////////////////////////
vector<short> Read::error_region_chop(vector<int> untrusted_subset) {
  // find intersection, or union
  vector<short> region;
  if(!untrusted_intersect(untrusted_subset, region))
    untrusted_union(untrusted_subset, region);

  // fix front
  int right_leftkmer = untrusted_subset.front()-1;
  if(right_leftkmer >= 0) {
    // erase all bp in rightmost left kmer
    vector<short> front_chop(region);
    region.clear();
    for(int i = 0; i < front_chop.size(); i++) {
      if(front_chop[i] > right_leftkmer+bithash::k-1)
    region.push_back(front_chop[i]);
    }

    // add back 1 base if it's low quality, or lower quality than the next base
    for(int er = 0; er < expand_region; er++) {
      int pre_region = region[0] - (er+1);
      if(pre_region >= 0 && (prob[pre_region] < .99 || prob[pre_region] <= prob[pre_region+1])) {
    vector<short>::iterator it;
    it = region.begin();
    region.insert(it, pre_region);
      }
    }
  } else {
    // extend to front
    for(int i = region[0]-1; i >= 0; i--)
      region.push_back(i);
  }

  // fix back
  int left_rightkmer = untrusted_subset.back()+1;
  if(left_rightkmer+bithash::k-1 < trim_length) {
    // erase all bp in leftmost right kmer
    vector<short> back_chop(region);
    region.clear();
    for(int i = 0; i < back_chop.size(); i++) {
      if(back_chop[i] < left_rightkmer)
    region.push_back(back_chop[i]);
    }

    // add back 1 base if it's low quality, or lower quality than the next base
    // Two issues with this:
    // 1. I think region could be empty, so there's a bug
    // 2. This won't help for errors in the middle of a read that are missing an untrusted kmer
    //    because the region will be empty, and we'll just try the intersection.
    /*
    for(int er = 0; er < expand_region; er++) {
      int post_region = region.back() + (er+1);
      if(post_region < trim_length && (prob[post_region] < .99 || prob[post_region] <= prob[post_region-1])) {
    region.push_back(post_region);
      }
    }
    */

  } else {
    // extend to back
    for(int i = region.back()+1; i < trim_length; i++)
      region.push_back(i);
  }

  return region;
}
    
////////////////////////////////////////////////////////////
// untrusted_intersect
//
// Compute the intersection of the untrusted kmers as
// start,end return true if it's non-empty or false
// otherwise
////////////////////////////////////////////////////////////
bool Read::untrusted_intersect(vector<int> untrusted_subset, vector<short> & region) {
  int start = 0;
  int end = read_length-1;

  int u;
  for(int i = 0; i < untrusted_subset.size(); i++) {
    u = untrusted_subset[i];

    // if overlap
    if(start <= u+bithash::k-1 && u <= end) {
      // take intersection
      start = max(start, u);
      end = min(end, u+bithash::k-1);
    } else {
      // intersection is empty
      return false;   
    }
  }
    
  // intersection is non-empty
  for(short i = start; i <= end; i++)
    region.push_back(i);
  return true;
}

////////////////////////////////////////////////////////////
// untrusted_union
//
// Compute the union of the untrusted kmers, though not
////////////////////////////////////////////////////////////
void Read::untrusted_union(vector<int> untrusted_subset, vector<short> & region) {
  short u;
  set<short> region_set;
  for(int i = 0; i < untrusted_subset.size(); i++) {
    u = untrusted_subset[i];

    for(short ui = u; ui < u+bithash::k; ui++)
      region_set.insert(ui);
  }

  set<short>::iterator it;
  for(it = region_set.begin(); it != region_set.end(); it++)
    region.push_back(*it);      
}

////////////////////////////////////////////////////////////
// quality_quicksort
//
// Sort the indexes from lowest probability of an accurate
// basecall to highest
////////////////////////////////////////////////////////////
void Read::quality_quicksort(vector<short> & indexes, int left, int right) {
  int i = left, j = right;  
  short tmp;
  float pivot = prob[indexes[(left + right) / 2]];
 
  /* partition */
  while (i <= j) {
    while (prob[indexes[i]] < pivot)
      i++;    
    while (prob[indexes[j]] > pivot)      
      j--;
    if (i <= j) {
      tmp = indexes[i];
      indexes[i] = indexes[j];
      indexes[j] = tmp;
      i++;
      j--;
    }
  }

  /* recursion */
  if (left < j)
    quality_quicksort(indexes, left, j);
  if (i < right)
    quality_quicksort(indexes, i, right);  
}

////////////////////////////////////////////////////////////
// check_trust
//
// Given a corrected read and data structure holding
// trusted kmers, update the corrected_reads's vector
// of untrusted kmers and return true if it's now empty
////////////////////////////////////////////////////////////
bool Read::check_trust(corrected_read *cr, bithash *trusted, unsigned int & check_count) {
  // original read HAS errors
  if(cr->corrections.empty())
    return false;

  // make corrections to sequence, saving nt's to fix later
  vector<int> seqsave;
  int i;
  for(i = 0; i < cr->corrections.size(); i++) {
    seqsave.push_back(seq[cr->corrections[i].index]);
    seq[cr->corrections[i].index] = cr->corrections[i].to;
  }

  int edit = cr->corrections.back().index;
  int kmer_start = max(0, edit-bithash::k+1);
  //int kmer_end = min(edit, read_length-k);
  int kmer_end = min(edit, trim_length-bithash::k);

  check_count += (kmer_end - kmer_start + 1);

  bool non_acgt = false;
  for(i = kmer_start; i < kmer_end+bithash::k; i++)
    if(seq[i] >=4)
      non_acgt = true;

  //non_acgt = true;
  if(non_acgt) {
    // easier to just check kmers one by one
    for(i = kmer_start; i <= kmer_end; i++)
      // check kmer
      cr->untrusted.set(i, !trusted->check(&seq[i]));

  } else {
    // check affected kmers
    Seq<kK> kmermap;
    // check first kmer and save map value  
    cr->untrusted.set(kmer_start, !trusted->check(&seq[kmer_start], kmermap));
    for(i = kmer_start+1; i <= kmer_end; i++) {
      // check kmer using map value
      cr->untrusted.set(i, !trusted->check(kmermap, seq[i-1], seq[i+bithash::k-1]));
    }
  }

  // fix sequence
  for(i = 0; i < cr->corrections.size(); i++)
    seq[cr->corrections[i].index] = seqsave[i];

  return(cr->untrusted.none());
}
