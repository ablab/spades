//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "bithash.h"
#include "Read.h"
#include "edit.h"
#include "gzstream.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <string.h>
#include <cstring>
#include <getopt.h>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_set_num_threads(x)
#define omp_get_max_threads()    1
#define omp_get_thread_num()     0
#define omp_get_num_threads()    0
#endif
#include <cstdlib>
#include <iomanip>
#include <sys/stat.h>

////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////
const static char* myopts = "r:f:k:m:b:c:a:t:q:l:p:v:zCuh";
static struct option  long_options [] = {
  {"headers", 0, 0, 1000},
  {"log", 0, 0, 1001},
  {0, 0, 0, 0}
};

// -r, fastq file of reads
//char* fastqf;
// -f, file of fastq files of reads
//char* file_of_fastqf;

// -z, zip output files
//bool zip_output = false;

// -k, kmer size
static int k = 0;

// -m, mer counts
static char* merf = NULL;
// -b, bithash
static char* bithashf = NULL;
// -v, good kmers from Hammer
static char* hammerf = NULL;

// -c, cutoff between trusted and untrusted mers
static double cutoff = 0;
// -a, AT cutoff
static char* ATcutf = NULL;

// -q
//int Read::quality_scale;
// -l
static int trim_t = 30;
// -t
//static int trimq = 3;

// -p, number of threads
//int threads;

// --headers, Print only normal headers
static bool orig_headers = false;

// -C, Contrail output
static bool contrail_out = false;
// -u, output uncorrected reads
static bool uncorrected_out = false;
// --log, output correction log
static bool out_log = false;

static bool overwrite_temp = true;

// Note: to not trim, set trimq=0 and trim_t>read_length-k

// constants
#define TESTING false
static char* nts = (char*)"ACGTN";
//unsigned int chunks_per_thread = 200;

 // to collect stats
struct stats {
  stats() {
    validated = 0;
    corrected = 0;
    removed = 0;
    trimmed = 0;
    trimmed_only = 0;
  }
  unsigned long long validated;
  unsigned long long corrected;
  unsigned long long removed;
  unsigned long long trimmed;
  unsigned long long trimmed_only;
};

static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

{
  fprintf (stderr,
           "USAGE:  correct [options]\n"
           "\n"
       "Correct sequencing errors in fastq file provided with -r\n"
       "and output trusted and corrected reads to\n"
       "<fastq-prefix>.cor.fastq.\n"
           "\n"
           "Options:\n"
           " -r <file>\n"
       "    Fastq file of reads\n"
       " -f <file>\n"
       "    File containing fastq file names, one per line or\n"
       "    two per line for paired end reads.\n"
       " -z\n"
       "    Write output files as gzipped.\n"
       " -m <file>\n"
       "    File containing kmer counts in format `seq\tcount`.\n"
       "    Can be gzipped.\n"
       " -b <file>\n"
       "    File containing saved bithash.\n"
       " -c <num>\n"
       "    Separate trusted/untrusted kmers at cutoff <num>\n"
       " -a <file>\n"
       "    Separate trusted/untrusted kmers as a function of AT\n"
       "    content, with cutoffs found in <file>, one per line\n"
       " -p <num>\n"
       "    Use <num> openMP threads\n"
       " -l <num>=30\n"
       "    Return only reads corrected and/or trimmed to >= <num>\n"
       "    bp\n"
       " -q <num>\n"
       "    Quality value ascii scale, generally 64 or 33. If not\n"
       "    specified, it will guess.\n"
       " -v <file>\n"
       "    File with good k-mers from Hammer.\n"
       " -t <num>=3\n"
       "    Use BWA trim parameter <num>\n"
       " -u\n"
       "    Output errors reads even if they can't be corrected,\n"
       "    maintaining paired end reads.\n"
       " --headers\n"
       "    Output only the original read headers without\n"
       "    correction messages\n"
       " --log\n"
       "    Output a log of all corrections into *.log as\n"
       "    'quality position new_nt old_nt'\n"
           "\n");

   return;
  }

////////////////////////////////////////////////////////////
// parse_command_line
////////////////////////////////////////////////////////////
static void parse_command_line(int argc, char **argv) {
  bool errflg = false;
  int ch;
  optarg = NULL;
  int option_index = 0;
  char* p;
  k = kK;
  // parse args
  while(!errflg && ((ch = getopt_long(argc, argv, myopts, long_options, &option_index)) != EOF)) {
  //while(!errflg && ((ch = getopt(argc, argv, myopts)) != EOF)) {
    switch(ch) {
    case 'r':
      fastqf = strdup(optarg);
      break;

    case 'f':
      file_of_fastqf = strdup(optarg);
      break;

    case 'z':
      zip_output = true;
      break;

    case 'm':
      merf = strdup(optarg);
      break;

    case 'b':
      bithashf = strdup(optarg);
      break;

    case 'c':
      cutoff = double(strtod(optarg, &p));
      if(p == optarg || cutoff < 0) {
    fprintf(stderr, "Bad mer cutoff value \"%s\"\n",optarg);
    errflg = true;
      }
      break;

    case 'a':
      ATcutf = strdup(optarg);
      break;

    case 't':
      trimq = int(strtol(optarg, &p, 10));
      if(p == optarg || trimq < 0) {
    fprintf(stderr, "Bad trim quality value \"%s\"\n",optarg);
    errflg = true;
      }
      break;

    case 'l':
      trim_t = int(strtol(optarg, &p, 10));
      if(p == optarg || trim_t < 1) {
    fprintf(stderr, "Bad trim threshold \"%s\"\n",optarg);
    errflg = true;
      }
      break;

    case 'q': 
      Read::quality_scale = int(strtol(optarg, &p, 10));
      if(p == optarg || Read::quality_scale < -1) {
    fprintf(stderr, "Bad quality value scale \"%s\"\n",optarg);
    errflg = true;
      }
      break;

    case 'C':
      contrail_out = true;
      break;

    case 'u':
      uncorrected_out = true;
      break;  

    case 'p':
      threads = int(strtol(optarg, &p, 10));
      if(p == optarg || threads <= 0) {
    fprintf(stderr, "Bad number of threads \"%s\"\n",optarg);
    errflg = true;
      }
      break;

    case 1000:
      orig_headers = true;
      break;

    case 1001:
      out_log = true;
      break;

    case 'v':
      hammerf = strdup(optarg);
      break;
      
    case 'h':
      Usage(argv[0]);
      exit(EXIT_FAILURE);

    case '?':
      fprintf (stderr, "Unrecognized option -%c\n", optopt);

    default:
      errflg = true;
    }
  }

  // for some reason, optind is not advancing properly so this
  // always returns an error

  // return errors
  /*
  if(errflg || optind != argc-1) {
    Usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  */

  ////////////////////////////////////////
  // correct user input errors
  ////////////////////////////////////////
  if(fastqf == NULL && file_of_fastqf == NULL) {
    cerr << "Must provide a fastq file of reads (-r) or a file containing a list of fastq files of reads (-f)" << endl;
    exit(EXIT_FAILURE);
  }

  if(k == 0) {
    cerr << "Must provide kmer size (-k)" << endl;
    exit(EXIT_FAILURE);
  }

  if(merf != NULL) {
    if(cutoff == 0 && ATcutf == NULL) {
      cerr << "Must provide a trusted/untrusted kmer cutoff (-c) or a file containing the cutoff as a function of the AT content (-a)" << endl;
      exit(EXIT_FAILURE);
    }
  } else if(bithashf == NULL && hammerf == NULL) {
    cerr << "Must provide a file of kmer counts (-m) or a saved bithash (-b) or solid kmers from Hammer (-v)" << endl;
    exit(EXIT_FAILURE);
  }
  
}


////////////////////////////////////////////////////////////
// regress_probs
//
// Use ntnt_counts to perform nonparametric regression
// on ntnt_prob across quality values.
////////////////////////////////////////////////////////////
void regress_probs(double ntnt_prob[Read::max_qual][4][4], unsigned int ntnt_counts[Read::max_qual][4][4]) {
  double sigma = 2.0;
  double sigma2 = pow(sigma, 2);
  
  // count # occurrences for each (quality=q,actual=a) tuple
  unsigned int actual_counts[Read::max_qual][4] = {0};
  for(int q = 1; q < Read::max_qual; q++)
    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++)
    actual_counts[q][i] += ntnt_counts[q][i][j];

  // regress
  double ntdsum;
  for(int q = 1; q < Read::max_qual; q++) {
    for(int i = 0; i < 4; i++) {
      //ntdsum = 0;
      for(int j = 0; j < 4; j++) {
    double pnum = 0;
    double pden = 0;
    for(int qr = 1; qr < Read::max_qual; qr++) {
      pnum += ntnt_counts[qr][i][j] * exp(-pow((double)(qr - q), 2)/(2*sigma2));
      pden += actual_counts[qr][i] * exp(-pow((double)(qr - q), 2)/(2*sigma2));
    }
    ntnt_prob[q][i][j] = pnum / pden;
    //ntdsum += ntnt_prob[q][i][j];
      }

      // re-normalize to sum to 1
      //for(int j = 0; j < 4; j++)
      //ntnt_prob[q][i][j] /= ntdsum;
    }
  }
}


////////////////////////////////////////////////////////////
// output_model
//
// Print the error model to the file error_model.txt
////////////////////////////////////////////////////////////
void output_model(double ntnt_prob[Read::max_qual][4][4], unsigned int ntnt_counts[Read::max_qual][4][4], string fqf) {
  string base = split(fqf,'/').back();

  int suffix_index = base.rfind(".");
  string prefix;
  if(suffix_index == -1) {
    prefix = base;
  } else {
    prefix = base.substr(0,suffix_index);
  }

  string outf = "error_model." + prefix + ".txt";

  ofstream mod_out(outf.c_str());

  unsigned int ntsum;
  for(int q = 1; q < Read::max_qual; q++) {
    mod_out << "Quality = " << q << endl;

    // counts
    mod_out << "\tA\tC\tG\tT" << endl;
    for(int i = 0; i < 4; i++) {
      mod_out << nts[i];

      ntsum = 0;
      for(int j = 0; j < 4; j++)
    ntsum += ntnt_counts[q][i][j];

      for(int j = 0; j < 4; j++) {
    if(i == j)
      mod_out << "\t-";
    else if(ntsum > 0)
      mod_out << "\t" << ((double)ntnt_counts[q][i][j] / (double)ntsum) << "(" << ntnt_counts[q][i][j] << ")";
    else
      mod_out << "\t0";
      }
      mod_out << endl;
    }

    // probs
    mod_out << "\tA\tC\tG\tT" << endl;
    for(int i = 0; i < 4; i++) {
      mod_out << nts[i];
      for(int j = 0; j < 4; j++) {
    if(i == j)
      mod_out << "\t-";
    else
      mod_out << "\t" << ntnt_prob[q][i][j];
      }
      mod_out << endl;
    }
    mod_out << endl;    
  }
}


////////////////////////////////////////////////////////////////////////////////
// output_read
//
// Output the given possibly corrected and/or trimmed
// read according to the given options.
////////////////////////////////////////////////////////////////////////////////
static void output_read(ofstream & reads_out, ofstream & corlog_out, int pe_code, string header, string ntseq, string mid, string strqual, string corseq, stats & tstats) {
  if(corseq.size() >= trim_t) {
    // check for changes
    bool corrected = false;
    for(int i = 0; i < corseq.size(); i++) {
      if(corseq[i] != ntseq[i]) {
    // log it
    if(corlog_out.good())
      corlog_out << (strqual[i]-Read::quality_scale) << "\t" << (i+1) << "\t" << corseq[i] << "\t" << ntseq[i] << endl;
    // note it
    corrected = true;
    // set qual to crap
    strqual[i] = (char)(Read::quality_scale+2);
      }
    }
    if(corrected)
      tstats.corrected++;

    // update header
    if(!orig_headers) {
      if(corrected)
    header += " correct";
      unsigned int trimlen = ntseq.size()-corseq.size();
      if(trimlen > 0) {
    stringstream trim_inter;
    trim_inter << trimlen;
    header += " trim=" + trim_inter.str();
    tstats.trimmed++;
    if(!corrected)
      tstats.trimmed_only++;
      } else {
    if(!corrected)
      tstats.validated++;
      }
    }
    // print
    if(contrail_out)
      reads_out << header << "\t" << corseq << endl;
    else
      reads_out << header << endl << corseq << endl << mid << endl << strqual.substr(0,corseq.size()) << endl;
    if(TESTING)
      cerr << header << "\t" << ntseq << "\t" << corseq << endl;
  } else {
    tstats.removed++;
    if(uncorrected_out || pe_code > 0) {
      // update header
      header += " error";

      //print
      if(contrail_out)
    reads_out << header << "\t" << ntseq << endl;
      else
    reads_out << header << endl << ntseq << endl << mid << endl << strqual << endl;      
    }
    if(TESTING)
      cerr << header << "\t" << ntseq << "\t-" << endl; // or . if it's only trimmed?
  }
}


////////////////////////////////////////////////////////////////////////////////
// correct_reads
//
// Correct the reads in the file 'fqf' using the data structure of trusted
// kmers 'trusted', matrix of nt->nt error rates 'ntnt_prob' and prior nt
// probabilities 'prior_prob'.  'starts' and 'counts' help openMP parallelize
// the read processing.  If 'pairedend_code' is 0, the reads are not paired;
// if it's 1, this file is the first of a pair so print all reads and withold
// combining; if it's 2, the file is the second of a pair so print all reads
// and then combine both 1 and 2.
////////////////////////////////////////////////////////////////////////////////
static void correct_reads(string fqf, int pe_code, bithash * trusted, vector<streampos> & starts, vector<unsigned long long> & counts, double ntnt_prob[Read::max_qual][4][4], double prior_prob[4]) {
  // output directory
  struct stat st_file_info;
  string path_suffix = split(fqf,'/').back();
  string out_dir("."+path_suffix);
  if(stat(out_dir.c_str(), &st_file_info) == 0) {
    cerr << "Hidden temporary directory " << out_dir << " already exists and will be used" << endl;
  } else {
    if(mkdir(out_dir.c_str(), S_IRWXU) == -1) {
      cerr << "Failed to create hidden temporary directory " << out_dir << endl;
      exit(EXIT_FAILURE);
    }
  }

  // collect stats
  stats * thread_stats = new stats[omp_get_max_threads()];

  unsigned int chunk = 0;
#pragma omp parallel //shared(trusted)
  {
    int tid = omp_get_thread_num();
    
    // input
    ifstream reads_in(fqf.c_str());
    
    unsigned int tchunk;
    string header,ntseq,mid,strqual,corseq;
    int trim_length;
    char* nti;
    Read *r;

    #pragma omp critical
    tchunk = chunk++;

    while(tchunk < starts.size()) {
      reads_in.seekg(starts[tchunk]);

      // output
      string toutf(out_dir+"/");
      stringstream tconvert;
      tconvert << tchunk;
      toutf += tconvert.str();

      if(overwrite_temp || stat(toutf.c_str(), &st_file_info) == -1) {
    ofstream reads_out(toutf.c_str());
    //cout << toutf << endl;

    // output log
    string tlogf = toutf + ".log";
    ofstream corlog_out;
    if(out_log) {
      corlog_out.open(tlogf.c_str());
    }

    unsigned long long tcount = 0;
    while(getline(reads_in, header)) {
      //cout << tid << " " << header << endl;
    
      // get sequence
      getline(reads_in, ntseq);
      //cout << ntseq << endl;
    
      // convert ntseq to iseq
      vector<unsigned int> iseq;
      for(int i = 0; i < ntseq.size(); i++) {
        nti = strchr(nts, ntseq[i]);    
        iseq.push_back(nti - nts);
      }

      // get quality values
      getline(reads_in,mid);
      //cout << mid << endl;
      getline(reads_in,strqual);
      //cout << strqual << endl;

      vector<int> untrusted;

      if(iseq.size() < trim_t)
        trim_length = 0;
      else {
        for(int i = 0; i < iseq.size()-k+1; i++) {
          if(!trusted->check(&iseq[i])) {
        untrusted.push_back(i);
          }
        }

        trim_length = quick_trim(strqual, untrusted);
        //trim_length = iseq.size();
      }

      // fix error reads
      if(untrusted.size() > 0) {
        r = new Read(header, &iseq[0], strqual, untrusted, trim_length);
        corseq = r->correct(trusted, ntnt_prob, prior_prob);

        // output read w/ trim and corrections
        output_read(reads_out, corlog_out, pe_code, header, ntseq, mid, strqual, corseq, thread_stats[tid]);
      
        delete r;
      } else {
        output_read(reads_out, corlog_out, pe_code, header, ntseq, mid, strqual, ntseq.substr(0,trim_length), thread_stats[tid]);
        // output read as trimmed
        /*
          if(contrail_out)
          reads_out << header << "\t" << ntseq.substr(0,trim_length) << endl;
          else
          reads_out << header << endl << ntseq.substr(0,trim_length) << endl << mid << endl << strqual.substr(0,trim_length) << endl;
        */
      }
    
      if(++tcount == counts[tchunk])
        break;
    }
    reads_out.close();
      }

#pragma omp critical
      tchunk = chunk++;
    }
    reads_in.close();
  }

  // combine stats
  for(int i = 1; i < omp_get_max_threads(); i++) {
    thread_stats[0].validated += thread_stats[i].validated;
    thread_stats[0].corrected += thread_stats[i].corrected;
    thread_stats[0].trimmed += thread_stats[i].trimmed;
    thread_stats[0].trimmed_only += thread_stats[i].trimmed_only;
    thread_stats[0].removed += thread_stats[i].removed;
  }

  // print stats
  int suffix_index = fqf.rfind(".");
  string outf;
  if(suffix_index == -1) {
    outf = fqf+".stats.txt";
  } else {
    outf = fqf.substr(0,suffix_index+1) + "stats.txt";
  }
  ofstream stats_out(outf.c_str());
  stats_out << "Validated: " << thread_stats[0].validated << endl;
  stats_out << "Corrected: " << thread_stats[0].corrected << endl;
  stats_out << "Trimmed: " << thread_stats[0].trimmed << endl;
  stats_out << "Trimmed only: " << thread_stats[0].trimmed_only << endl;
  stats_out << "Removed: " << thread_stats[0].removed << endl;
  stats_out.close();
}


////////////////////////////////////////////////////////////
// learn_errors
//
// Correct reads using a much stricter filter in order
// to count the nt->nt errors and learn the errors
// probabilities
////////////////////////////////////////////////////////////
//static void learn_errors(string fqf, bithash * trusted, vector<streampos> & starts, vector<unsigned long long> & counts, double (&ntnt_prob)[4][4], double prior_prob[4]) {
static void learn_errors(string fqf, bithash * trusted, vector<streampos> & starts, vector<unsigned long long> & counts, double ntnt_prob[Read::max_qual][4][4], double prior_prob[4]) {
  unsigned int ntnt_counts[Read::max_qual][4][4] = {0};
  unsigned int samples = 0;

  unsigned int chunk = 0;
#pragma omp parallel //shared(trusted)
  {    
    unsigned int tchunk;
    string header,ntseq,strqual,corseq;
    int trim_length;
    char* nti;
    Read *r;    
    ifstream reads_in(fqf.c_str());
    
    while(chunk < threads*chunks_per_thread) {
#pragma omp critical
      tchunk = chunk++;     
      
      reads_in.seekg(starts[tchunk]);
      
      unsigned long long tcount = 0;
      while(getline(reads_in, header)) {
    //cout << header << endl;
    
    // get sequence
    getline(reads_in, ntseq);
    //cout << ntseq << endl;
    
    // convert ntseq to iseq
    vector<unsigned int> iseq;
    for(int i = 0; i < ntseq.size(); i++) {
      nti = strchr(nts, ntseq[i]);
      iseq.push_back(nti - nts);
    }
        
    // get quality values
    getline(reads_in,strqual);
    //cout << strqual << endl;
    getline(reads_in,strqual);
    //cout << strqual << endl;

    vector<int> untrusted;

    if(iseq.size() < trim_t)
      trim_length = 0;
    else {
      for(int i = 0; i < iseq.size()-k+1; i++) {
        if(!trusted->check(&iseq[i])) {
          untrusted.push_back(i);
        }
      }
      
      trim_length = quick_trim(strqual, untrusted);
    }

    // fix error reads
    if(untrusted.size() > 0) {
      // correct
      r = new Read(header, &iseq[0], strqual, untrusted, trim_length);
      corseq = r->correct(trusted, ntnt_prob, prior_prob, true);
        
      // if trimmed to long enough
      if(corseq.size() >= trim_t) {
        if(r->trusted_read != 0) { // else no guarantee there was a correction
          for(int c = 0; c < r->trusted_read->corrections.size(); c++) {
        correction cor = r->trusted_read->corrections[c];
        if(iseq[cor.index] < 4) {
          // P(obs=o|actual=a,a!=o) for Bayes
          ntnt_counts[strqual[cor.index]-Read::quality_scale][cor.to][iseq[cor.index]]++;
          
          // P(actual=a|obs=o,a!=o)
          //ntnt_counts[iseq[cor.index]][cor.to]++;
          samples++;
        }
          }
        }
      }
      delete r;
    }
    
    if(++tcount == counts[tchunk] || samples > 200000)
      break;
      }
    }
    reads_in.close();
  }

  regress_probs(ntnt_prob, ntnt_counts);

  output_model(ntnt_prob, ntnt_counts, fqf);
}


////////////////////////////////////////////////////////////
// load_AT_cutoffs
//
// Load AT cutoffs from file
////////////////////////////////////////////////////////////
vector<double> load_AT_cutoffs() {
  vector<double> cutoffs;
  ifstream cut_in(ATcutf);
  string line;
  double cut;
  
  while(getline(cut_in, line)) {
    stringstream ss(stringstream::in | stringstream::out);
    ss << line;
    ss >> cut;
    cutoffs.push_back(cut);
  }

  if(cutoffs.size() != (k+1)) {
    cerr << "Must specify " << (k+1) << " AT cutoffs in " << ATcutf << endl;
    exit(EXIT_FAILURE);
  }

  return cutoffs;
}


////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  parse_command_line(argc, argv);

  // prepare AT and GC counts
  unsigned long long atgc[2] = {0};

  // make trusted kmer data structure
  bithash *trusted = new bithash(k);

  // get good kmers from Hammer
  if (hammerf != NULL) {
    string hammerf_str(hammerf);
    if (hammerf_str.substr(hammerf_str.size()-3) == ".gz") {
    igzstream hammerf_in(hammerf);
    trusted->hammer_file_load(hammerf_in, atgc);
    } else {
    ifstream hammerf_in(hammerf);
    trusted->hammer_file_load(hammerf_in, atgc);
    }   
  }
  
  // get kmer counts
  if(merf != NULL) {
    string merf_str(merf);
    if(ATcutf != NULL) {
      if(merf_str.substr(merf_str.size()-3) == ".gz") {
    igzstream mer_in(merf);
    trusted->tab_file_load(mer_in, load_AT_cutoffs(), atgc);
      } else {
    ifstream mer_in(merf);
    trusted->tab_file_load(mer_in, load_AT_cutoffs(), atgc);
      }
    } else {
      if(merf_str.substr(merf_str.size()-3) == ".gz") {
    igzstream mer_in(merf);
    trusted->tab_file_load(mer_in, cutoff, atgc);
      } else {
    ifstream mer_in(merf);
    trusted->tab_file_load(mer_in, cutoff, atgc);
      }
    }

  // saved bithash
  } else if(bithashf != NULL) {
    if(strcmp(bithashf,"-") == 0) {
      cerr << "Saved bithash cannot be piped in.  Please specify file." << endl;
      exit(EXIT_FAILURE);
    } else
      trusted->binary_file_input(bithashf, atgc);
  }  
  cout << trusted->num_kmers() << " trusted kmers" << endl;

  double prior_prob[4];
  prior_prob[0] = (double)atgc[0] / (double)(atgc[0]+atgc[1]) / 2.0;
  prior_prob[1] = .5 - prior_prob[0];
  prior_prob[2] = prior_prob[1];
  prior_prob[3] = prior_prob[0];
  
  //cout << "AT: " << atgc[0] << " GC: " << atgc[1] << endl;
  cout << "AT% = " << (2*prior_prob[0]) << endl;

  // make list of files
  vector<string> fastqfs;
  vector<int> pairedend_codes;
  parse_fastq(fastqfs, pairedend_codes);

  // process each file
  string fqf;
  bool zip;
  for(int f = 0; f < fastqfs.size(); f++) {
    fqf = fastqfs[f];
    cout << fqf << endl;

    // unzip
    if(fqf.substr(fqf.size()-3) == ".gz") {
      zip = true;
      unzip_fastq(fqf);
    } else
      zip = false;

    // determine quality value scale
    if(Read::quality_scale == -1)
     guess_quality_scale(fqf);

    // split file
    vector<streampos> starts;
    vector<unsigned long long> counts;
    chunkify_fastq(fqf, starts, counts);

    // learn nt->nt transitions
    double ntnt_prob[Read::max_qual][4][4] = {0};
    for(int q = 0; q < Read::max_qual; q++)
      for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      if(i != j)
        ntnt_prob[q][i][j] = 1.0/3.0;

    if(!TESTING)
      learn_errors(fqf, trusted, starts, counts, ntnt_prob, prior_prob);

    // correct
    correct_reads(fqf, pairedend_codes[f], trusted, starts, counts, ntnt_prob, prior_prob);
    
    // combine
    if(pairedend_codes[f] == 0) {
      combine_output(fqf, string("cor"), uncorrected_out);
    }

    // combine paired end
    if(pairedend_codes[f] == 2) {
      if(!zip) {
    combine_output_paired(fastqfs[f-1], fqf, string("cor"), uncorrected_out);
      } else {
    combine_output_paired(fastqfs[f-1].substr(0,fastqfs[f-1].size()-3), fqf, string("cor"), uncorrected_out);
      }
    }

    if(zip)
      zip_fastq(fqf);
  }

  return 0;
}
