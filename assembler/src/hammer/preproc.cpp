/*
 * preproc.cpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */
#include "hammer_config.hpp" 
#include <omp.h>
#include <string>
#include <cstdlib>
#include <vector>
#include <set>
#include <utility>
#include "hammer/defs.hpp"
#include "hammer/kmer_functions.hpp"
#include "common/read/read.hpp"
#include "common/read/ireadstream.hpp"

using std::string;
using std::vector;
using std::set;
using std::pair;
using std::make_pair;

namespace {

char message[100];
const int kStep = 1e5;

struct Options {
  int qvoffset;
  string ifile;
  string ofile;
  size_t nthreads;
  size_t read_batch_size;
  size_t file_number;
  bool valid;
  Options() : nthreads(1), read_batch_size(1e6), file_number(2), valid(true) {}
};


void PrintHelp() {
  printf("Usage: ./preproc qvoffset ifile.fastq ofile.kmer [nthreads]\n");
  printf("Where:\n");
  printf("\tqvoffset\tan offset of fastq quality data\n");
  printf("\tifile.fastq\tan input file with reads in fastq format\n");
  printf("\tofile.kmer\ta filename where k-mer statistics will be outputed\n");
  printf("\tnthreads\ta number of threads (one by default)\n");
}
  
Options ParseOptions(int argc, char * argv[]) {
  Options ret;
  if (argc != 4 && argc != 5) {
    ret.valid =  false;
  } else {
    ret.qvoffset = atoi(argv[1]);  
    ret.valid &= (ret.qvoffset >= 0 && ret.qvoffset <= 255);
    ret.ifile = argv[2];
    ret.ofile = argv[3];
    if (argc == 5) {
      ret.nthreads = atoi(argv[4]);
    }
  }
  return ret;
}

void Log(const string &message) {
  printf("%s", message.c_str());
}
  
void SplitToFiles(const string &ifile, size_t qvoffset, size_t file_number) {
  ireadstream ifs(ifile.c_str(), qvoffset);
  vector<FILE*> files(file_number);
  for (size_t i = 0; i < file_number; ++i) {
    char filename[50];
    sprintf(filename, "%d.kmer.part", (int)i);
    files[i] = fopen(filename, "w");
  }
  size_t read_number = 0;
  while (!ifs.eof()) {
    // reading a batch of reads
    ++read_number;
    if (read_number % kStep == 0) {
      sprintf(message, "Reading read %d.\n", (int)read_number);
      Log(message);
    }
    Read r;      
    ifs >> r; 
    vector<KMer> kmers = GetKMers(r);
    KMer::hash hash_function;
    for (size_t i = 0; i < kmers.size(); ++i) {
      int file_id = hash_function(kmers[i]) % file_number;
      fprintf(files[file_id], "%s\n", kmers[i].str().c_str());
    }
  }    
  for (size_t i = 0; i < file_number; ++i) {
    fclose(files[i]);
  }
  ifs.close();
  Log("Reads wroten to separate files.\n");  
}

void EvalFile(const string &ifile, const string &ofile) {
  sprintf(message, "Processing %s.\n", ifile.c_str());
  Log(message);
  FILE *fin = fopen(ifile.c_str(), "r");  
  FILE *fout = fopen(ofile.c_str(), "w");
  char buffer[K + 1];
  KMerStatMap stat_map;
  while (fscanf(fin, "%s", buffer) != EOF) {   
#pragma message("Warning about uninitialized _M_instance looks like a fake")
    // Next line produces a misterious warning saying that _M_instance is undefined
    // Looks like it is in some way connected to the line
    //       int file_id = hash_function(kmers[i]) % file_number;
    // line in SplitToFiles
    KMer kmer(buffer); 
    ++stat_map[kmer].count;
  }
  for (KMerStatMap::iterator it = stat_map.begin(); it != stat_map.end(); ++it) {
    fprintf(fout, "%s %u\n", it->first.str().c_str(), (unsigned int)it->second.count);
  }
  fclose(fin);
  fclose(fout);
  sprintf(message, "Processed %s. You can find results in %s\n", ifile.c_str(), ofile.c_str());
  Log(message);
}


class KMertPartJoiner {
public:
  KMertPartJoiner (const vector<FILE*> &ifiles) {
    for (size_t i = 0; i < ifiles.size(); ++i) {
      KMerPartParser kpp(ifiles[i]);
      if (!kpp.eof()) {
	kmer_parsers_.insert(kpp);	
      }
    }
  }

  pair<string, int> Next() {
    KMerPartParser kpp (*kmer_parsers_.begin());
    pair<string, int> ret = make_pair(kpp.last_string(), kpp.last_count());
    kmer_parsers_.erase(kpp);
    kpp.Next();
    if (!kpp.eof()) {
      kmer_parsers_.insert(kpp);
    }
    return ret;
  }

  bool IsEmpty() {
    return kmer_parsers_.size() == 0;
  }

private:

  class KMerPartParser {
  public: 
    KMerPartParser(FILE *file) {
      file_ = file;
      eof_ = false;
      Next();
    }
    
    bool operator<(const KMerPartParser &other) const {
      return last_string_ < other.last_string_;
    }

    KMerPartParser(const KMerPartParser &other) {
      file_ = other.file_;
      last_string_ = other.last_string_;
      last_count_ = other.last_count_;
      eof_ = other.eof_;      
    }
    
    void Next() {
      char buf[K + 1];
      eof_ = (fscanf(file_, "%s %d", buf, &last_count_) == EOF);           
      last_string_ = buf;
    }
    
    bool eof() {
      return eof_;
    }
    
    string last_string() {
      return last_string_;
    }
    
    int last_count() {
      return last_count_;
    }

  private:
    string last_string_; 
    int last_count_;
    FILE *file_;
    bool eof_;
  };

  set<KMerPartParser> kmer_parsers_;
};

void MergeAndSort(const vector<FILE*> &ifiles, FILE *ofile) {
  Log("Starting merge.\n");
  KMertPartJoiner joiner(ifiles);
  while (!joiner.IsEmpty()) {
    pair<string, int> kmer_stat = joiner.Next();
    fprintf(ofile, "%s %d\n", kmer_stat.first.c_str(), kmer_stat.second);
  }
}

}

int main(int argc, char * argv[]) {
  Options opts = ParseOptions(argc, argv);
  if (!opts.valid) {
    PrintHelp();
    return 1;
  }

  sprintf(message, "Starting preproc: evaluating %s in %d threads.\n", opts.ifile.c_str(), (int)opts.nthreads);
  Log(message);

  SplitToFiles(opts.ifile, opts.qvoffset, opts.file_number);  
  for (size_t i = 0; i < opts.file_number; ++i) {
    char ifile[50];
    char ofile[50];
    sprintf(ifile, "%d.kmer.part", (int)i);
    sprintf(ofile, "%d.result.part", (int)i);
    EvalFile(ifile, ofile);
  }
  
  vector<FILE*> ifiles;
  for (size_t i = 0; i < opts.file_number; ++i) {
    char ifile_name[50];
    sprintf(ifile_name, "%d.result.part", (int)i);
    FILE *ifile = fopen(ifile_name, "r");
    ifiles.push_back(ifile);
  }

  FILE *ofile = fopen(opts.ofile.c_str(), "w");
  MergeAndSort(ifiles, ofile);
  for (size_t i = 0; i < opts.file_number; ++i) {
    char ifile_name[50];
    sprintf(ifile_name, "%d.result.part", (int)i);
    FILE *ifile = fopen(ifile_name, "r");
    fclose(ifile);
  }
  fclose(ofile);
  sprintf(message, "Preprocessing done. You can find results in %s.\n", opts.ofile.c_str());
  Log(message);
  return 0;
}


